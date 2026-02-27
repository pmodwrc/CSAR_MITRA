import math
import os
import pathlib
import warnings

import numpy
from scipy.interpolate import interpolate


class DiffractionLookupTable:
    """
        A class to retrieve the diffraction correction and its uncertainty based on a provided lookup-table

        ...

        Attributes
        ----------
        lookup_table_4d : str
            The name of the lookup table, a table  previously generated


        Methods
        -------
        get_diffraction_correction_4d(correction_input)
            Returns the correction factor to correct for the diffraction effect
        get_combined_uncertainty_based_on_absolute_uncertainty(correction_coordinates, delta_abs_sza=0,
                                                               delta_abs_iwv=0, delta_abs_aod=0, delta_abs_angstrom=0)
            Returns the uncertainty of the diffraction correction factor based on the absolute uncertainty of the input
            parameters
        get_combined_uncertainty_based_on_relative_uncertainty(correction_coordinates, delta_rel_sza=0,
                                                               delta_rel_iwv=0,
                                                               delta_rel_aod=0, delta_rel_angstrom=0)
            Returns the uncertainty of the diffraction correction factor based on the relative uncertainty of the input
            parameters
        load_new_lookup_table_4d(table_file_name)
            Loads a new lookup table file
        """

    def __init__(self, lookup_table):
        """
            Parameters
            ----------
            lookup_table : str
                The name of the lookup table, a table  previously generated
            """
        self.lookup_table_4d = os.path.join(pathlib.Path(__file__).parent, lookup_table)
        try:
            self.__set_up_lookup_4d()
        except:
            warnings.warn('Could not load requested lookup table')
            print('Could not load lookup table')

    def load_new_lookup_table_4d(self, table_file_name):
        """ Loads a new lookup table file

                Parameters
                ----------
                table_file_name: str
                    The name of the lookup table file to load
                """
        self.lookup_table_4d = table_file_name
        self.__set_up_lookup_4d()

    def get_diffraction_correction_4d(self, correction_input):
        """ Returns the correction factor to correct for the diffraction effect

            Parameters
            ----------
            correction_input : CorrectionInput
                CorrectionInput Class, containing the input parameters Solar Zenith Angle, Integrated Water Vapour,
                AOD Alpha and Beta.
            """
        return self.__evaluate_diffraction_table_4d(correction_input)

    def get_combined_uncertainty_based_on_absolute_uncertainty(self, correction_coordinates, delta_abs_sza=0.0,
                                                               delta_abs_iwv=0.0, delta_abs_aod=0.0,
                                                               delta_abs_angstrom_alpha=0.0):
        """ Returns the uncertainty of the diffraction correction factor based on the absolute uncertainty of the input
            parameters

            Parameters
            ----------
            correction_coordinates : CorrectionInput
                CorrectionInput Class, containing the input parameters Solar Zenith Angle, Integrated Water Vapour,
                AOD Alpha and Beta.
            delta_abs_sza : float, optional
                Standard uncertainty of the input parameter 'Solar Zenith Angle' in degrees (default is 0)
            delta_abs_iwv: float, optional
                Standard uncertainty of the input parameter 'Integrated Water Vapor' in mm (default is 0)
            delta_abs_aod: float, optional
                Standard uncertainty of the input parameter 'AOD@1000nm or Angström Beta' (default is 0)
            delta_abs_angstrom_alpha: float, optional
                Standard uncertainty of the input parameter 'AOD Angström Exponent Alpha' (default is 0)
            """
        return ((self.__get_partial_derivative_sza(correction_coordinates,
                                                   delta_sza=delta_abs_sza) * delta_abs_sza) ** 2
                + (self.__get_partial_derivative_iwv(correction_coordinates,
                                                     delta_iwv=delta_abs_iwv) * delta_abs_iwv) ** 2
                + (self.__get_partial_derivative_aod(correction_coordinates,
                                                     delta_aod=delta_abs_aod) * delta_abs_aod) ** 2
                + (self.__get_partial_derivative_angstrom_alpha(correction_coordinates,
                                                                delta_angstrom=delta_abs_angstrom_alpha)
                   * delta_abs_angstrom_alpha) ** 2) ** 0.5

    def get_combined_uncertainty_based_on_relative_uncertainty(self, correction_coordinates, delta_rel_sza=0.0,
                                                               delta_rel_iwv=0.0,
                                                               delta_rel_aod=0.0, delta_rel_angstrom_alpha=0.0):
        """Returns the uncertainty of the diffraction correction factor based on the absolute uncertainty of the input
            parameters

            Parameters
            ----------
            correction_coordinates : CorrectionInput
                CorrectionInput Class, containing the input parameters Solar Zenith Angle, Integrated Water Vapour,
                AOD Alpha and Beta.
            delta_rel_sza : float, optional
                Relative standard uncertainty of the input parameter 'Solar Zenith Angle' in degrees (default is 0)
            delta_rel_iwv: float, optional
                Relative standard uncertainty of the input parameter 'Integrated Water Vapor' in mm (default is 0)
            delta_rel_aod: float, optional
                Relative standard uncertainty of the input parameter 'AOD@1000nm or Angström Beta' (default is 0)
            delta_rel_angstrom_alpha: float, optional
                Relative standard uncertainty of the input parameter 'AOD Angström Exponent Alpha' (default is 0)
                        """
        delta_sza = correction_coordinates.solar_zenith_angle * delta_rel_sza
        delta_iwv = correction_coordinates.integrated_water_vapour_in_mm * delta_rel_iwv
        delta_aod = correction_coordinates.aerosol_optical_depth_at_1000nm * delta_rel_aod
        delta_angstrom = correction_coordinates.angstrom_alpha * delta_rel_angstrom_alpha
        return ((self.__get_partial_derivative_sza(correction_coordinates, delta_sza=delta_sza) * delta_sza) ** 2
                + (self.__get_partial_derivative_iwv(correction_coordinates, delta_iwv=delta_iwv) * delta_iwv) ** 2
                + (self.__get_partial_derivative_aod(correction_coordinates, delta_aod=delta_aod) * delta_aod) ** 2
                + (self.__get_partial_derivative_angstrom_alpha(correction_coordinates, delta_angstrom=delta_angstrom)
                   * delta_angstrom) ** 2) ** 0.5

    def __evaluate_diffraction_table_4d(self, correction_coordinates, delta_sza=0, delta_iwv=0, delta_aod=0,
                                        delta_angstrom=0):
        if (isinstance(correction_coordinates.solar_zenith_angle, numpy.ndarray) and
                isinstance(correction_coordinates.integrated_water_vapour_in_mm, numpy.ndarray) and
                isinstance(correction_coordinates.aerosol_optical_depth_at_1000nm, numpy.ndarray) and
                isinstance(correction_coordinates.angstrom_alpha, numpy.ndarray)):
            x = numpy.array([correction_coordinates.solar_zenith_angle + delta_sza,
                             correction_coordinates.integrated_water_vapour_in_mm + delta_iwv,
                             correction_coordinates.aerosol_optical_depth_at_1000nm + delta_aod,
                             correction_coordinates.angstrom_alpha + delta_angstrom])
            return 1 / interpolate.interpn([self.sza_4d, self.water_vapour_4d, self.aod_4d, self.angstrom_4d],
                                           self.diffraction_4d,
                                           x.T,
                                           bounds_error=False, fill_value=math.nan)
        else:
            return 1 / interpolate.interpn([self.sza_4d, self.water_vapour_4d, self.aod_4d, self.angstrom_4d],
                                           self.diffraction_4d,
                                           [correction_coordinates.solar_zenith_angle + delta_sza,
                                            correction_coordinates.integrated_water_vapour_in_mm + delta_iwv,
                                            correction_coordinates.aerosol_optical_depth_at_1000nm + delta_aod,
                                            correction_coordinates.angstrom_alpha + delta_angstrom],
                                           bounds_error=False, fill_value=math.nan)

    def __get_partial_derivative_sza(self, correction_coordinates, delta_sza):
        small_delta = 0.01
        derivative_calculated_with_order_of_magnitude = self.__partial_derivative_sza(correction_coordinates, delta_sza)
        derivative_calculated_with_small_delta = self.__partial_derivative_sza(correction_coordinates, small_delta)
        return self.__return_alternative_value_if_first_is_not_finite(derivative_calculated_with_order_of_magnitude,
                                                                      derivative_calculated_with_small_delta)

    def __get_partial_derivative_iwv(self, correction_coordinates, delta_iwv):
        small_delta = 0.01
        derivative_calculated_with_order_of_magnitude = self.__partial_derivative_iwv(correction_coordinates,
                                                                                      delta_iwv)
        derivative_calculated_with_small_delta = self.__partial_derivative_iwv(correction_coordinates, small_delta)
        return self.__return_alternative_value_if_first_is_not_finite(derivative_calculated_with_order_of_magnitude,
                                                                      derivative_calculated_with_small_delta)

    def __get_partial_derivative_aod(self, correction_coordinates, delta_aod):
        small_delta = 0.0005
        derivative_calculated_with_order_of_magnitude = self.__partial_derivative_aod(correction_coordinates,
                                                                                      delta_aod)
        derivative_calculated_with_small_delta = self.__partial_derivative_aod(correction_coordinates,
                                                                               small_delta)
        return self.__return_alternative_value_if_first_is_not_finite(derivative_calculated_with_order_of_magnitude,
                                                                      derivative_calculated_with_small_delta)

    def __get_partial_derivative_angstrom_alpha(self, correction_coordinates, delta_angstrom):
        small_delta = 0.01
        derivative_calculated_with_order_of_magnitude = self.__partial_derivative_angstrom_alpha(correction_coordinates,
                                                                                                 delta_angstrom)
        derivative_calculated_with_small_delta = self.__partial_derivative_angstrom_alpha(correction_coordinates,
                                                                                          small_delta)
        return self.__return_alternative_value_if_first_is_not_finite(derivative_calculated_with_order_of_magnitude,
                                                                      derivative_calculated_with_small_delta)

    def __partial_derivative_sza(self, correction_coordinates, delta_sza):
        if numpy.all(delta_sza):
            return (self.__evaluate_diffraction_table_4d(correction_coordinates, delta_sza=delta_sza) -
                    self.__evaluate_diffraction_table_4d(correction_coordinates, delta_sza=-delta_sza)) / (
                           2 * delta_sza)
        else:
            return delta_sza * 0

    def __partial_derivative_iwv(self, correction_coordinates, delta_iwv):
        if numpy.all(delta_iwv):
            return (self.__evaluate_diffraction_table_4d(correction_coordinates, delta_iwv=delta_iwv) -
                    self.__evaluate_diffraction_table_4d(correction_coordinates, delta_iwv=-delta_iwv)) / (
                           2 * delta_iwv)
        else:
            return delta_iwv * 0

    def __partial_derivative_aod(self, correction_coordinates, delta_aod):
        if numpy.all(delta_aod):
            return (self.__evaluate_diffraction_table_4d(correction_coordinates, delta_aod=delta_aod) -
                    self.__evaluate_diffraction_table_4d(correction_coordinates, delta_aod=-delta_aod)) / (
                           2 * delta_aod)
        else:
            return delta_aod * 0

    def __partial_derivative_angstrom_alpha(self, correction_coordinates, delta_angstrom):
        if numpy.all(delta_angstrom):
            return (self.__evaluate_diffraction_table_4d(correction_coordinates, delta_angstrom=delta_angstrom, ) -
                    self.__evaluate_diffraction_table_4d(correction_coordinates, delta_angstrom=-delta_angstrom, )) / \
                   (2 * delta_angstrom)
        else:
            return delta_angstrom * 0

    @staticmethod
    def __return_alternative_value_if_first_is_not_finite(value1, alternative):
        if numpy.isfinite(value1).all():
            return value1
        else:
            return alternative


    def __set_up_lookup_4d(self):
        npzfile = numpy.load(self.lookup_table_4d)
        self.diffraction_4d = npzfile['arr_0']
        self.sza_4d = npzfile['arr_1']
        self.water_vapour_4d = npzfile['arr_2']
        self.aod_4d = npzfile['arr_3']
        self.angstrom_4d = npzfile['arr_4']


class CorrectionInput:
    """
        Input (data) Class for methods within the CorrectionCsar and similar classes, includes Solar Zenith Angle,
        Integrated Water Vapour and AOD parameter Alpha and Beta

        The described attributes can also be numpy arrays of the given type. This allows to call the methods from
        the CorrectionCsar class directly for multiple outputs. In this case all of the attribute arrays must have
        the same length.

        Attributes
        ----------
        solar_zenith_angle: float
            The solar zenith angle in degrees
        integrated_water_vapour_in_mm: float
            The integrated water vapour column in mm
        aerosol_optical_depth_at_1000nm: float
            The aerosol optical depth (AOD) Beta parameter according to the Angström formula tau = beta * lambda ^ - alpha
        angstrom_alpha: float
                The aerosol optical depth (AOD) Alpha parameter according to the Angström formula tau = beta * lambda ^ - alpha

        Methods
        -------
        set_solar_zenith_angle_in_degrees(self, solar_zenith_angle_in_degrees)

        set_integrated_water_vapour_in_mm(self, integrated_water_vapour_in_mm)

        set_aerosol_optical_depth_at_1000nm(self, aod_beta)

        set_angstrom_alpha(self, angstrom_alpha)
        """

    def __init__(self):
        self.solar_zenith_angle = 45
        self.integrated_water_vapour_in_mm = 10
        self.aerosol_optical_depth_at_1000nm = 0.005
        self.angstrom_alpha = 1.7

    def set_solar_zenith_angle_in_degrees(self, solar_zenith_angle_in_degrees):
        """ Sets the solar zenith angle parameter in degrees (0-90)

                Parameters
                ----------
                solar_zenith_angle_in_degrees: float
                """
        self.solar_zenith_angle = solar_zenith_angle_in_degrees

    def set_integrated_water_vapour_in_mm(self, integrated_water_vapour_in_mm):
        """ Sets the integrated water vapour parameter

                 Parameters
                 ----------
                 integrated_water_vapour_in_mm: float
                 """
        self.integrated_water_vapour_in_mm = integrated_water_vapour_in_mm

    def set_aerosol_optical_depth_at_1000nm(self, aod_beta):
        """ Sets the AOD Beta parameter (AOD at 1000 nm) according to Angstroms parametrisation.

                 Parameters
                 ----------
                 aod_beta: float
                 """
        self.aerosol_optical_depth_at_1000nm = aod_beta

    def set_angstrom_alpha(self, angstrom_alpha):
        """ Sets the AOD Alpha parameter (Angstrom exponent) according to Angstroms parametrisation.

                 Parameters
                 ----------
                 angstrom_alpha: float
                 """
        self.angstrom_alpha = angstrom_alpha
