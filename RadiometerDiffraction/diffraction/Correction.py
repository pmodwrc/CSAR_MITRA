import os
import pathlib

import numpy

from diffraction.Diffraction import DiffractionCsar, DiffractionPMO6, DiffractionPMO8Fliana, DiffractionCsarWindowless
from diffraction.LibRadtran import LibRadtranParameter, LibRadtranSpectrum
from diffraction.LookupEvaluation import CorrectionInput
from diffraction.LookupEvaluation import DiffractionLookupTable
from diffraction.Spectra import DiffractionResult, SpectralIntegration


class DiffractionCorrection(DiffractionLookupTable):
    """
        A superclass to retrieve the diffraction correction for different radiometers

        ...

        Attributes
        ----------
        lookup_table_4d : str
            The name of the lookup table, a table  previously generated with the method:
            'generate_new_lookup_table_4d'


        Methods
        -------
        generate_lookup_4d(lookup_table_name, ozone_coulmn_in_du=372.5, model_atmosphere='afglms',
                           zenith_angle_range_limit_davos=True,
                           integrated_water_vapour_in_mm_array=[2, 4, 7, 10, 13, 17, 21, 25, 30, 35],
                           aod_at_1000_nm_array=[0.001, 0.005, 0.007, 0.01, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3])
            Generates a new lookup table for the diffraction correction, with a given name and given interpolation points

        """

    def __init__(self, diffraction_instance, lookup_table='default_lookup_csar.npz'):
        """
            Parameters
            ----------
            diffraction_instance: GenericRadiometer from Diffraction module
                An object of a class that extends GenericRadiometer
            lookup_table : str, optional
                The name of the lookup table, a table  previously generated with the method:
                'load_new_lookup_table_4d' (default is 'default_lookup_csar.npz')
            """
        self.diffraction = diffraction_instance
        self.lookup_table_4d = os.path.join(pathlib.Path(__file__).parent, lookup_table)
        super().__init__(self.lookup_table_4d)

    def generate_lookup_4d(self, lookup_table_name, ozone_coulmn_in_du=315, model_atmosphere='afglms',
                           zenith_angle_range_limit_davos=True,
                           integrated_water_vapour_in_mm_array=[0.99, 2, 4, 7, 10, 13, 17, 21, 25, 30, 35, 42.01],
                           aod_at_1000_nm_array=[0.001, 0.005, 0.007, 0.01, 0.05, 0.07, 0.1, 0.15, 0.2, 0.31],
                           altitude_in_km=1.59):
        """ Generates a new lookup table for the diffraction correction, with a given name and given interpolation points

                Parameters
                ----------
                lookup_table_name: str
                    name of the lookup table to be generated
                ozone_coulmn_in_du: float, optional
                    Amount of ozone [column in dobsun units] to be applied in the model calculation (default is 315)
                model_atmosphere: string, optional
                    Model atmosphere to be applied in the libRadtran model (default is 'afglms')
                zenith_angle_range_limit_davos: bool, optional
                    If true, limits the model calculations to zenith angles larger than 20 degrees (default is True)
                integrated_water_vapour_in_mm_array: list object, optional
                    specifies the interpolation points for integrated water vapour [column of mm] as a list of float
                    (default is [2, 4, 7, 10, 13, 17, 21, 25, 30, 35])
                aod_at_1000_nm_array: list object, optional
                    specifies the interpolation points for the AOD@1000n (Angström Beta) parameter as a list of float
                    (default is [0.001, 0.005, 0.007, 0.01, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3])
                    """
        diffraction = self.diffraction
        diffraction_result = self.__prepare_diffraction_correction(diffraction)
        libradtran_spectrum = LibRadtranSpectrum()
        spectral_integration = SpectralIntegration()
        if zenith_angle_range_limit_davos:
            solar_zenith_angle = numpy.array([20, 30, 40, 50, 55, 60, 65, 70, 75, 80, 85.01], dtype=float)
        else:
            solar_zenith_angle = numpy.array([0, 10, 20, 30, 40, 50, 55, 60, 65, 70, 75, 80, 85.01], dtype=float)

        integrated_water_vapour = numpy.array(integrated_water_vapour_in_mm_array, dtype=float)
        aod_at_1000nm = numpy.array(aod_at_1000_nm_array, dtype=float)
        angstrom = numpy.array([0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.41], dtype=float)
        input_parameter = LibRadtranParameter()
        input_parameter.aerosol_season = 2
        input_parameter.atmosphere = model_atmosphere
        input_parameter.ozone_column_in_du = ozone_coulmn_in_du
        input_parameter.altitude_in_km = altitude_in_km

        window_transmission = diffraction.get_simple_window_transmission_function()
        diffraction_table = numpy.zeros(
            [len(solar_zenith_angle), len(integrated_water_vapour), len(aod_at_1000nm), len(angstrom)])
        for i in range(len(solar_zenith_angle)):
            input_parameter.set_solar_zenith_angle(solar_zenith_angle[i])
            for j in range(len(integrated_water_vapour)):
                input_parameter.set_integrated_water_vapour_in_mm(integrated_water_vapour[j])
                for k in range(len(aod_at_1000nm)):
                    input_parameter.set_aod_at1000nm(aod_at_1000nm[k])
                    for l in range(len(angstrom)):
                        input_parameter.set_angstrom_alpha(angstrom[l])
                        libradtran_filename = "working_file"
                        spectrum_file = libradtran_spectrum.create_spectrum(input_parameter, libradtran_filename)
                        spectrum = numpy.genfromtxt(spectrum_file)
                        diffraction = spectral_integration.integrate_spectrum_with_window_function(spectrum,
                                                                                                   diffraction_result,
                                                                                                   window_transmission)
                        diffraction_table[i, j, k, l] = diffraction
                        print(i, j, k, l)

        numpy.savez(lookup_table_name, diffraction_table, solar_zenith_angle, integrated_water_vapour,
                    aod_at_1000nm,
                    angstrom)

    @staticmethod
    def __prepare_diffraction_correction(diffraction_instance):
        wave_step = 100
        wavelength_vector = numpy.arange(200, 5000 + wave_step, wave_step, dtype=float)
        diffraction_result = DiffractionResult(wavelength_vector,
                                               diffraction_instance.get_diffraction_vector_wolfs_formula(
                                                   wavelength_vector))
        return diffraction_result


class CorrectionCsar(DiffractionCorrection):
    """
        A class to retrieve the diffraction correction for the Cryogenic Solar Absolute Radiometer (CSAR)

        ...
    """

    def __init__(self, lookup_table='default_lookup_csar.npz'):
        self.lookup_table = lookup_table
        self.diffraction = DiffractionCsar()
        super().__init__(self.diffraction, lookup_table=lookup_table)


class CorrectionCsarWindowless(DiffractionCorrection):
    """
        A class to retrieve the diffraction correction for the Cryogenic Solar Absolute Radiometer (CSAR) disregarding
        the entrance window.

        ...
    """

    def __init__(self, lookup_table='lookup_csar_windowless.npz'):
        self.lookup_table = lookup_table
        self.diffraction = DiffractionCsarWindowless()
        super().__init__(self.diffraction, lookup_table=lookup_table)


class CorrectionPMO6(DiffractionCorrection):
    """
        A class to retrieve the diffraction correction for the PMO6 radiometer

        ...
    """

    def __init__(self, lookup_table='lookup_pmo8_example_davos.npz'):
        self.lookup_table = lookup_table
        self.diffraction = DiffractionPMO6()
        super().__init__(self.diffraction, lookup_table=lookup_table)


class CorrectionPMO8Fliana(DiffractionCorrection):
    """
        A class to retrieve the diffraction correction for the PMO8 radiometer

        ...
    """

    def __init__(self, lookup_table='lookup_pmo8_example_davos.npz'):
        self.lookup_table = lookup_table
        self.diffraction = DiffractionPMO8Fliana()
        super().__init__(self.diffraction, lookup_table=lookup_table)


class CorrectionInput(CorrectionInput):
    def __init__(self):
        super().__init__()
