import math

import numpy.polynomial
import scipy.integrate


class GenericRadiometer:
    """
        A superclass to calculate the wavelength dependent Diffraction Effect for Solar Radiometers with
        DARA/CSAR style geometry, as well as PMO6 geometry
        ...

        Attributes
        ----------
        geometry : diffraction.Geometry
            The geometrical properties of the instrument


        Methods
        -------
        get_diffraction_wolfs_formula(wavelength_in_nm)
            Returns the relative diffraction effect, using Wolf's formula
        get_diffraction_revised_formula(wavelength_in_nm)
            Returns the relative diffraction effect, using Shirley's revised formula
            (not available for inverted [PMO6-like] geometries)
        get_diffraction_vector_wolfs_formula(wavelength_vector_in_nm)
            Returns the relative diffraction effects, using Wolf's formula, supporting array input
        get_diffraction_vector_revised_formula(wavelength_vector_in_nm):
            Returns the relative diffraction effects, using Shirley's revised formula, supporting array input
            (not available for inverted [PMO6-like] geometries)
            """

    def __init__(self, geometry):
        self.geometry = geometry
        if geometry.inverted:
            self.wolf = DiffractionWolfInverted(self.geometry)
        else:
            self.wolf = DiffractionWolf(self.geometry)

            self.revised = DiffractionRevised(self.geometry)

    def get_diffraction_wolfs_formula(self, wavelength_in_nm):
        """ Returns the relative effect of diffraction effect, using wolfs formula.
        Or in other words the actual flux (considering diffraction effects) divided by the flux expected from pure
        geometry.
        (E. Wolf. Light Distribution Near Focus in an Error Free Diffraction Image. Royal Society of London Proceedings
        Series A, 204:533-548, January 1951)

                    Parameters
                    ----------
                    wavelength_in_nm : float
                        the wavelength in nm
        """
        return self.wolf.get_diffraction(wavelength_in_nm)

    def get_diffraction_revised_formula(self, wavelength_in_nm):
        """ Returns the relative effect of diffraction effect, using the revised formula from E. Shirley. Or in other
        words the actual flux (considering diffraction effects) divided by the flux expected from pure geometry.
        (E.L. Shirley. Revised Formulas for Diffraction Effects with Point and Extended Sources. Applied Optics,
        37:6581:6590, October 1998)

                    Parameters
                    ----------
                    wavelength_in_nm : float
                        the wavelength in nm
        """
        return self.revised.get_diffraction(wavelength_in_nm)

    def get_diffraction_vector_wolfs_formula(self, wavelength_vector_in_nm):
        """ Returns an array of floats, representing the relative effect of diffraction effect, using wolfs formula.
        Or in other words the actual flux (considering diffraction effects) divided by the flux expected from pure
        geometry.
        (E. Wolf. Light Distribution Near Focus in an Error Free Diffraction Image. Royal Society of London Proceedings
        Series A, 204:533-548, January 1951)

                    Parameters
                    ----------
                    wavelength_vector_in_nm : array
                        an array of floats the wavelength in nm
        """
        integrated_diffraction = numpy.zeros(len(wavelength_vector_in_nm))
        for i in range(len(wavelength_vector_in_nm)):
            integrated_diffraction[i] = self.get_diffraction_wolfs_formula(wavelength_vector_in_nm[i])
        return integrated_diffraction

    def get_diffraction_vector_revised_formula(self, wavelength_vector_in_nm):
        """ Returns the relative effect of diffraction effect, using the revised formula from E. Shirley. Or in other
        words the actual flux (considering diffraction effects) divided by the flux expected from pure geometry.
        (not available for inverted [PMO6-like] geometry.)
        (E.L. Shirley. Revised Formulas for Diffraction Effects with Point and Extended Sources. Applied Optics,
        37:6581:6590, October 1998)

                    Parameters
                    ----------
                    wavelength_vector_in_nm : array
                        an array of floats the wavelength in nm
        """
        integrated_diffraction = numpy.zeros(len(wavelength_vector_in_nm))
        for i in range(len(wavelength_vector_in_nm)):
            integrated_diffraction[i] = self.get_diffraction_revised_formula(wavelength_vector_in_nm[i])
        return integrated_diffraction


class DiffractionCsar(GenericRadiometer):
    """
            A class to calculate the wavelength dependent Diffraction Effect for the CSAR Aperture Geometry
            ...

            Attributes
            ----------
            geometry : diffraction.Geometry
                The geometrical properties of the CSAR instrument


            Methods
            -------
            get_simple_window_transmission_function():
                Returns the simple window transmission function for the CSAR Suprasil entrance window, an array
                containing wavelength in nm and the corresponding window transmission.
            get_advanced_window_transmission_function():
                Returns the more detailed window transmission function for the CSAR Suprasil entrance window, an array
                containing wavelength in nm and the corresponding window transmission.
            """

    def __init__(self):
        self.geometry = Geometry()
        self.geometry.radiusPrecisionAperture_in_mm = 2.5
        self.geometry.distanceBetweenApertures_in_mm = 103.98
        self.geometry.radiusDefiningAperture_in_mm = 5
        self.geometry.inverted = False
        self.__window_transmission_function_wavelength_in_nm = numpy.arange(1, 10001, 1, dtype=float)
        self.__window_transmission_function_relative_transmission = self.__simple_window_transmission_function()
        super().__init__(self.geometry)

    def get_simple_window_transmission_function(self):
        return WindowTransmission(self.__window_transmission_function_wavelength_in_nm,
                                  self.__simple_window_transmission_function())

    def get_advanced_window_transmission_function(self):
        return WindowTransmission(self.__window_transmission_function_wavelength_in_nm,
                                  self.__advanced_window_transmission_function())

    @staticmethod
    def __simple_window_transmission_function():
        transmission_function = numpy.zeros(10000, dtype=float)
        transmission_function[240:2760] = 0.95
        transmission_function[2760:4000] = 0.9
        transmission_function[3500:4400] = 0.9 - numpy.arange(0, 900, 1) * 0.9 / 900
        return transmission_function

    @staticmethod
    def __advanced_window_transmission_function():
        transmission_function = numpy.zeros(10000, dtype=float)
        transmission_function[240:2630] = 0.955
        transmission_function[2630:3000] = 0.955 - numpy.arange(0, 370, 1) * ((0.955 - 0.889) / 370)
        transmission_function[3000:3400] = 0.889
        transmission_function[3400:3770] = 0.889 - numpy.arange(0, 370, 1) * ((0.889 - 0.415) / 370)
        transmission_function[3770:3950] = 0.415
        transmission_function[3950:4395] = 0.415 - numpy.arange(0, 445, 1) * ((0.414 - 0.025) / 445)
        transmission_function[4395:4700] = 0.025
        return transmission_function


class DiffractionCsarWindowless(GenericRadiometer):
    def __init__(self):
        """All numbers in mm"""
        self.geometry = Geometry()
        self.geometry.radiusPrecisionAperture_in_mm = 2.5
        self.geometry.distanceBetweenApertures_in_mm = 103.98
        self.geometry.radiusDefiningAperture_in_mm = 5
        self.geometry.inverted = False
        self.__window_transmission_function_wavelength_in_nm = numpy.arange(1, 10001, 1, dtype=float)
        self.__window_transmission_function_relative_transmission = numpy.ones(10000)
        super().__init__(self.geometry)

    def get_simple_window_transmission_function(self):
        return WindowTransmission(self.__window_transmission_function_wavelength_in_nm,
                                  self.__window_transmission_function_relative_transmission)

    def get_advanced_window_transmission_function(self):
        return self.get_simple_window_transmission_function()


class DiffractionDara(GenericRadiometer):
    def __init__(self):
        """All numbers in mm"""
        self.geometry = Geometry()
        self.geometry.radiusPrecisionAperture_in_mm = 2.5
        self.geometry.distanceBetweenApertures_in_mm = 54.1
        self.geometry.radiusDefiningAperture_in_mm = 6.9 / 2
        self.geometry.inverted = False
        super().__init__(self.geometry)


class DiffractionPMO6(GenericRadiometer):
    def __init__(self):
        """All numbers in mm"""
        self.geometry = Geometry()
        self.geometry.radiusPrecisionAperture_in_mm = 8.5 / 2
        self.geometry.distanceBetweenApertures_in_mm = 95.4
        self.geometry.radiusDefiningAperture_in_mm = 2.5
        self.geometry.inverted = True
        self.__window_transmission_function_wavelength_in_nm = numpy.arange(1, 10001, 1, dtype=float)
        self.__window_transmission_function_relative_transmission = numpy.ones(10000)
        super().__init__(self.geometry)

    def get_simple_window_transmission_function(self):
        return WindowTransmission(self.__window_transmission_function_wavelength_in_nm,
                                  self.__window_transmission_function_relative_transmission)

    def get_advanced_window_transmission_function(self):
        return self.get_simple_window_transmission_function()


class DiffractionPMO8Fliana(GenericRadiometer):
    def __init__(self):
        """All numbers in mm"""
        self.geometry = Geometry()
        self.geometry.radiusPrecisionAperture_in_mm = 8.3 / 2
        self.geometry.distanceBetweenApertures_in_mm = 95
        self.geometry.radiusDefiningAperture_in_mm = 2.5
        self.geometry.inverted = True
        self.__window_transmission_function_wavelength_in_nm = numpy.arange(1, 10001, 1, dtype=float)
        self.__window_transmission_function_relative_transmission = numpy.ones(10000)
        super().__init__(self.geometry)

    def get_simple_window_transmission_function(self):
        return WindowTransmission(self.__window_transmission_function_wavelength_in_nm,
                                  self.__window_transmission_function_relative_transmission)

    def get_advanced_window_transmission_function(self):
        return self.get_simple_window_transmission_function()


class DiffractionWolf:
    def __init__(self, geometry):
        sun = Sun()
        self.radiusSun = sun.radius_in_mm
        self.distanceSun = sun.distance_to_Earth_in_mm
        self.radiusPrecisionAperture = geometry.radiusPrecisionAperture_in_mm
        self.distanceBetweenApertures = geometry.distanceBetweenApertures_in_mm
        self.radiusDefiningAperture = geometry.radiusDefiningAperture_in_mm

    def get_diffraction(self, wavelength_in_nm):
        """Calculates the Diffraction Effect according to Wolfs Formula"""
        wavelength_in_mm = wavelength_in_nm * 1e-6
        q = 2 * math.pi / wavelength_in_mm
        u = q * self.radiusPrecisionAperture ** 2 * (1 / self.distanceSun + 1 / self.distanceBetweenApertures)
        vs = q * (self.radiusSun / self.distanceSun) * self.radiusPrecisionAperture
        vd = q * (self.radiusDefiningAperture / self.distanceBetweenApertures) * self.radiusPrecisionAperture
        vm = max(vs, vd)
        sig = min(vs, vd) / max(vs, vd)
        l_geom = u ** 2 / vm ** 2 / (self.distanceBetweenApertures + self.distanceSun) ** 2
        pref = 4 * math.pi * self.radiusPrecisionAperture ** 4 / (
                self.distanceSun ** 2 * self.distanceBetweenApertures ** 2 * (wavelength_in_mm * vm) ** 2)
        f = lambda x: (((1 - x ** 2) * ((2 + sig * x) ** 2 - sig ** 2)) ** (1 / 2) * self._l_wolf(u, vm * (
                1 + sig * x))) / (1 + sig * x)
        phi = scipy.integrate.quad(f, -1, 1)
        return (phi[0] * pref) / l_geom

    def _l_wolf(self, u, v):
        w = u / v
        ll = 2 * self._sig_k(0, w) / (math.pi * v) \
             - (self._sig_k(0, w) * math.cos(2 * v)) / (math.pi * v ** 2)
        - (16 * self._sig_k(4, w) + 32 * self._sig_k(3, w) + 8 * self._sig_k(2, w) - 8 * self._sig_k(1, w) + 9
           * self._sig_k(0, w)) / (12 * v ** 3 * math.pi)
        + ((8 * self._sig_k(2, w) + 8 * self._sig_k(1, w) - self._sig_k(0, w)) / (4 * math.pi * v ** 3)) \
        * math.sin(2 * v) + ((64 * self._sig_k(4, w) + 128 * self._sig_k(3, w) - 16 * self._sig_k(2, w) - 80
                              * self._sig_k(1, w) + 9 * self._sig_k(0, w))
                             / (32 * math.pi * v ** 4)) * math.cos(2 * v)
        return 1 - ll

    def _sig_k(self, k, w):
        return self._w(k, w ** 2) / (1 - w ** 2) ** (k + 1)

    def _w(self, m, x):
        switcher = {
            0: float(1),
            1: float(x),
            2: numpy.polyval([1, 1, 0], x),
            3: numpy.polyval([1, 4, 1, 0], x),
            4: numpy.polyval([1, 11, 11, 1, 0], x),
            5: numpy.polyval([1, 26, 66, 26, 1, 0], x),
        }
        return switcher.get(m, float("nan"))


class DiffractionWolfInverted:
    def __init__(self, geometry):
        sun = Sun()
        self.radiusSun = sun.radius_in_mm
        self.distanceSun = sun.distance_to_Earth_in_mm
        self.radiusPrecisionAperture = geometry.radiusPrecisionAperture_in_mm
        self.distanceBetweenApertures = geometry.distanceBetweenApertures_in_mm
        self.radiusDefiningAperture = geometry.radiusDefiningAperture_in_mm

    def get_diffraction(self, wavelength_in_nm):
        """Calculates the Diffraction Effect according to Wolfs Formula Inverted"""
        wavelength_in_mm = wavelength_in_nm * 1e-6
        q = 2 * math.pi / wavelength_in_mm
        u = q * self.radiusPrecisionAperture ** 2 * (1 / self.distanceSun + 1 / self.distanceBetweenApertures)
        vs = q * (self.radiusSun / self.distanceSun) * self.radiusPrecisionAperture
        vd = q * (self.radiusDefiningAperture / self.distanceBetweenApertures) * self.radiusPrecisionAperture
        vm = max(vs, vd)
        sig = min(vs, vd) / max(vs, vd)
        l_geom = u ** 2 / vm ** 2 / (self.distanceBetweenApertures + self.distanceSun) ** 2
        pref = 4 * math.pi * self.radiusPrecisionAperture ** 4 / (
                self.distanceSun ** 2 * self.distanceBetweenApertures ** 2 * (wavelength_in_mm * vm) ** 2)
        f = lambda x: (((1 - x ** 2) * ((2 + sig * x) ** 2 - sig ** 2)) ** (1 / 2) * self._l_wolf_inverted(u, vm * (
                1 + sig * x))) / (1 + sig * x)
        phi = scipy.integrate.quad(f, -1, 1)
        return (phi[0] * pref) / l_geom

    def _l_wolf_inverted(self, u, v):
        w = v / u
        ll = 2 * self._sig_k(0, w) / (math.pi * v) \
             - (self._sig_k(0, w) * math.cos(2 * v)) / (math.pi * v ** 2)
        - (16 * self._sig_k(4, w) + 32 * self._sig_k(3, w) + 8 * self._sig_k(2, w) - 8 * self._sig_k(1, w) + 9
           * self._sig_k(0, w)) / (12 * v ** 3 * math.pi)
        + ((8 * self._sig_k(2, w) + 8 * self._sig_k(1, w) - self._sig_k(0, w)) / (4 * math.pi * v ** 3)) \
        * math.sin(2 * v) + ((64 * self._sig_k(4, w) + 128 * self._sig_k(3, w) - 16 * self._sig_k(2, w) - 80
                              * self._sig_k(1, w) + 9 * self._sig_k(0, w))
                             / (32 * math.pi * v ** 4)) * math.cos(2 * v)
        return 1 + ll - 4 / u * (self._Y1(u, v) * math.cos(1 / 2 * (u + v ** 2 / u)) + self._Y2(u, v) * math.sin(
            1 / 2 * (u + v ** 2 / u)))

    def _sig_k(self, k, w):
        return self._w(k, w ** 2) / (1 - w ** 2) ** (k + 1)

    def _w(self, m, x):
        switcher = {
            0: float(1),
            1: float(x),
            2: numpy.polyval([1, 1, 0], x),
            3: numpy.polyval([1, 4, 1, 0], x),
            4: numpy.polyval([1, 11, 11, 1, 0], x),
            5: numpy.polyval([1, 26, 66, 26, 1, 0], x),
        }
        return switcher.get(m, float("nan"))

    def _Y1(self, u, v):
        w = v / u
        return v ** 2 / 2 / u * (2 * math.pi / v) ** 5 * (
                (2 * self._sig_k(0, w) + 4 * self._sig_k(1, w)) / v * math.sin(v - math.pi / 4)
                + ((3 * self._sig_k(0, w) + 22 * self._sig_k(1, w) + 48 * self._sig_k(2, w) + 32 * self._sig_k(3,
                                                                                                               w)) / 4 / v ** 2)
                * math.cos(v - math.pi / 4) + (
                        15 * self._sig_k(0, w) + 62 * self._sig_k(1, w) - 160 * self._sig_k(2, w)
                        - 960 * self._sig_k(3, w) - 1280 * self._sig_k(4, w) - 512 * self._sig_k(5,
                                                                                                 w)) / 64 / v ** 3 * math.sin(
            v - math.pi / 4))

    def _Y2(self, u, v):
        w = v / u
        return v / 2 / u * (2 * math.pi / v) ** 5 * ((4 * self._sig_k(1, w) / v) * math.sin(v + math.pi / 4)
                                                     + ((-self._sig_k(1, w) + 16 * self._sig_k(3, w)) / (
                        2 * v ** 2)) * math.cos(v + math.pi / 4)
                                                     + ((-9 * self._sig_k(1, w) + 160 * self._sig_k(3,
                                                                                                    w) - 256 * self._sig_k(
                    5, w)) / (32 * v ** 3))
                                                     * math.sin(v - math.pi / 4))


class DiffractionRevised:
    def __init__(self, geometry):
        sun = Sun()
        self.radiusSun = sun.radius_in_mm
        self.distanceSun = sun.distance_to_Earth_in_mm
        self.radiusPrecisionAperture = geometry.radiusPrecisionAperture_in_mm
        self.distanceBetweenApertures = geometry.distanceBetweenApertures_in_mm
        self.radiusDefiningAperture = geometry.radiusDefiningAperture_in_mm

    def get_diffraction(self, wavelength_in_nm):
        """Calculates the diffraction effect according to Steel et al, revised formula by E. Shirley"""
        wavelength_in_mm = wavelength_in_nm * 1e-6
        q = 2 * math.pi / wavelength_in_mm
        u = q * self.radiusPrecisionAperture ** 2 * (1 / self.distanceSun + 1 / self.distanceBetweenApertures)
        vs = q * (self.radiusSun / self.distanceSun) * self.radiusPrecisionAperture
        vd = q * (self.radiusDefiningAperture / self.distanceBetweenApertures) * self.radiusPrecisionAperture
        vm = max(vs, vd)
        sig = min(vs, vd) / max(vs, vd)
        f_revised = lambda sig_x: (1 - (2 * vm * (1 + sig_x)) / (math.pi * (vm ** 2 * (1 + sig_x) ** 2 - u ** 2)) + (
            math.cos(2 * vm * (1 + sig_x))) / (math.pi * (vm ** 2 * (1 + sig_x) ** 2 - u ** 2)))

        i_revised = lambda x: (((1 - x ** 2) * ((2 + sig * x) ** 2 - sig ** 2)) ** .5 * f_revised(sig * x)) / (
                1 + sig * x)

        integral = scipy.integrate.quad(i_revised, -1, 1)
        return integral[0] / math.pi


class Geometry:
    """
           A class describing the geometry of thr aperture set-up of a radiometer
           ...

           Attributes
           ----------
           radiusPrecisionAperture_in_mm : float
                The radius of the precision aperture in mm
           distanceBetweenApertures_in_mm : float
                The distance between the apertures in mm
           radiusDefiningAperture_in_mm : float
                The radius of the view limiting aperture in mm
           inverted : boolean
                Describes whether the set-up is a DARA/CSAR style geometry with the precision aperture in front
                (inverted = False)
                or a PMO6/PMO8 style geometry with the view limiting aperture in front (inverted = True)
    """

    def __init__(self):
        self.radiusPrecisionAperture_in_mm = 0
        self.distanceBetweenApertures_in_mm = 0
        self.radiusDefiningAperture_in_mm = 0
        self.inverted = False


class Sun:
    def __init__(self):
        self.radius_in_mm = 6.75e11
        self.distance_to_Earth_in_mm = 1.5e14


class WindowTransmission:
    def __init__(self, wavelength_vector_in_nm, transmission_vector):
        self.wavelength_vector_in_nm = wavelength_vector_in_nm
        self.transmission_vector = transmission_vector
