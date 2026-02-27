import os

import numpy
import pathlib


class KuruczSpectrum:
    def __init__(self):
        filename = os.path.join("spectra", "kurudz_1.0nm.dat")
        data = numpy.genfromtxt(filename, skip_header=11)
        self.wavelength_in_mm = data[:, 0] * 1e-6
        self.irradiance_in_mW_per_m2_per_nm = data[:, 1]


class MeanWaterSpectraSza55:
    def __init__(self):
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "sza_55_year_cycle")
        self.spectra = []
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy001m.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy051m.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy101m.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy151m.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy201m.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy251m.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy301m.out")))
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "sza_55_year_cycle", "v2_10000")
        self.spectra_to_10000_nm = []
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy001m.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy051m.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy101m.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy151m.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy201m.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy251m.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy301m.out")))

        self.day_of_year = _day_of_year()
        self.water_column_in_dobson = numpy.array([0.8, 0.8, 1.05, 1.8, 2.25, 1.95, 1.3]) * 1e6


class LowWaterSpectraSza55:
    def __init__(self):
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "sza_55_year_cycle")
        self.spectra = []
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy001l.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy051l.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy101l.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy151l.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy201l.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy251l.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy301l.out")))
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "sza_55_year_cycle", "v2_10000")
        self.spectra_to_10000_nm = []
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy001l.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy051l.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy101l.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy151l.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy201l.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy251l.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy301l.out")))
        self.day_of_year = _day_of_year()
        self.water_column_in_dobson = numpy.array([0.4, 0.4, 0.7, 1.2, 1.7, 1.4, 0.9]) * 1e6


class HighWaterSpectraSza55:
    def __init__(self):
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "sza_55_year_cycle")
        self.spectra = []
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy001h.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy051h.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy101h.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy151h.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy201h.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy251h.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55doy301h.out")))
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "sza_55_year_cycle", "v2_10000")
        self.spectra_to_10000_nm = []
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy001h.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy051h.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy101h.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy151h.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy201h.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy251h.out")))
        self.spectra_to_10000_nm.append(numpy.genfromtxt(os.path.join(directory, "out_iwv_55_doy301h.out")))
        self.day_of_year = _day_of_year()
        self.water_column_in_dobson = numpy.array([1.2, 1.2, 3.0, 2.2, 2.7, 2.2, 1.8]) * 1e6


class SummerStandardAtmosphereSpectra:
    def __init__(self):
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "v2_10000")
        self.spectra = []
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_30v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_31v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_32v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_33v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_34v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_35v2.out")))
        self.solar_zenith_angle = numpy.array([75, 65, 55, 45, 35, 25])


class SummerHighWaterAtmosphereSpectra:
    def __init__(self):
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "v2_10000")
        self.spectra = []
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_40v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_41v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_42v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_43v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_44v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_45v2.out")))
        self.solar_zenith_angle = numpy.array([75, 65, 55, 45, 35, 25])


class WinterStandardAtmosphereSpectra:
    def __init__(self):
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "v2_10000")
        self.spectra = []
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_23v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_22v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_21v2.out")))
        self.solar_zenith_angle = numpy.array([75, 65, 55])


class WinterLowWaterAtmosphereSpectra:
    def __init__(self):
        directory = os.path.join(pathlib.Path(__file__).parent, "spectra", "libratran", "v2_10000")
        self.spectra = []
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_50v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_51v2.out")))
        self.spectra.append(numpy.genfromtxt(os.path.join(directory, "out_spec_52v2.out")))
        self.solar_zenith_angle = numpy.array([75, 65, 55])


class SpectralIntegration:

    def integrate_multiple_spectra(self, spectra, diffraction):
        integrated_diffraction = numpy.zeros(len(spectra))
        for i in range(len(spectra)):
            integrated_diffraction[i] = self.integrate_spectrum(spectra[i], diffraction)
        return integrated_diffraction

    def integrate_multiple_spectra_with_window_function(self, spectra, diffraction, window):
        integrated_diffraction = numpy.zeros(len(spectra))
        for i in range(len(spectra)):
            integrated_diffraction[i] = self.integrate_spectrum_with_window_function(spectra[i], diffraction, window)
        return integrated_diffraction

    @staticmethod
    def integrate_spectrum(spectrum, diffraction):
        wavelength_in_nm = spectrum[:, 0]
        diffraction_regrided = numpy.interp(wavelength_in_nm, diffraction.wavelength_vector_in_nm,
                                            diffraction.diffraction_vector)
        integrated_diffraction = numpy.trapz(spectrum[:, 0], diffraction_regrided * spectrum[:, 1]) \
                                 / numpy.trapz(spectrum[:, 0], spectrum[:, 1])
        return integrated_diffraction

    @staticmethod
    def integrate_spectrum_with_window_function(spectrum, diffraction, window):
        wavelength_in_nm = spectrum[:, 0]
        diffraction_regrided = numpy.interp(wavelength_in_nm, diffraction.wavelength_vector_in_nm,
                                            diffraction.diffraction_vector)
        window_transmission_regrided = numpy.interp(wavelength_in_nm, window.wavelength_vector_in_nm,
                                                    window.transmission_vector)
        integrated_diffraction = numpy.trapz(spectrum[:, 0],
                                             diffraction_regrided * spectrum[:, 1] * window_transmission_regrided) \
                                 / numpy.trapz(spectrum[:, 0], spectrum[:, 1] * window_transmission_regrided)
        return integrated_diffraction


class DiffractionResult:
    def __init__(self, wavelength_in_nm, diffraction):
        self.wavelength_vector_in_nm = wavelength_in_nm
        self.diffraction_vector = diffraction


def _day_of_year():
    return numpy.array([1, 51, 101, 151, 201, 251, 301])
