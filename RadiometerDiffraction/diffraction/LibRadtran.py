import os
import sys


class LibRadtranSpectrum:
    def __init__(self):
        self.work_directory = os.path.join("spectra", "libradtran_work")
        os.makedirs(self.work_directory, exist_ok=True)
        self.docker_skript_name = "runLibRadtran.sh"
        self.docker_skript_path = os.path.join(self.work_directory, self.docker_skript_name)
        self.rte_solver = 'disort'

    def create_spectrum(self, libradtran_parameter, file_name):
        input_file_name = "{name}.inp".format(name=file_name)
        output_file_name = "{name}.out".format(name=file_name)
        input_file_path = os.path.join(self.work_directory, input_file_name)
        output_file_path = os.path.join(self.work_directory, output_file_name)
        fid = open(input_file_path, "w")
        default_string = "output_user lambda  edir  \n" \
                         "atmosphere_file ../data/atmmod/{atmosphere}.dat  # Location of atmospheric profile file. \n" \
                         "source solar ../data/solar_flux/kurudz_0.1nm.dat  # Location of the extraterrestrial " \
                         "spectrum\n" \
                         "mol_modify O3 {ozone} DU    # Set ozone column\n" \
                         "albedo 0.2               # Surface albedo\n" \
                         "sza {sza}               # Solar zenith angle\n" \
                         "rte_solver {solver}        # Radiative transfer equation solver\n" \
                         "number_of_streams  6     # Number of streams\n" \
                         "wavelength 250.0 5000.0   # Wavelength range [nm]\n" \
                         "altitude {altitude}\n" \
                         "aerosol_default\n" \
                         "aerosol_season {season}\n" \
                         "aerosol_angstrom {a} {aod}\n" \
                         "mol_modify H2O {iwv} MM\n".format(atmosphere=libradtran_parameter.atmosphere,
                                                            sza=libradtran_parameter.solar_zenith_angle,
                                                            altitude=libradtran_parameter.altitude_in_km,
                                                            season=libradtran_parameter.aerosol_season,
                                                            iwv=libradtran_parameter.integrated_water_vapour_in_mm,
                                                            a=libradtran_parameter.angstrom_alpha,
                                                            aod=libradtran_parameter.aerosol_optical_depth,
                                                            ozone=libradtran_parameter.ozone_column_in_du,
                                                            solver=self.rte_solver)
        fid.write(default_string)
        fid.close()
        self.__prepare_docker_skript(input_file_name, output_file_name)
        self.__run_libradtran_in_docker()
        return output_file_path

    def __run_libradtran_in_docker(self):
        libradtran_working_dir = os.path.join(os.getcwd(), self.work_directory)
        if sys.platform == 'linux':
            docker_command = r"docker run -u ${UID}" + \
                             r" -v {working_directory}:/opt/libRadtran/examples:Z" \
                             " -w /opt/libRadtran/examples  siarhei/libradtran " \
                             " /bin/bash {skript}" \
                                 .format(working_directory=libradtran_working_dir, skript=self.docker_skript_name)
        else:
            docker_command = r"docker run -v {working_directory}:" \
                             "/opt/libRadtran/examples" \
                             " -w /opt/libRadtran/examples  siarhei/libradtran  /bin/bash {skript}" \
                .format(working_directory=libradtran_working_dir, skript=self.docker_skript_name)
        print(docker_command)
        os.system(docker_command)

    def __prepare_docker_skript(self, input_file_name, output_file_name):
        fid = open(self.docker_skript_path, 'w')
        fid.write("uvspec < {inp} > {out}".format(inp=input_file_name, out=output_file_name))
        fid.close()


class LibRadtranParameter:
    def __init__(self):
        self.solar_zenith_angle = 45
        self.integrated_water_vapour_in_mm = 10
        self.angstrom_alpha = 1.5
        self.aerosol_optical_depth = 0.0048
        self.aerosol_season = 2
        self.atmosphere = 'afglms'
        self.ozone_column_in_du = 315
        self.altitude_in_km = 1.59

    def set_parameter(self, solar_zenith_angle, integrated_vater_vapour_in_mm, angstrom_alpha, aerosol_optical_depth,
                      aerosol_season=2, atmosphere='afglms', ozone_column_in_du=315, altitude_in_km=1.59):
        self.solar_zenith_angle = solar_zenith_angle
        self.integrated_water_vapour_in_mm = integrated_vater_vapour_in_mm
        self.angstrom_alpha = angstrom_alpha
        self.aerosol_optical_depth = aerosol_optical_depth
        self.aerosol_season = aerosol_season
        self.atmosphere = atmosphere
        self.ozone_column_in_du = ozone_column_in_du
        self.altitude_in_km = altitude_in_km

    def set_solar_zenith_angle(self, sza):
        self.solar_zenith_angle = sza

    def set_integrated_water_vapour_in_mm(self, iwv):
        self.integrated_water_vapour_in_mm = iwv

    def set_aod_at1000nm(self, aod):
        self.aerosol_optical_depth = aod

    def set_angstrom_alpha(self, angstrom_alpha):
        self.angstrom_alpha = angstrom_alpha
