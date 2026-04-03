#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include "probable.h"

int main() {
probable::SuperlatticeParams sl_params;
sl_params.delta = 0.015 * units::eV;           // Ширина минизоны
sl_params.period = 5e-6 * units::cm;           // Период сверхрешётки
sl_params.mass = 0.29 * consts::me;            // Эффективная масса
sl_params.sound_vel = 5e5 * units::cm / units::s;  // Скорость звука
sl_params.deform_pot_ac = 10 * units::eV;      // Деформационный потенциал акустических
sl_params.deform_pot_op = 1e8 * units::eV / units::cm;  // Деформационный потенциал оптических
sl_params.optical_phonon_energy = 0.035 * units::eV;  // Энергия оптического фонона
sl_params.density = 5.0 * units::g / pow(units::cm, 3);  // Плотность

probable::Material superlattice{sl_params};

double temperature = 300 * units::K;
double time_step = 1e-14 * units::s;
double all_time = 50e-12 * units::s;
size_t ensemble_size = 2000;

std::vector<double> E_values = {100, 200, 300, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000};

std::ofstream vah_file("VAH_data.txt");
vah_file << "# E(V/m) \t v_y(m/s) \t j(A/m^2)" << std::endl;

std::cout << "Calculating V-A characteristic..." << std::endl;
std::cout << "Temperature: " << temperature / units::K << " K" << std::endl;
std::cout << "Ensemble size: " << ensemble_size << std::endl;
std::cout << "Simulation time: " << all_time / units::ps << " ps" << std::endl;
std::cout << std::endl;

for (double E : E_values) {
Vec3 electric_field{0, E * units::V / units::m, 0};
Vec3 magnetic_field{0, 0, 0};

std::vector<probable::Scattering *> scattering_mechanisms{
new probable::AcousticScattering(superlattice, temperature),
new probable::OpticalScattering(superlattice, temperature, false),  // поглощение
new probable::OpticalScattering(superlattice, temperature, true)   // испускание
};

auto results = probable::simulate(superlattice,
scattering_mechanisms,
temperature,
electric_field,
magnetic_field,
time_step,
all_time,
ensemble_size,
probable::DumpFlags(probable::DumpFlags::none));

Vec3 average_velocity{0, 0, 0};
for (std::size_t i = 0; i < results.size(); ++i) {
average_velocity += (results[i].average_velocity - average_velocity) / (i + 1);
}

// j = q * n * v (n = 1e17 см^-3)
double n = 1e17 * pow(units::cm, -3);

double j_y = consts::e * n * average_velocity.y;

vah_file << std::setprecision(6)
<< E << "\t"
<< average_velocity.y / units::m * units::s << "\t"
<< j_y << std::endl;

std::cout << "E = " << E << " V/m, v_y = " << average_velocity.y / units::m * units::s << " m/s" << std::endl;

for (auto &sc : scattering_mechanisms) delete sc;
}

vah_file.close();
std::cout << "\nData saved to VAH_data.txt" << std::endl;

return 0;
}