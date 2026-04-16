#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "probable.h"

using namespace probable;

int main() {
    // распараллеливание
    omp_set_num_threads(omp_get_max_threads());
    std::cout << "Using " << omp_get_max_threads() << " OpenMP threads.\n";

    SuperlatticeParams sl_params;
    sl_params.delta = 0.015 * units::eV;
    sl_params.period = 5e-6 * units::cm;
    sl_params.mass = 0.29 * consts::me;
    sl_params.sound_vel = 5e5 * units::cm / units::s;
    sl_params.deform_pot_ac = 10 * units::eV;
    sl_params.deform_pot_op = 1e8 * units::eV / units::cm;
    sl_params.optical_phonon_energy = 0.035 * units::eV;
    sl_params.density = 5.0 * units::g / pow(units::cm, 3);

    Material superlattice{sl_params};
    double temperature = 300 * units::K;
    double time_step = 1e-14 * units::s;
    double all_time = 200e-12 * units::s;
    size_t ensemble_size = 5000;
    double frequency = 1e12;
    double omega = 2 * math::pi * frequency;

    std::cout << "=== Circular Hall Effect in Superlattices ===\n";
    std::cout << "Temperature: " << temperature / units::K << " K\n";
    std::cout << "Frequency: " << frequency / 1e12 << " THz\n";
    std::cout << "Ensemble size: " << ensemble_size << "\n\n";

    std::vector<Scattering *> mechanisms{
        new AcousticScattering(superlattice, temperature),
        new OpticalScattering(superlattice, temperature, false),
        new OpticalScattering(superlattice, temperature, true)
    };

    // === ГРАФИК 1: DC поле ===
    {
        std::cout << "[1/3] Calculating DC field dependence...\n";
        std::ofstream file("dc_field_dependence.txt");
        file << "# E_dc_x(V/m) \t E_dc_y(V/m) \t v_x(m/s) \t v_y(m/s)\n";
        std::vector<double> E_dc_values = {100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000};
        ACField ac_field(500 * units::V / units::m, 500 * units::V / units::m, omega, 0);
        Vec3 B{0, 0, 0};

        for (double E_dc : E_dc_values) {
            DCField dc_field(E_dc * units::V / units::m, E_dc * units::V / units::m);
            auto results = simulate(superlattice, mechanisms, temperature, dc_field, ac_field, B, time_step, all_time, ensemble_size, DumpFlags::none);

            double sum_x = 0, sum_y = 0, sum_z = 0;
            #pragma omp parallel for reduction(+:sum_x, sum_y, sum_z)
            for (size_t i = 0; i < results.size(); ++i) {
                sum_x += results[i].average_velocity.x;
                sum_y += results[i].average_velocity.y;
                sum_z += results[i].average_velocity.z;
            }
            Vec3 avg_v{sum_x / results.size(), sum_y / results.size(), sum_z / results.size()};

            file << E_dc << "\t" << E_dc << "\t" << avg_v.x / units::m * units::s << "\t" << avg_v.y / units::m * units::s << "\n";
            std::cout << "  E_dc=" << E_dc << " V/m | v_x=" << avg_v.x / units::m * units::s << " | v_y=" << avg_v.y / units::m * units::s << "\n";
        }
        file.close();
    }

    // === ГРАФИК 2: AC амплитуда ===
    {
        std::cout << "\n[2/3] Calculating AC amplitude dependence...\n";
        std::ofstream file("ac_field_dependence.txt");
        file << "# E0_x(V/m) \t E0_y(V/m) \t v_x(m/s) \t v_y(m/s)\n";
        std::vector<double> E0_values = {100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000};
        DCField dc_field(500 * units::V / units::m, 500 * units::V / units::m);
        Vec3 B{0, 0, 0};

        for (double E0 : E0_values) {
            ACField ac_field(E0 * units::V / units::m, E0 * units::V / units::m, omega, math::pi / 4);
            auto results = simulate(superlattice, mechanisms, temperature, dc_field, ac_field, B, time_step, all_time, ensemble_size, DumpFlags::none);

            double sum_x = 0, sum_y = 0, sum_z = 0;
            #pragma omp parallel for reduction(+:sum_x, sum_y, sum_z)
            for (size_t i = 0; i < results.size(); ++i) {
                sum_x += results[i].average_velocity.x;
                sum_y += results[i].average_velocity.y;
                sum_z += results[i].average_velocity.z;
            }
            Vec3 avg_v{sum_x / results.size(), sum_y / results.size(), sum_z / results.size()};

            file << E0 << "\t" << E0 << "\t" << avg_v.x / units::m * units::s << "\t" << avg_v.y / units::m * units::s << "\n";
            std::cout << "  E0=" << E0 << " V/m | v_x=" << avg_v.x / units::m * units::s << " | v_y=" << avg_v.y / units::m * units::s << "\n";
        }
        file.close();
    }

    // === ГРАФИК 3: Фаза ===
    {
        std::cout << "\n[3/3] Calculating phase dependence...\n";
        std::ofstream file("phase_dependence.txt");
        file << "# phi(rad) \t phi(deg) \t v_x(m/s) \t v_y(m/s)\n";
        size_t num_phases = 25;
        DCField dc_field(500 * units::V / units::m, 500 * units::V / units::m);
        double E0 = 1000;
        Vec3 B{0, 0, 0};

        for (size_t i = 0; i <= num_phases; ++i) {
            double phi = 2 * math::pi * i / num_phases;
            ACField ac_field(E0 * units::V / units::m, E0 * units::V / units::m, omega, phi);
            auto results = simulate(superlattice, mechanisms, temperature, dc_field, ac_field, B, time_step, all_time, ensemble_size, DumpFlags::none);

            double sum_x = 0, sum_y = 0, sum_z = 0;
            #pragma omp parallel for reduction(+:sum_x, sum_y, sum_z)
            for (size_t j = 0; j < results.size(); ++j) {
                sum_x += results[j].average_velocity.x;
                sum_y += results[j].average_velocity.y;
                sum_z += results[j].average_velocity.z;
            }
            Vec3 avg_v{sum_x / results.size(), sum_y / results.size(), sum_z / results.size()};

            file << phi << "\t" << phi * 180 / math::pi << "\t" << avg_v.x / units::m * units::s << "\t" << avg_v.y / units::m * units::s << "\n";
            std::cout << "  phi=" << phi*180/math::pi << "° | v_x=" << avg_v.x / units::m * units::s << " | v_y=" << avg_v.y / units::m * units::s << "\n";
        }
        file.close();
    }

    for (auto &sc : mechanisms) delete sc;
    std::cout << "\n=== All calculations completed ===\n";
    return 0;
}