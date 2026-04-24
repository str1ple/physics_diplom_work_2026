#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "probable.h"

using namespace probable;

double transverse_velocity(const Vec3& v, const DCField& dc) {
    double E_mag = sqrt(dc.E_x*dc.E_x + dc.E_y*dc.E_y);
    if (E_mag < 1e-12) return 0.0;
    // ед векторь вдоль E_dc
    Vec3 e_par = {dc.E_x / E_mag, dc.E_y / E_mag, 0.0};
    // перп ед вект (+90 град)
    Vec3 e_perp = {-e_par.y, e_par.x, 0.0};
    return v.x * e_perp.x + v.y * e_perp.y;
}

int main() {
    // распараллеливание
    omp_set_num_threads(omp_get_max_threads());
    std::cout << "Using " << omp_get_max_threads() << " OpenMP threads.\n";

    // парам сверхрешетки
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
    double frequency = 1e12; // 1 ТГц
    double omega = 2 * math::pi * frequency;

    std::cout << "=== Circular Hall Effect in Superlattices ===\n";
    std::cout << "Temperature: " << temperature / units::K << " K\n";
    std::cout << "Frequency: " << frequency / 1e12 << " THz\n";
    std::cout << "Ensemble size: " << ensemble_size << "\n\n";

    // мех. расс-я
    std::vector<Scattering *> mechanisms{
        new AcousticScattering(superlattice, temperature),
        new OpticalScattering(superlattice, temperature, false),
        new OpticalScattering(superlattice, temperature, true)
    };

    // только пост поле
    {
        std::cout << "[1/3] DC field only, transverse velocity...\n";
        std::ofstream file("dc_only_transverse.txt");
        file << "# E_dc (V/m) \t v_perp (m/s)\n";
        std::vector<double> E_vals = {100, 200, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000};
        ACField ac_field(0.0, 0.0, 0.0, 0.0);  // переменное поле отключено
        Vec3 B{0, 0, 0};

        for (double E : E_vals) {
            DCField dc_field(E * units::V / units::m, E * units::V / units::m);
            auto results = simulate(superlattice, mechanisms, temperature,
                                    dc_field, ac_field, B,
                                    time_step, all_time, ensemble_size,
                                    DumpFlags::none);

            Vec3 avg_v{0, 0, 0};
            for (const auto& res : results)
                avg_v += res.average_velocity;
            avg_v = avg_v / static_cast<double>(results.size());

            double v_perp = transverse_velocity(avg_v, dc_field);
            file << E << "\t" << v_perp / (units::m / units::s) << "\n";
            std::cout << "  E=" << E << " V/m, v_perp="
                      << v_perp / (units::m / units::s) << " m/s\n";
        }
        file.close();
    }

    // + перем поле
    {
        std::cout << "\n[2/3] AC amplitude effect on transverse velocity...\n";
        std::ofstream file("ac_amplitude_transverse.txt");
        file << "# E0 (V/m) \t v_perp (m/s)\n";
        std::vector<double> E0_vals = {100, 200, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000};
        DCField dc_field(500 * units::V / units::m, 500 * units::V / units::m);
        double phase = math::pi / 4.0;   // фиксированная фаза 45°
        Vec3 B{0, 0, 0};

        for (double E0 : E0_vals) {
            ACField ac_field(E0 * units::V / units::m, E0 * units::V / units::m,
                             omega, phase);
            auto results = simulate(superlattice, mechanisms, temperature,
                                    dc_field, ac_field, B,
                                    time_step, all_time, ensemble_size,
                                    DumpFlags::none);

            Vec3 avg_v{0, 0, 0};
            for (const auto& res : results)
                avg_v += res.average_velocity;
            avg_v = avg_v / static_cast<double>(results.size());

            double v_perp = transverse_velocity(avg_v, dc_field);
            file << E0 << "\t" << v_perp / (units::m / units::s) << "\n";
            std::cout << "  E0=" << E0 << " V/m, v_perp="
                      << v_perp / (units::m / units::s) << " m/s\n";
        }
        file.close();
    }

    // фаза
    {
        std::cout << "\n[3/3] Phase dependence of transverse velocity...\n";
        std::ofstream file("phase_transverse.txt");
        file << "# phi (rad) \t phi (deg) \t v_perp (m/s)\n";
        DCField dc_field(500 * units::V / units::m, 500 * units::V / units::m);
        double E0 = 1000.0;   // фиксированная амплитуда переменного поля
        Vec3 B{0, 0, 0};
        size_t num_phases = 18;

        for (size_t i = 0; i <= num_phases; ++i) {
            double phi = 2 * math::pi * i / num_phases;
            ACField ac_field(E0 * units::V / units::m, E0 * units::V / units::m,
                             omega, phi);
            auto results = simulate(superlattice, mechanisms, temperature,
                                    dc_field, ac_field, B,
                                    time_step, all_time, ensemble_size,
                                    DumpFlags::none);

            Vec3 avg_v{0, 0, 0};
            for (const auto& res : results)
                avg_v += res.average_velocity;
            avg_v = avg_v / static_cast<double>(results.size());

            double v_perp = transverse_velocity(avg_v, dc_field);
            file << phi << "\t" << phi * 180.0 / math::pi << "\t"
                 << v_perp / (units::m / units::s) << "\n";
            std::cout << "  phi=" << phi * 180.0 / math::pi << "°, v_perp="
                      << v_perp / (units::m / units::s) << " m/s\n";
        }
        file.close();
    }

    for (auto sc : mechanisms)
        delete sc;

    std::cout << "\n=== All calculations completed ===\n";
    return 0;
}