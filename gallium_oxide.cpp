#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <sstream>
#include "probable.h"

using namespace probable;

// Функция для разложения скорости на продольную и поперечную относительно DC поля
void decompose_velocity(const Vec3& v, const DCField& dc, double& v_long, double& v_trans) {
    double E_mag = sqrt(dc.E_x*dc.E_x + dc.E_y*dc.E_y);
    if (E_mag < 1e-12) {
        v_long = 0.0;
        v_trans = 0.0;
        return;
    }
    // Единичный вектор вдоль E_dc (направлен под 45° к осям, так как E_x = E_y)
    Vec3 e_par = {dc.E_x / E_mag, dc.E_y / E_mag, 0.0};
    // Единичный вектор перпендикулярно
    Vec3 e_perp = {-e_par.y, e_par.x, 0.0};

    v_long = v.x * e_par.x + v.y * e_par.y;
    v_trans = v.x * e_perp.x + v.y * e_perp.y;
}

int main() {
    omp_set_num_threads(omp_get_max_threads());
    std::cout << "Используется " << omp_get_max_threads() << " OpenMP ядер.\n";

    // Параметры сверхрешетки
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
    Vec3 B{0, 0, 0};

    // Механизмы рассеяния
    std::vector<Scattering *> mechanisms{
        new AcousticScattering(superlattice, temperature),
        new OpticalScattering(superlattice, temperature, false),
        new OpticalScattering(superlattice, temperature, true)
    };

    // НАЧАЛО ЗАПОЛНЕНИЯ E_DC_VALS

    std::vector<double> E_dc_vals = {0, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000};
    // int n_points = 100;
    // double start = 0.0;
    // double end = 10000.0;
    // std::vector<double> E_dc_vals(n_points);
    // double step = (end - start) / (n_points - 1);
    // for (int i = 0; i < n_points; ++i) {
    //     E_dc_vals[i] = start + i * step;
    // }

    // КОНЕЦ ЗАПОЛНЕНИЯ E_DC_VALS

    // Вспомогательная функция для запуска симуляции
    auto run_sim = [&](const DCField& dc, const ACField& ac, Vec3& out_v) {
        auto results = simulate(superlattice, mechanisms, temperature,
                                dc, ac, B,
                                time_step, all_time, ensemble_size,
                                DumpFlags::none);
        Vec3 avg_v{0, 0, 0};
        for (const auto& res : results) avg_v += res.average_velocity;
        avg_v = avg_v / static_cast<double>(results.size());
        out_v = avg_v;
    };

    // 1. Сравнение при изменении амплитуды E0
    std::cout << "\n[1/3] СЕРИЯ 1: Изменение амплитуды переменного поля...\n";
    std::vector<double> E0_vals = {500, 1500, 3000, 6000};
    double freq1 = 1.0e12;
    double phase1 = math::pi / 4.0; // 45 градусов

    for (double E0 : E0_vals) {
        std::string fname = "comp_s1_E0_" + std::to_string(static_cast<int>(E0)) + ".txt";
        std::ofstream file(fname);
        file << "# E_dc (V/m) \t v_trans_DC_only (m/s) \t v_trans_DC_plus_AC (m/s)\n";
        std::cout << "  Вычисление для E0 = " << E0 << " В/м...\n";

        for (double E_dc_val : E_dc_vals) {
            // Поле под 45 градусов
            DCField dc_field(E_dc_val * units::V / units::m, E_dc_val * units::V / units::m);

            // 1. Только постоянное поле
            ACField ac_zero(0, 0, 0, 0);
            Vec3 v1;
            run_sim(dc_field, ac_zero, v1);
            double v1_long, v1_trans;
            decompose_velocity(v1, dc_field, v1_long, v1_trans);

            // 2. Постоянное + Переменное
            ACField ac_field(E0 * units::V / units::m, E0 * units::V / units::m,
                           2 * math::pi * freq1, phase1);
            Vec3 v2;
            run_sim(dc_field, ac_field, v2);
            double v2_long, v2_trans;
            decompose_velocity(v2, dc_field, v2_long, v2_trans);

            file << E_dc_val << "\t"
                 << v1_trans / (units::m / units::s) << "\t"
                 << v2_trans / (units::m / units::s) << "\n";
        }
        file.close();
    }

    // 2. Сравнение при изменении частоты f
    std::cout << "\n[2/3] СЕРИЯ 2: Изменение частоты переменного поля...\n";
    std::vector<double> freq_vals_THz = {0.1, 0.5, 1.0, 2.0};
    double E0_fixed = 3000.0;
    double phase2 = math::pi / 4.0;

    for (double f_THz : freq_vals_THz) {
        std::string fname = "comp_s2_f_" + std::to_string(static_cast<int>(f_THz*10)) + ".txt";
        std::ofstream file(fname);
        file << "# E_dc (V/m) \t v_trans_DC_only (m/s) \t v_trans_DC_plus_AC (m/s)\n";
        std::cout << "  Вычисление для f = " << f_THz << " ТГц...\n";

        for (double E_dc_val : E_dc_vals) {
            DCField dc_field(E_dc_val * units::V / units::m, E_dc_val * units::V / units::m);

            // 1. Только постоянное поле
            ACField ac_zero(0, 0, 0, 0);
            Vec3 v1;
            run_sim(dc_field, ac_zero, v1);
            double v1_trans;
            { double l; decompose_velocity(v1, dc_field, l, v1_trans); }

            // 2. Постоянное + Переменное
            ACField ac_field(E0_fixed * units::V / units::m, E0_fixed * units::V / units::m,
                           2 * math::pi * f_THz * 1e12, phase2);
            Vec3 v2;
            run_sim(dc_field, ac_field, v2);
            double v2_trans;
            { double l; decompose_velocity(v2, dc_field, l, v2_trans); }

            file << E_dc_val << "\t"
                 << v1_trans / (units::m / units::s) << "\t"
                 << v2_trans / (units::m / units::s) << "\n";
        }
        file.close();
    }

    // 3. Сравнение при изменении фазы phi
    std::cout << "\n[3/3] СЕРИЯ 3: Изменение фазы переменного поля...\n";
    std::vector<double> phase_vals_deg = {0, 45, 90, 150, 240, 360};
    double E0_fixed3 = 3000.0;
    double freq3 = 1.0e12;

    for (double phi_deg : phase_vals_deg) {
        double phi = phi_deg * math::pi / 180.0;
        std::string fname = "comp_s3_phi_" + std::to_string(static_cast<int>(phi_deg)) + ".txt";
        std::ofstream file(fname);
        file << "# E_dc (V/m) \t v_trans_DC_only (m/s) \t v_trans_DC_plus_AC (m/s)\n";
        std::cout << "  Вычисление для phi = " << phi_deg << " градусов...\n";

        for (double E_dc_val : E_dc_vals) {
            DCField dc_field(E_dc_val * units::V / units::m, E_dc_val * units::V / units::m);

            // 1. Только постоянное поле
            ACField ac_zero(0, 0, 0, 0);
            Vec3 v1;
            run_sim(dc_field, ac_zero, v1);
            double v1_trans;
            { double l; decompose_velocity(v1, dc_field, l, v1_trans); }

            // 2. Постоянное + Переменное
            ACField ac_field(E0_fixed3 * units::V / units::m, E0_fixed3 * units::V / units::m,
                           2 * math::pi * freq3, phi);
            Vec3 v2;
            run_sim(dc_field, ac_field, v2);
            double v2_trans;
            { double l; decompose_velocity(v2, dc_field, l, v2_trans); }

            file << E_dc_val << "\t"
                 << v1_trans / (units::m / units::s) << "\t"
                 << v2_trans / (units::m / units::s) << "\n";
        }
        file.close();
    }

    for (auto sc : mechanisms) delete sc;
    std::cout << "\n=== Все вычисления успешно завершены. Можно строить графики ===\n";
    return 0;
}