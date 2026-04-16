#include <iostream>
#include <random>
#include <cmath>
#include <omp.h>
#include "probable.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace probable {

double ellint_1(double k) {
    if (k <= 0) return M_PI / 2.0;
    if (k >= 1) return 10.0;
    double result = 0;
    int n = 100;
    for (int i = 0; i < n; ++i) {
        double theta = (i + 0.5) * M_PI / n;
        double denom = sqrt(1 - k * sin(theta) * sin(theta));
        if (denom > 1e-10) result += 1.0 / denom;
    }
    return result * M_PI / n;
}

double uniform() {
    // Каждый поток получает свой независимый генератор
    thread_local std::mt19937 generator(std::random_device{}());
    thread_local std::uniform_real_distribution<double> distribution(0, 1);
    return distribution(generator);
}

Vec3 Material::create_particle(double temperature) const {
    double p_max = 5 * sqrt(2 * params.mass * consts::kB * temperature);
    while (true) {
        double prob = uniform();
        double p1 = p_max * cbrt(uniform());
        if (prob < exp(-p1 * p1 / (2 * params.mass * consts::kB * temperature))) {
            double cos_theta = 1 - 2 * uniform();
            double sin_theta = sqrt(1 - cos_theta * cos_theta);
            double phi = 2 * math::pi * uniform();
            Vec3 p = {sin_theta * cos(phi), sin_theta * sin(phi), 0};
            p *= p1;
            return p;
        }
    }
}

AcousticScattering::AcousticScattering(const Material &m, double T)
    : Scattering(m, 0), temperature(T) {
    constant = pow(m.params.deform_pot_ac, 2) * m.params.mass * consts::kB * T /
               (2 * math::pi * consts::hbar * consts::hbar * m.params.density *
                m.params.sound_vel * m.params.sound_vel * m.params.delta);
}

double AcousticScattering::rate(const Vec3 &p) const {
    double D_ac = m.params.deform_pot_ac;
    double rho = m.params.density;
    double s = m.params.sound_vel;
    double prefactor = pow(D_ac, 2) * m.params.mass * consts::kB * temperature /
                       (2 * math::pi * consts::hbar * consts::hbar * rho * s * s * m.params.delta);
    double gamma = 2 * m.params.mass * m.params.delta * m.params.period * m.params.period / (consts::hbar * consts::hbar);
    double a = gamma / (m.params.delta);
    double k2 = 2.0 / (1.0 - a);
    if (k2 < 0 || k2 > 1) k2 = 0.5;
    return prefactor * ellint_1(k2);
}

double AcousticScattering::total_rate(double temperature) const { return constant; }

Vec3 AcousticScattering::scatter(const Vec3 &p) const {
    double e = m.energy(p);
    double p_new = sqrt(2 * m.params.mass * e);
    double cos_theta = 1 - 2 * uniform();
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    double phi = 2 * math::pi * uniform();
    return {p_new * sin_theta * cos(phi), p_new * sin_theta * sin(phi), 0};
}

OpticalScattering::OpticalScattering(const Material &m, double T, bool em)
    : Scattering(m, em ? m.params.optical_phonon_energy : -m.params.optical_phonon_energy),
      temperature(T), emission(em) {
    double n_phonon = 1.0 / (exp(m.params.optical_phonon_energy / (consts::kB * T)) - 1.0);
    double factor = em ? (n_phonon + 1.0) : n_phonon;
    constant = pow(m.params.deform_pot_op, 2) * m.params.mass * factor /
               (2 * math::pi * consts::hbar * consts::hbar * m.params.density * m.params.optical_phonon_energy);
}

double OpticalScattering::rate(const Vec3 &p) const {
    double D_op = m.params.deform_pot_op;
    double rho = m.params.density;
    double omega = m.params.optical_phonon_energy;
    double n_phonon = 1.0 / (exp(omega / (consts::kB * temperature)) - 1.0);
    double factor = emission ? (n_phonon + 1.0) : n_phonon;
    double prefactor = pow(D_op, 2) * m.params.mass / (2 * math::pi * consts::hbar * consts::hbar * rho * omega);
    if (emission && m.energy(p) < omega) return 0;
    return prefactor * factor;
}

double OpticalScattering::total_rate(double temperature) const { return constant; }

Vec3 OpticalScattering::scatter(const Vec3 &p) const {
    double e = m.energy(p) - energy_change;
    if (e < 0) return p;
    double p_new = sqrt(2 * m.params.mass * e);
    double cos_theta = 1 - 2 * uniform();
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    double phi = 2 * math::pi * uniform();
    return {p_new * sin_theta * cos(phi), p_new * sin_theta * sin(phi), 0};
}

void Results::append(uint32_t n, double t, const Vec3 &p, const Vec3 &v, double e, size_t s) {
    average_velocity += (v - average_velocity) / (n + 1);
    if (s) scattering_count[s - 1] += 1;
    if (!(flags & DumpFlags::on_scatterings) || s) {
        if (flags & DumpFlags::scattering) scatterings.push_back(s);
        if (flags & DumpFlags::number) ns.push_back(n);
        if (flags & DumpFlags::time) ts.push_back(t);
        if (flags & DumpFlags::momentum) momentums.push_back(p);
        if (flags & DumpFlags::velocity) velocities.push_back(v);
        if (flags & DumpFlags::energy) energies.push_back(e);
        size += 1;
    }
}

std::ostream &operator<<(std::ostream &s, const Results &r) {
    if (r.flags == DumpFlags::none) return s;
    for (size_t i = 0; i < r.size; ++i) {
        if (r.flags & DumpFlags::number) s << r.ns[i] << " ";
        if (r.flags & DumpFlags::time) s << r.ts[i] / units::s << " ";
        if (r.flags & DumpFlags::momentum) s << r.momentums[i].x << " " << r.momentums[i].y << " " << r.momentums[i].z << " ";
        if (r.flags & DumpFlags::velocity) s << r.velocities[i].x / units::m * units::s << " " << r.velocities[i].y / units::m * units::s << " " << r.velocities[i].z / units::m * units::s << " ";
        if (r.flags & DumpFlags::energy) s << r.energies[i] / units::eV << " ";
        if (r.flags & DumpFlags::scattering) s << r.scatterings[i];
        s << "\n";
    }
    return s;
}

// === РАСПАРАЛЛЕЛЕННАЯ ВЕРСИЯ ===
std::vector<Results> simulate(const Material &material,
                              const std::vector<Scattering *> mechanisms,
                              double temperature,
                              const DCField &dc_field,
                              const ACField &ac_field,
                              const Vec3 &magnetic_field,
                              double time_step,
                              double all_time,
                              size_t ensemble_size,
                              DumpFlags flags) {
    std::vector<Results> results(ensemble_size);
    size_t steps = all_time / time_step + 1;
    size_t alloc = (flags & DumpFlags::on_scatterings) ? steps / 10 : steps;

    // Распараллеливание по ансамблю частиц
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < ensemble_size; ++i) {
        Results result;
        if (flags != DumpFlags::none) result = Results(alloc, flags);
        result.average_velocity = {0, 0, 0};
        result.scattering_count.assign(mechanisms.size(), 0);

        Vec3 p = material.create_particle(temperature);
        std::vector<double> free_flight(mechanisms.size(), 0);
        for (double &l : free_flight) l = -log(uniform());

        size_t steady_steps = (10e-12 * units::s) / time_step;
        Vec3 v_sum{0, 0, 0};
        size_t v_count = 0;

        for (size_t j = 0; j < steps; ++j) {
            double t = j * time_step;
            Vec3 total_field = dc_field.field() + ac_field.field_at_time(t);
            Vec3 p_ = p;
            Vec3 v = material.velocity(p_);
            double e = material.energy(p_);

            if (j > steady_steps) {
                v_sum += v;
                v_count++;
            }

            size_t scattering_mechanism = 0;
            for (size_t k = 0; k < mechanisms.size(); ++k) {
                free_flight[k] -= mechanisms[k]->rate(p_) * time_step;
                if (free_flight[k] < 0) {
                    p = mechanisms[k]->scatter(p_);
                    free_flight[k] = -log(uniform());
                    scattering_mechanism = k + 1;
                    break;
                }
            }

            result.append(j, t, p_, v, e, scattering_mechanism);
            Vec3 lorentz_force = magnetic_field.cross(v);
            p += -consts::e * (total_field + lorentz_force) * time_step;
        }

        if (v_count > 0) result.average_velocity = v_sum / v_count;
        results[i] = std::move(result); // Эффективное перемещение без копирования
    }

    return results;
}

} // namespace probable