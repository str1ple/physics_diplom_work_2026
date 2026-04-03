#pragma once
#include <iostream>
#include <cxxabi.h>
#include <vector>
#include <cmath>
#include "vec3.h"

namespace units {
const double eV = 1;
const double J = 1 / 1.6e-19;
const double um = 1;
const double cm = um * 1e4;
const double m = um * 1e6;
const double nm = um * 1e-3;
const double ps = 1;
const double s = ps * 1e12;
const double kg = J * s * s / m / m;
const double g = 1e-3 * kg;
const double V = 1;
const double T = 3e8 * (V / m) / (m / s);
const double C = J / V;
const double K = 1.38e-23 / 1.6e-19;
}

namespace consts {
const double hbar = 1.0545718e-34 * units::J * units::s;
const double e = 1.6021766e-19 * units::C;
const double c = 3e8 * units::m / units::s;
const double me = 9.1093835e-31 * units::kg;
const double kB = 1.3806485e-23 * units::J / units::K;
const double eps0 = 8.8541878e-12 * units::C / units::V / units::m;
}

namespace math {
const double pi = acos(-1);
}

namespace probable {
struct SuperlatticeParams {
double delta;              // Ширина минизоны (эВ)
double period;             // Период сверхрешётки d (см)
double mass;               // Эффективная масса
double sound_vel;          // Скорость звука (см/с)
double deform_pot_ac;      // Деформационный потенциал акустических (эВ)
double deform_pot_op;      // Деформационный потенциал оптических (эВ/см)
double optical_phonon_energy; // Энергия оптического фонона (эВ)
double density;            // Плотность (г/см³)
};

struct Material {
SuperlatticeParams params;

// 2.57
double energy(const Vec3 &p) const {
    double p_y_d = p.y * params.period / consts::hbar;
    double kinetic_x = p.x * p.x / (2.0 * params.mass);
    double miniband = params.delta * (1.0 - cos(p_y_d));
    return kinetic_x + miniband;
}

Vec3 velocity(const Vec3 &p) const {
    double v_x = p.x / params.mass;
    double p_y_d = p.y * params.period / consts::hbar;
    double v_y = params.delta * params.period / consts::hbar * sin(p_y_d);
    return {v_x, v_y, 0};
}

Vec3 create_particle(double temperature) const;
};

struct Scattering {
const Material &m;
const double energy_change;
Scattering(const Material &m, double e) : m(m), energy_change(e) {}
virtual double rate(const Vec3 &p) const = 0;
virtual double total_rate(double temperature) const { return 0; }
virtual Vec3 scatter(const Vec3 &p) const = 0;
virtual ~Scattering() {}
};

struct AcousticScattering : public Scattering {
double constant;
double temperature;
AcousticScattering(const Material &m, double T);
double rate(const Vec3 &p) const override;
double total_rate(double temperature) const override;
Vec3 scatter(const Vec3 &p) const override;
};

struct OpticalScattering : public Scattering {
double constant;
double temperature;
bool emission;
OpticalScattering(const Material &m, double T, bool em);
double rate(const Vec3 &p) const override;
double total_rate(double temperature) const override;
Vec3 scatter(const Vec3 &p) const override;
};

inline std::ostream &operator<<(std::ostream &s, const Scattering &sc) {
int status;
char *realname = abi::__cxa_demangle(typeid(sc).name(), 0, 0, &status);
s << realname << " dE=" << sc.energy_change;
free(realname);
return s;
}

enum DumpFlags {
none = 0,
number = 1,
time = number << 1,
momentum = time << 1,
energy = momentum << 1,
velocity = energy << 1,
scattering = velocity << 1,
all = number | time | momentum | energy | velocity | scattering,
on_scatterings = scattering << 1,
};

struct Results {
size_t size;
DumpFlags flags;
Vec3 average_velocity;
std::vector<uint32_t> scattering_count;
std::vector<uint32_t> ns;
std::vector<double> ts;
std::vector<Vec3> momentums;
std::vector<Vec3> velocities;
std::vector<double> energies;
std::vector<uint32_t> scatterings;

Results() : size(0), flags(DumpFlags::none) {}
Results(size_t cap, DumpFlags flags = DumpFlags::none)
    : size(0), flags(flags), average_velocity{0,0,0} {
    ns.reserve(cap);
    ts.reserve(cap);
    momentums.reserve(cap);
    velocities.reserve(cap);
    energies.reserve(cap);
    scatterings.reserve(cap);
}

void append(uint32_t n, double t, const Vec3 &p, const Vec3 &v, double e, size_t s);
friend std::ostream &operator<<(std::ostream &s, const Results &r);
};

std::vector<Results> simulate(const Material &material,
const std::vector<Scattering *> mechanisms,
double temperature,
const Vec3 &electric_field,
const Vec3 &magnetic_field,
double time_step,
double all_time,
size_t ensemble_size,
DumpFlags flags = DumpFlags::all);

double uniform();

double ellint_1(double k);

} // namespace probable