#pragma once
#include <cmath>
#include <iostream>

// вект
struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    double dot(const Vec3 &other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    double length() const {
        return std::sqrt(this->dot(*this));
    }

    Vec3 cross(const Vec3 &other) const {
        return {
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.z
        };
    }

    Vec3 sqrt() const {
        return {std::sqrt(x), std::sqrt(y), std::sqrt(z)};
    }

    Vec3 &operator+=(const Vec3 &other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    friend Vec3 operator+(Vec3 a, const Vec3 &b) { return a += b; }

    Vec3 &operator-=(const Vec3 &other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    friend Vec3 operator-(Vec3 a, const Vec3 &b) { return a -= b; }

    Vec3 &operator*=(const Vec3 &other) {
        x *= other.x; y *= other.y; z *= other.z;
        return *this;
    }

    friend Vec3 operator*(Vec3 a, const Vec3 &b) { return a *= b; }

    Vec3 &operator*=(double other) {
        x *= other; y *= other; z *= other;
        return *this;
    }

    friend Vec3 operator*(Vec3 a, double b) { return a *= b; }
    friend Vec3 operator*(double b, Vec3 a) { return a *= b; }

    Vec3 &operator/=(double other) {
        x /= other; y /= other; z /= other;
        return *this;
    }

    friend Vec3 operator/(Vec3 a, double b) { return a /= b; }

    friend std::ostream &operator<<(std::ostream &s, const Vec3 &v) {
        s << "{" << v.x << ", " << v.y << ", " << v.z << "}";
        return s;
    }
};