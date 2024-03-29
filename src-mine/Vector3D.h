#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>

class Vector3D
{
public:
    double x, y, z;
    Vector3D()
    {
        this->x = 0.0;
        this->y = 0.0;
        this->z = 0.0;
    }

    Vector3D(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // 两个向量相加
    Vector3D operator+(const Vector3D &other) const
    {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }

    // 两个向量相减
    Vector3D operator-(const Vector3D &other) const
    {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }

    // 向量乘以标量
    Vector3D operator*(double scalar) const
    {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    // 向量的模
    double magnitude() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    // 向量点乘向量
    double dot(const Vector3D &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }
};

#endif
