#include "Vector3D.h"

Vector3D::Vector3D()
{
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
}

Vector3D::Vector3D(double x, double y, double z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector3D Vector3D::operator+(const Vector3D &other) const
{
    return Vector3D(x + other.x, y + other.y, z + other.z);
}

Vector3D Vector3D::operator-(const Vector3D &other) const
{
    return Vector3D(x - other.x, y - other.y, z - other.z);
}

Vector3D Vector3D::operator*(double scalar) const
{
    return Vector3D(x * scalar, y * scalar, z * scalar);
}

double Vector3D::magnitude() const
{
    return sqrt(x * x + y * y + z * z); // 计算向量模
}

Vector3D Vector3D::normalized() const
{
    double mag = magnitude();
    return Vector3D(x / mag, y / mag, z / mag); // 归一化
}
