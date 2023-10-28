#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>
// 定义三维向量类
class Vector3D
{
public:
    double x, y, z;                                  // 坐标 x, y, z
    Vector3D();                                      // 默认构造函数
    Vector3D(double x, double y, double z);          // 构造函数
    Vector3D operator+(const Vector3D &other) const; // 重载运算符+
    Vector3D operator-(const Vector3D &other) const; // 重载运算符-
    Vector3D operator*(double scalar) const;         // 重载运算符*
    double magnitude() const;                        // 计算向量模
    Vector3D normalized() const;                     // 向量归一化
};

#endif