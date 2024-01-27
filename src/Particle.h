#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vector3D.h"

class Particle
{
public:
    Vector3D position;
    Particle();                         // 默认构造函数的声明
    Particle(const Vector3D &position); // 带参数的构造函数的声明
};

#endif