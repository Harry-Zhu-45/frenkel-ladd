#include "Particle.h"

Particle::Particle() : position() {}                                 // 默认构造函数的定义
Particle::Particle(const Vector3D &position) : position(position) {} // 带参数的构造函数的定义
