#ifndef RAY_H
#define RAY_H

#include "vector.h"

class Ray
{
public:
	Ray() {};
	Ray(const Vector& o, const Vector& dir, float t = 0.0f ) : origin(o), direction(dir), time(t) {};

	Vector origin;
	Vector direction;
	float time;
};
#endif