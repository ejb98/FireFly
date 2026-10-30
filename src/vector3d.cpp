#include "vector3d.hpp"

Vector3D::Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

void Vector3D::operator+=(const Vector3D& other)
{
	x += other.x;
	y += other.y;
	z += other.z;
}
