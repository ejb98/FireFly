#pragma once

struct Vector3D
{
	double x{};
	double y{};
	double z{};

	Vector3D() = default;
	Vector3D(double x = 0.0, double y = 0.0, double z = 0.0);
	void operator+=(const Vector3D& other);
};
