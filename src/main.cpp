#include <iostream>
#include "vector3d.hpp"

int main()
{
	Vector3D a{ 1.1, 3.3, 2.6 };
	Vector3D b{ 2.1, 3.9, 4.7 };

	a += b;
	
	std::cout << a.x << ' ' << a.y << ' ' << a.z << '\n';
}
