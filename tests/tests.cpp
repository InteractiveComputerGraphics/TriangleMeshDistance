#include "catch.hpp"
#include <array>
#include <vector>

#include <tmd/TriangleMeshDistance.h>

double point_AABB_signed(const tmd::Vec3d& point, const tmd::Vec3d& bottom, const tmd::Vec3d& top)
{
	const tmd::Vec3d dx = { std::max(bottom[0] - point[0], point[0] - top[0]),
									std::max(bottom[1] - point[1], point[1] - top[1]),
									std::max(bottom[2] - point[2], point[2] - top[2]) };

	const double max_dx = std::max(dx[0], std::max(dx[1], dx[2]));
	if (max_dx < 0.0) { // Inside
		return max_dx;
	}
	else { // Outside
		double dist_sq = 0.0;
		for (int i = 0; i < 3; i++) {
			if (dx[i] > 0.0) {
				dist_sq += dx[i] * dx[i];
			}

		}
		return std::sqrt(dist_sq);
	}
}


TEST_CASE("TriangleMeshDistance", "")
{
	SECTION("SDF to a cube")
	{
		std::vector<tmd::Vec3d> vertices = { { 1, -1, -1 }, { 1, 0, -1 }, { 1, 1, -1 }, { 1, -1, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 1, -1, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { -1, -1, -1 }, { -1, 0, -1 }, { -1, 1, -1 }, { -1, -1, 0 }, { -1, 0, 0 }, { -1, 1, 0 }, { -1, -1, 1 }, { -1, 0, 1 }, { -1, 1, 1 }, { 0, 1, -1 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, -1, -1 }, { 0, -1, 0 }, { 0, -1, 1 }, { 0, 0, 1 }, { 0, 0, -1 } };
		std::vector<std::array<int, 3>> connectivity = { { 0, 1, 3 }, { 1, 4, 3 }, { 1, 2, 4 }, { 2, 5, 4 }, { 3, 4, 6 }, { 4, 7, 6 }, { 4, 5, 7 }, { 5, 8, 7 }, { 12, 10, 9 }, { 12, 13, 10 }, { 13, 11, 10 }, { 13, 14, 11 }, { 15, 13, 12 }, { 15, 16, 13 }, { 16, 14, 13 }, { 16, 17, 14 }, { 14, 18, 11 }, { 14, 19, 18 }, { 19, 2, 18 }, { 19, 5, 2 }, { 17, 19, 14 }, { 17, 20, 19 }, { 20, 5, 19 }, { 20, 8, 5 }, { 9, 21, 12 }, { 21, 22, 12 }, { 21, 0, 22 }, { 0, 3, 22 }, { 12, 22, 15 }, { 22, 23, 15 }, { 22, 3, 23 }, { 3, 6, 23 }, { 15, 23, 16 }, { 23, 24, 16 }, { 23, 6, 24 }, { 6, 7, 24 }, { 16, 24, 17 }, { 24, 20, 17 }, { 24, 7, 20 }, { 7, 8, 20 }, { 10, 21, 9 }, { 10, 25, 21 }, { 25, 0, 21 }, { 25, 1, 0 }, { 11, 25, 10 }, { 11, 18, 25 }, { 18, 1, 25 }, { 18, 2, 1 } };

		tmd::TriangleMeshDistance mesh_distance(vertices, connectivity);

		for (double x = -2; x < 2; x += 0.13) {
			for (double y = -2; y < 2; y += 0.13) {
				for (double z = -2; z < 2; z += 0.13) {
					const auto result = mesh_distance.signed_distance({ x, y, z });
					const double exact = point_AABB_signed({ x, y, z }, { -1, -1, -1 }, { 1, 1, 1 });
					REQUIRE(std::abs(result.distance - exact) < 1e-10);
				}
			}
		}
	};
}