# TriangleMeshDistance
Header only, single file, simple and efficient C++11 library to compute the signed distance function (SDF) to a triangle mesh.

The distance computation to the triangle collection is accelerated with a sphere bounding volume hierarchy. The sign of the distance is resolved with the method presented in *"Generating Signed Distance Fields From Triangle Meshes"* by Bærentzen, Andreas & Aanæs, Henrik. (2002).

Assuming triangle normals point outwards from the enclosed volume, the sign of the distance will be positive outside and negative inside.

## Example
```cpp
// Declare mesh vertices and triangles
std::vector<std::array<double, 3>> vertices;
std::vector<std::array<int, 3>> triangles;

// (... fill the `vertices` and `triangles` with the mesh data ...)

// Initialize TriangleMeshDistance
tmd::TriangleMeshDistance mesh_distance(vertices, triangles);

// Query TriangleMeshDistance
tmd::Result result = mesh_distance.signed_distance({ x, y, z });

// Print result
std::cout << "Signed distance: " << result.distance << std::endl;
std::cout << "Nearest point: " << result.nearest_point << std::endl;
std::cout << "Nearest entity: " << result.nearest_entity << std::endl;
std::cout << "Nearest triangle index: " << result.triangle_id << std::endl;
```

## What you need to know about TriangleMeshDistance
- The input triangle mesh must be fully connected and watertight. Triangle soups and meshes with holes will return the correct distance but the sign will be undefined.
- Triangle winding (consistent normals orientation) is not verified. The input mesh is required to have consistent normals.
- `TriangleMeshDistance` keeps a copy of the vertex and triangle data.
- The pseudonormals required to compute signed distances are calculated and stored at building time.
- `TriangleMeshDistance` can be declared empty and constructed multiple times with different meshes. If the new mesh needs less memory than the curent one, memory allocations will be avoided.

## Projects using TriangleMeshDistance

- [Discregrid](https://github.com/InteractiveComputerGraphics/Discregrid) - A static C++ library for the generation of discrete functions on a box-shaped domain. This is especially suited for the discretization of signed distance fields.
- [PBD](https://github.com/InteractiveComputerGraphics/PositionBasedDynamics) - A C++ library for physically-based simulation of rigid bodies, deformables, cloth and fluids using Position-Based Dynamics.
- [SPlisHSPlasH](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH) - A C++ library for the physically-based simulation of fluids using Smoothed Particle Hydrodynamics.
