# TriangleMeshDistance
Header only, single file, simple and efficient C++ library to compute the signed distance function to a triangle mesh.

The distance computation to the triangle collection is accelerated with a sphere bounding volume hierarchy. The signed of the distance is resolved with the method presented in *"Generating Signed Distance Fields From Triangle Meshes"* by Bærentzen, Andreas & Aanæs, Henrik. (2002).

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
- `TriangleMeshDistance` keeps a copy of the vertex and triangle data passed as input.
- Additionally, the pseudonormals required to compute signed distances are calculated and stored at building time.
- `TriangleMeshDistance` can be declared empty and constructed multiple times with different meshes. This can potentially reuse memory allocations.