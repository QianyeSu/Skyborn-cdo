// modified code from:
// https://schneide.wordpress.com/2016/07/15/generating-an-icosphere-in-c

#include "vector3d.h"
#include "varray.h"
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include <utility>
#include <limits>
#include <array>
#include <vector>

using Index = uint64_t;
using Vertex = Vector3d;
using Triangle = std::array<size_t, 3>;
using TriangleList = std::vector<Triangle>;
using VertexList = std::vector<Vertex>;

// flat_map needs g++ version 15 or clang++ version 20
// #define USE_FLAT_MAP 1
#ifdef USE_FLAT_MAP
#include <flat_map>
using Lookup = std::vector<std::flat_map<Index, Index>>;
#else
#include <map>
using Lookup = std::vector<std::map<Index, Index>>;
#endif

// icosahedron from ICON
namespace icosahedron
{

// Northern hemisphere are the first 6 elements of vertices[0:5]
// Southern hemisphere are the other 6 elements of vertices[6:11]
// 12 vertices
static VertexList vertices(12);

static void
init(void)
{
  // first define the vertices of the icosahedron
  // set poles first
  vertices[0] = Vertex{ 0.0, 0.0, 1.0 };
  vertices[11] = Vertex{ 0.0, 0.0, -1.0 };
  // now set the vertices on the two latitude rings
  int i_mdist[10];
  for (int j = 1; j < 11; ++j) { i_mdist[(j % 2 == 0) ? (j / 2 + 4) : ((j + 1) / 2 - 1)] = -1 + (j - 1) - 10 * ((j - 1) / 7); }

  constexpr auto pi_5 = M_PI * 0.2;
  auto z_w = 2.0 * std::acos(1.0 / (2.0 * std::sin(pi_5)));
  for (int j = 1; j < 11; ++j)
  {
    // toggle the hemisphere
    auto i_msgn = (j >= 6) ? -1.0 : 1.0;
    // compute the meridian angle for the base vertex.
    auto z_rlon = (1.0 + i_mdist[j - 1]) * pi_5;
    // now initialize the coordinates
    vertices[j] = Vertex{ std::sin(z_w) * std::cos(z_rlon), std::sin(z_w) * std::sin(z_rlon), std::cos(z_w) * i_msgn };
  }
}

// 20 triangles
static const TriangleList triangles
    = { { { 0, 1, 2 } },  { { 0, 2, 3 } },  { { 0, 3, 4 } },  { { 0, 4, 5 } },   { { 0, 5, 1 } },
        { { 6, 2, 1 } },  { { 7, 3, 2 } },  { { 8, 4, 3 } },  { { 9, 5, 4 } },   { { 10, 1, 5 } },
        { { 2, 6, 7 } },  { { 3, 7, 8 } },  { { 4, 8, 9 } },  { { 5, 9, 10 } },  { { 1, 10, 6 } },
        { { 11, 7, 6 } }, { { 11, 8, 7 } }, { { 11, 9, 8 } }, { { 11, 10, 9 } }, { { 11, 6, 10 } } };

}  // namespace icosahedron

static Index
vertex_for_edge(Lookup &lookup, VertexList &vertices, Index first, Index second)
{
  if (first > second) std::swap(first, second);

  auto [it, success] = lookup[first].insert({ second, vertices.size() });
  if (success) vertices.push_back((vertices[first] + vertices[second]).normalised());

  return it->second;
}

static TriangleList
subdivide(VertexList &vertices, TriangleList const &triangles)
{
  auto n = triangles.size();
  Lookup lookup(n / 2 + 1);
  TriangleList result(4 * n);
  Triangle mid;
  for (size_t i = 0; i < n; ++i)
  {
    auto const &each = triangles[i];
    for (int edge = 0; edge < 3; edge++) { mid[edge] = vertex_for_edge(lookup, vertices, each[edge], each[(edge + 1) % 3]); }

    result[i * 4 + 0] = { each[0], mid[0], mid[2] };
    result[i * 4 + 1] = { each[1], mid[1], mid[0] };
    result[i * 4 + 2] = { each[2], mid[2], mid[1] };
    result[i * 4 + 3] = mid;
  }

  return result;
}

size_t
gen_icosphere_coords(int subdivisions, bool withBounds, Varray<double> &xvals, Varray<double> &yvals, Varray<double> &xbounds,
                     Varray<double> &ybounds)
{
  icosahedron::init();
  auto triangles = icosahedron::triangles;
  auto vertices = icosahedron::vertices;

  size_t numTriangles = std::pow(4, subdivisions) * 20;
  if (numTriangles > (size_t) std::numeric_limits<Index>::max())
  {
    std::fprintf(stderr, "Too many grid cells:%zu (limit=%zu)!\n", numTriangles, (size_t) std::numeric_limits<Index>::max());
    std::exit(EXIT_FAILURE);
  }

  size_t numVerticies = std::pow(4, subdivisions) * 10 + 2;
  vertices.reserve(numVerticies);
  // printf("numTriangles %zu, numVerticies %zu\n", numTriangles, numVerticies);

  for (int i = 0; i < subdivisions; ++i) triangles = subdivide(vertices, triangles);

  auto numCells = triangles.size();
  xvals.resize(numCells);
  yvals.resize(numCells);
  if (withBounds)
  {
    xbounds.resize(3 * numCells);
    ybounds.resize(3 * numCells);
  }

#ifdef _OPENMP
#pragma omp parallel for if (numCells > 999999) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < numCells; ++i)
  {
    auto const &t = triangles[i];
    auto center = circum_center_mean(vertices[t[0]], vertices[t[1]], vertices[t[2]]);
    xvals[i] = center.longitude();
    yvals[i] = center.latitude();
    if (withBounds)
      for (size_t k = 0; k < 3; ++k)
      {
        xbounds[i * 3 + k] = vertices[t[k]].longitude();
        ybounds[i * 3 + k] = vertices[t[k]].latitude();
      }
  }

  return numCells;
}
