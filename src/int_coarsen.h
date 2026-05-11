#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/quadric_error_simplification.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/boundary_first_flattening.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/utilities/utilities.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include "args/args.hxx"
#include "imgui.h"

#include <Eigen/Dense>
#include <fmt/base.h>
#include <fmt/format.h>

namespace gc  = geometrycentral;
namespace gcs = geometrycentral::surface;

struct VertexCost;

bool flattenVertexCETM(gcs::IntrinsicTriangulation& tri, gcs::Vertex i);

bool flatVertexRemoval(gcs::IntrinsicTriangulation& tri, gcs::Vertex i); 

void updateCurvatureChannels(gcs::IntrinsicTriangulation& tri, gcs::VertexData<double>& posKVal, gcs::VertexData<double>& negKVal);

// void singleChannelUpdate(IntrinsicTriangulation& tri, Vertex i, VertexData<double>& k, VertexData<Vector2> tanVect, VertexData<double>& removalCost);

// void singleChannelCost(IntrinsicTriangulation& tri, Vertex i, VertexData<double>& k, VertexData<Vector2>& tanVect, VertexData<double>& removalCost);

void computeRemovalCost(gcs::IntrinsicTriangulation& tri, gcs::Vertex i, gcs::VertexData<double>& posKVal, gcs::VertexData<double>& negKVal, gcs::VertexData<gc::Vector2>& posKTan, gcs::VertexData<gc::Vector2>& negKTan, gcs::VertexData<double>& removalCost, bool updateTans, bool revertFlattening);

bool removeVertAndUpdate(gcs::IntrinsicTriangulation& tri, gcs::Vertex i, gcs::VertexData<double>& posKVal, gcs::VertexData<double>& negKVal, gcs::VertexData<gc::Vector2>& posKTan, gcs::VertexData<gc::Vector2>& negKTan, gcs::VertexData<double>& removalCost);

void intrinsicallyCoarsen(gcs::IntrinsicTriangulation& tri, size_t count);