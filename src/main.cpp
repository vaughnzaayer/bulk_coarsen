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
#include "geometrycentral/surface/transfer_functions.h"
#include "geometrycentral/utilities/utilities.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"

#include "happly.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include <filesystem>
#include <iostream>

#include "args/args.hxx"
#include "imgui.h"

#include <Eigen/Dense>
#include <fmt/base.h>
#include <fmt/format.h>

#include "int_coarsen.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace fs = std::filesystem;

// Mesh Parameters
const int V_COUNT_THRESHOLD = 100;  // Mesh must have at least this many vertices
const int QEM_THRESHOLD = 1000000;     // Meshes above this are reduced via QEM first
const int ICE_TARGET_V_COUNT = 100; // ICE will attempt to reduce to this many vertices


// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<IntrinsicTriangulation> intrinsic;

// Some algorithm parameters
float param1 = 1.0;
double tolerance = 1.0;
int target;

double inf = std::numeric_limits<double>::max();

// Override GC's QEM implementation to have a target vertex count
void quadricErrorSimplify(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, size_t count) {
  double tol = 5.0;
  MutationManager mm(mesh, geo);
  if (mesh.nVertices() < count) return;
  
  auto toEigen = [](Vector3 v) -> Eigen::Vector3d {
    Eigen::Vector3d ret;
    ret << v.x, v.y, v.z;
    return ret;
  };
  auto fromEigen = [](Eigen::Vector3d v) -> Vector3 { return Vector3{v(0), v(1), v(2)}; };

  VertexData<Quadric> Q(mesh, Quadric());

  geo.requireFaceNormals();

  for (Face f : mesh.faces()) {
    Eigen::Vector3d n = toEigen(geo.faceNormals[f]);
    Eigen::Matrix3d M = n * n.transpose();
    for (Vertex v : f.adjacentVertices()) {
      Eigen::Vector3d q = toEigen(geo.inputVertexPositions[v]);
      double d = -n.dot(q);

      Q[v] += Quadric(M, d * n, d * d);
    }
  }

  using PotentialEdge = std::tuple<double, Edge>;

  auto cmp = [](const PotentialEdge& a, const PotentialEdge& b) -> bool { return std::get<0>(a) > std::get<0>(b); };

  std::priority_queue<PotentialEdge, std::vector<PotentialEdge>, decltype(cmp)> edgesToCheck(cmp);

  for (Edge e : mesh.edges()) {
    Quadric Qe = Q[e.halfedge().tailVertex()] + Q[e.halfedge().tipVertex()];
    Eigen::Vector3d q = Qe.optimalPoint();
    double cost = Qe.cost(q);
    edgesToCheck.push(std::make_tuple(cost, e));
  }

  while (!edgesToCheck.empty() && mesh.nVertices()>count) {
    PotentialEdge best = edgesToCheck.top();
    edgesToCheck.pop();

    // Stop when collapse becomes too expensive
    double cost = std::get<0>(best);
    if (cost > tol) break;
    if (!std::isfinite(cost)) continue; // numerical safety

    Edge e = std::get<1>(best);
    if (e.isDead()) continue; // edge no longer exists

    Vertex v1 = e.halfedge().tailVertex();
    Vertex v2 = e.halfedge().tipVertex();

    // Get edge quadric
    Quadric Qe(Q[v1], Q[v2]);
    Eigen::Vector3d q = Qe.optimalPoint();

    // If either vertex has been collapsed since the edge was pushed
    // onto the queue, the old cost was wrong. In that case, give up
    if (abs(cost - Qe.cost(q)) > 1e-8) continue;
    if (!q.array().isFinite().all()) continue; // numerical safety

    Vertex v = mm.collapseEdge(e, fromEigen(q));
    if (v == Vertex()) continue;
    Q[v] = Qe;

    for (Edge f : v.adjacentEdges()) {
      Quadric Qf(Q[f.halfedge().tailVertex()], Q[f.halfedge().tipVertex()]);
      Eigen::Vector3d q = Qf.optimalPoint();
      double cost = Qf.cost(q);
      edgesToCheck.push(std::make_tuple(cost, f));
    }
  }

  mesh.compress();
  return;
}

void loadAndProcessMesh(fs::directory_entry inputFile) {
  // Read the mesh and load it into GC
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geo;
  std::unique_ptr<IntrinsicTriangulation> itri;

  try {
    std::tie(mesh, geo) = readManifoldSurfaceMesh(inputFile.path().string());
  } catch(const std::runtime_error& e) {
    std::cerr << "Runtime exception: " << e.what() << std::endl;
    return;
  }
  

  // Check if the threshold is met 
  if (mesh->nVertices() < V_COUNT_THRESHOLD) {
    return;
  }

  // If above the QEM threshold, use QEM 
  if (mesh->nVertices() > QEM_THRESHOLD) {
    fmt::println("Simplifying via QEM");
    try {
      quadricErrorSimplify(*mesh, *geo, size_t(QEM_THRESHOLD));
    } catch(const std::runtime_error& e) {
      std::cerr << "Runtime exception: " << e.what() << std::endl;
      return;
    }
  }

  mesh->compress();
  geo->refreshQuantities();

  itri.reset(new IntegerCoordinatesIntrinsicTriangulation(*mesh, *geo));
  fmt::println("Mesh Vertex Count: {}", int(mesh->nVertices()));
  fmt::println("InTri Vertex Count: {}", int(itri->mesh.nVertices()));
  itri->flipToDelaunay();
  // Intrinsically coarsen the mesh to the target vertex count
  try {
    intrinsicallyCoarsen(*itri, size_t(ICE_TARGET_V_COUNT));
    itri->refreshQuantities();
  } catch(const std::runtime_error& e) {
    std::cerr << "Runtime exception: " << e.what() << std::endl;
    return;
  }
  

  // Clean output of dead vertices, then record curvature data to .ply file
  itri->mesh.compress();
  std::vector<Vector3> cleanPositions;
  std::vector<double> cleanCurvatures;
  std::unordered_map<Vertex, size_t> vMap;
  
  size_t idx = 0;
  for (Vertex v : itri->mesh.vertices()) {
      if (!v.isDead()) {
          cleanPositions.push_back(itri->vertexLocations[v].interpolate(geo->inputVertexPositions));
          cleanCurvatures.push_back(itri->vertexGaussianCurvatures[v]);
          vMap[v] = idx++;
      }
  }
  std::vector<std::vector<size_t>> cleanFaces;
  for (Face f : itri->mesh.faces()) {
      if (!f.isDead()) {
          std::vector<size_t> fVerts;
          for (Vertex v : f.adjacentVertices()) {
              fVerts.push_back(vMap[v]);
          }
          cleanFaces.push_back(fVerts);
      }
  }
  std::unique_ptr<ManifoldSurfaceMesh> exportMesh;
  std::unique_ptr<VertexPositionGeometry> exportGeo;

  // try {
  //   std::tie(exportMesh, exportGeo) = makeManifoldSurfaceMeshAndGeometry(cleanFaces, cleanPositions);
  //   itri.reset(new IntegerCoordinatesIntrinsicTriangulation(*exportMesh, *exportGeo));
  // } catch(const std::runtime_error& e) {
  //   std::cerr << "Runtime exception: " << e.what() << std::endl;
  //   return;
  // }
  
  mesh->compress();
  geo->refreshQuantities();
  itri->refreshQuantities();

  RichSurfaceMeshData richData(itri->mesh);
  richData.outputFormat = happly::DataFormat::ASCII; 

  itri->requireVertexGaussianCurvatures();
  richData.addVertexProperty("vertex_gaussian_curvature", itri->vertexGaussianCurvatures);
  richData.addIntrinsicGeometry(*itri);
  richData.addMeshConnectivity();
  richData.addGeometry(*geo);

  auto outFile = inputFile.path();
  outFile.replace_extension(".ply");
  richData.write(outFile.string());

  return;

}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("Takes a directory and recursively computes QEM and ICE coarsened data");
  args::Positional<std::string> inputDirectory(parser, "input", "Directory to read .obj files from.");
  args::Positional<std::string> outputDirectory(parser, "output", "Directory to write .ply files to.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure input and output directories are given
  if (!inputDirectory) {
    std::cerr << "Please specify an input directory." << std::endl;
    return EXIT_FAILURE;
  }
  if (!outputDirectory) {
    std::cerr << "Please specify an output directory." << std::endl;
    return EXIT_FAILURE;
  }


  // Read through input all inputs
  std::string in_dir = args::get(inputDirectory);
  std::string out_dir = args::get(outputDirectory);

  // Check if the output directory exists -- if not, create it
  if (!fs::exists(fs::path(out_dir))) {
    fs::create_directory(out_dir);
  }

  for (const auto& entry : fs::recursive_directory_iterator(in_dir)) {
    fmt::println("Traversing: {}", entry.path().string());

    if (fs::is_regular_file(entry) && (entry.path().extension().string() == ".obj")) {
      // If we find an obj file, check if the directory already exists
      // auto dir = fs::path()
      loadAndProcessMesh(entry);
    }
  }
  
  
  // std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(some_input_mesh));
  // intrinsic.reset(new SignpostIntrinsicTriangulation(*mesh, *geometry));
  // intrinsic->flipToDelaunay();

  return EXIT_SUCCESS;
}
