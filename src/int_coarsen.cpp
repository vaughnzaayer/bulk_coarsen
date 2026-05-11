#include "int_coarsen.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct VertexCost {
    public:
    Vertex v;
    double cost;
    bool operator>(const VertexCost& other) const {
        return cost > other.cost;
    }
};

bool flattenVertexCETM(IntrinsicTriangulation& tri, Vertex i) {
    const int MAX_ITERS = 10;
    const double TOLERANCE = 1e-6;

    double ui = 0.0;
    std::vector<Halfedge> halfedges;
    std::vector<double> originalLens;
    for (Halfedge he : i.outgoingHalfedges()) {
        halfedges.push_back(he);
        originalLens.push_back(tri.edgeLengths[he.edge()]);
    }

    for (int iter = 0; iter < MAX_ITERS; iter++) {
        double scale = std::exp(ui / 2.0);
        for (size_t k = 0; k < halfedges.size(); k++) {
            tri.edgeLengths[halfedges[k].edge()] = originalLens[k] * scale;
        }


        bool valid = true;
        double currentAngleSum = 0;
        double cotanSum = 0;

        for (Corner c : i.adjacentCorners()) {
            Halfedge he = c.halfedge();

            double l1 = tri.edgeLengths[he.edge()];
            double l2 = tri.edgeLengths[he.next().edge()];
            double l3 = tri.edgeLengths[he.next().next().edge()];

            // Check triangle inequality
            if (l1 + l2 <= l3 || l2 + l3 <= l1 || l3 + l1 <= l2) {
                valid = false;
                break;
            }

            double cosTheta_i = (l1*l1 + l3*l3 - l2*l2) / (2.0 * l1 * l3);
            currentAngleSum += std::acos(std::clamp(cosTheta_i, -1.0, 1.0));

            double cosTheta_j = (l1*l1 + l2*l2 - l3*l3) / (2.0 * l1 * l2);
            double cosTheta_k = (l3*l3 + l2*l2 - l1*l1) / (2.0 * l3 * l2);
                
            double theta_j = std::acos(std::clamp(cosTheta_j, -1.0, 1.0));
            double theta_k = std::acos(std::clamp(cosTheta_k, -1.0, 1.0));

            cotanSum += (1.0 / std::tan(theta_j)) + (1.0 / std::tan(theta_k));
        }

        if (!valid) {
            for (size_t k = 0; k < halfedges.size(); k++) {
                tri.edgeLengths[halfedges[k].edge()] = originalLens[k];
            }
            return false;
        }

        double error = (2 * PI) - currentAngleSum;
        if (std::abs(error) < TOLERANCE) {
            return true;
        }

        double step = error / (0.5 * cotanSum);
        ui -= step;
    }

    for (size_t k = 0; k < halfedges.size(); k++) {
        tri.edgeLengths[halfedges[k].edge()] = originalLens[k];
    }
    return false;
}

bool flatVertexRemoval(IntrinsicTriangulation& tri, Vertex i) {
    int candidates;
    while (i.degree() > 3) {
        candidates = i.degree();
        for (Edge e : i.adjacentEdges()) {
            if (tri.flipEdgeIfPossible(e)) {
                break;
            }
        }
        if (candidates == i.degree()) { break; } 
    }

    if (i.degree() == 3) {
        std::vector<Vertex> neighs;
        for (Vertex n : i.adjacentVertices()) {
            neighs.push_back(n);
        }

        if (neighs[0] == neighs[1] || neighs[1] == neighs[2] || neighs[0] == neighs[2]) {
            return false;
        }

        bool sharedFace = false;
        for (Face f : neighs[0].adjacentFaces()) {
            bool has1 = false, has2 = false, hasI = false;
            for (Vertex fv : f.adjacentVertices()) {
                if (fv == neighs[1]) has1 = true;
                if (fv == neighs[2]) has2 = true;
                if (fv == i) hasI = true; // Ignore the faces we are about to delete
            }
            if (has1 && has2 && !hasI) {
                sharedFace = true;
                break;
            }
        }

        if (sharedFace) {
            return false; 
        }

        tri.intrinsicMesh->removeVertex(i);
        tri.refreshQuantities();
        tri.flipToDelaunay();
        return true;
    } else {
        return false;
    }
}

void updateCurvatureChannels(IntrinsicTriangulation& tri, VertexData<double>& posKVal, VertexData<double>& negKVal) {
    // For each vertex, update the curvature value channels
    tri.requireVertexGaussianCurvatures();
    for (Vertex v : tri.mesh.vertices()) {
        posKVal[v] =  std::max(tri.vertexGaussianCurvatures[v], 0.0);
        negKVal[v] = -std::min(tri.vertexGaussianCurvatures[v], 0.0);
    }
}

void computeRemovalCost(IntrinsicTriangulation& tri, Vertex i, VertexData<double>& posKVal, VertexData<double>& negKVal, VertexData<Vector2>& posKTan, VertexData<Vector2>& negKTan, VertexData<double>& removalCost, bool updateTans, bool revertFlattening) {
    // Save the pre-flattening edge lengths
    std::unordered_map<Edge, double> storedEdgeLens;
    for (Edge e : i.adjacentEdges()) {
        storedEdgeLens[e] = tri.edgeLengths[e];
    }

    // Save the pre-flattening vertex curvatures
    std::unordered_map<Vertex, double> oldPosVertK;
    std::unordered_map<Vertex, double> oldNegVertK;
    double oldPosKTotal = 0.0;
    double oldNegKTotal = 0.0;
    for (Vertex j : i.adjacentVertices()) {
        oldPosVertK[j] =  std::max(0.0, tri.vertexGaussianCurvatures[j]);
        oldNegVertK[j] = -std::min(0.0, tri.vertexGaussianCurvatures[j]);
        oldPosKTotal +=  std::max(0.0, tri.vertexGaussianCurvatures[j]);
        oldNegKTotal += -std::min(0.0, tri.vertexGaussianCurvatures[j]);
    }
    oldPosVertK[i] =  std::max(0.0, tri.vertexGaussianCurvatures[i]);
    oldNegVertK[i] = -std::min(0.0, tri.vertexGaussianCurvatures[i]);
    
    // Attempt to flatten i via CETM
    if (!flattenVertexCETM(tri, i)) {
        // If flattening is unsuccessful, set max cost and return 
        removalCost[i] = std::numeric_limits<double>::max();
        return;
    }

    // Compute the transfer weight (\alpha_ij) for each neighbor vertex
    std::unordered_map<Vertex, double> posVertTransferWeight;
    std::unordered_map<Vertex, double> newPosVertK;
    std::unordered_map<Vertex, double> negVertTransferWeight;
    std::unordered_map<Vertex, double> newNegVertK;
    double newPosKTotal = 0.0;
    double newNegKTotal = 0.0;
    for (Vertex j : i.adjacentVertices()) {
        // We first have to find each new curvature and compute the new curvature sum over \cN_i
        auto k = tri.vertexGaussianCurvature(j);
        // pos
        newPosVertK[j] = std::max(k, 0.0);
        newPosKTotal += std::max(k, 0.0);
        // neg
        newNegVertK[j] = -std::min(k, 0.0);
        newNegKTotal += -std::min(k, 0.0);
    }

    // Compute cost
    double cost = 0.0;

    double posDiff = std::abs(newPosKTotal - oldPosKTotal);
    double negDiff = std::abs(newNegKTotal - oldNegKTotal);

    for (Halfedge he : i.outgoingHalfedges()) {
        // Compute values for finding transport cost from i to neighbor j
        Vertex j = he.tipVertex();
        auto posWeight = (posDiff > 1e-8) ? std::abs(newPosVertK[j]-oldPosVertK[j]) / posDiff : 0.0;
        auto negWeight = (negDiff > 1e-8) ? std::abs(newNegVertK[j]-oldNegVertK[j]) / negDiff : 0.0;

        auto edgeVect = tri.halfedgeVectorsInVertex[he.twin()];
        auto rot = tri.transportVectorsAlongHalfedge[he.twin()];

        // Compute the new tangent vectors
        double posDenom = posWeight * oldPosVertK[i] + oldPosVertK[j];
        Vector2 newPosTanVect = (posDenom > 1e-8) ? (posWeight * oldPosVertK[i] * (rot * posKTan[i] + edgeVect) + (oldPosVertK[j] * posKTan[j])) / posDenom : Vector2{0.0, 0.0};
        double negDenom = negWeight * oldNegVertK[i] + oldNegVertK[j];
        Vector2 newNegTanVect = (negDenom > 1e-8) ? (negWeight * oldNegVertK[i] * (rot * negKTan[i] + edgeVect) + (oldNegVertK[j] * negKTan[j])) / negDenom : Vector2{0.0, 0.0};

        // Compute cost and add to total
        auto posCost = newPosVertK[j] * newPosTanVect.norm();
        auto negCost = newNegVertK[j] * newNegTanVect.norm();
        cost += posCost + negCost;

        // If updating, apply the new tangent vectors
        if (updateTans) {
            posKTan[j] = newPosTanVect;
            negKTan[j] = newNegTanVect;
        }
    }
    removalCost[i] = cost;

    // If updating, update the curvature channels and return
    if (revertFlattening) {
        // Restore edge lengths and return
        for (Edge e : i.adjacentEdges()) {
            tri.edgeLengths[e] = storedEdgeLens[e];
        }
    } 
    else {
        tri.refreshQuantities();
        updateCurvatureChannels(tri, posKVal, negKVal);
    }
    
    return;
}

bool removeVertAndUpdate(IntrinsicTriangulation& tri, Vertex i, VertexData<double>& posKVal, VertexData<double>& negKVal, VertexData<Vector2>& posKTan, VertexData<Vector2>& negKTan, VertexData<double>& removalCost) {
    // Save edges in case some operation fails
    tri.requireEdgeLengths();
    std::unordered_map<Edge, double> storedEdgeLens;
    for (Edge e : i.adjacentEdges()) {
        storedEdgeLens[e] = tri.edgeLengths[e];
    }

    // Attempt local flattening
    if (!flattenVertexCETM(tri, i)) {
        return false;
    }
}


void intrinsicallyCoarsen(IntrinsicTriangulation& tri, size_t count) {
    // Get vertex count info for console output
    uint totalInitVs = tri.mesh.nVertices();

    // Flip to Delaunay
    tri.flipToDelaunay();

    // Initialize weight system values
    VertexData<double> removalCost(tri.mesh, 0.0);
    VertexData<double> posKVal(tri.mesh);
    VertexData<double> negKVal(tri.mesh);
    VertexData<Vector2> posTanVect(tri.mesh, Vector2{0.0,0.0});
    VertexData<Vector2> negTanVect(tri.mesh, Vector2{0.0,0.0});

    tri.requireVertexGaussianCurvatures();
    tri.requireEdgeLengths();
    tri.requireHalfedgeVectorsInVertex();
    tri.requireTransportVectorsAlongHalfedge();

    updateCurvatureChannels(tri, posKVal, negKVal);

    // Initialize a priority queue for vertex costs
    std::priority_queue<VertexCost, std::vector<VertexCost>, std::greater<VertexCost>> cost_pq; 

    // Initialize the cost value of each vertex, adding it to the queue as we go
    fmt::println("Initializing...");
    uint processed = 0;
    for (Vertex v : tri.mesh.vertices()) {
        computeRemovalCost(tri, v, posKVal, negKVal, posTanVect, negTanVect, removalCost, false, true);
        cost_pq.push(VertexCost{v, removalCost[v]});
        ++processed;
        if (processed % 1000 == 0) {
            fmt::println("Initialized: {} / {}", processed, totalInitVs);
        }
    }
    fmt::println("Initializing done.");

    // Main loop: Until target vertex count is reached or there are no more vertices to remove
    fmt::println("Coarsening to {} vertices.", int(count));
    processed = 0;
    while (tri.mesh.nVertices() > count && !cost_pq.empty()) {
        // Pop the top vertex from the queue
        Vertex i = cost_pq.top().v;
        double queuedCost = cost_pq.top().cost;
        cost_pq.pop();
        // std::cout << "BEGINNING REMOVAL OF VERTEX: " << i << std::endl;

        // Check if i is dead or stale
        if (i.isDead() || queuedCost != removalCost[i]) {
            // fmt::println("{} =/= {}", queuedCost, removalCost[i]);
            continue;
        }

        // Temporarily store the neighboring vertices of i
        std::vector<Vertex> neighbors;

        for (Vertex j : i.adjacentVertices()) {
            neighbors.push_back(j);
        }

        computeRemovalCost(tri, i, posKVal, negKVal, posTanVect, negTanVect, removalCost, true, false);
        if (removalCost[i] == std::numeric_limits<double>::max()) {
            continue;
        }

        // Remove i
        // std::cout << "ATTEMPTING TO REMOVE VERTEX: " << i << std::endl;
        if (!flatVertexRemoval(tri, i)) {
          removalCost[i] = std::numeric_limits<double>::max();
        //   fmt::println("Vertex could not be removed");
          continue;
        } else {
            ++processed;
        }
        // std::cout << "REMOVED VERTEX: " << i << std::endl;

        // Flip remaining edges to intrinsic Delaunay as needed
        // TODO: Don't flip all edges, just those remaining from the neighborhood
        // std::cout << "FLIPPING BACK TO DELAUNAY" << std::endl;
        tri.flipToDelaunay();

        // Update the cost of the neighboring vertices
        // std::cout << "BEGINNING UPDATE NEIGHBORS" << std::endl;
        for (Vertex j : neighbors) {
            if (!j.isDead() && removalCost[j] != std::numeric_limits<double>::max()) {
                // std::cout << "ATTEMPTING TO UPDATE NEIGHBOR: " << j << std::endl;
                computeRemovalCost(tri, j, posKVal, negKVal, posTanVect, negTanVect, removalCost, false, true);
                cost_pq.push(VertexCost{j, removalCost[j]});
            }
        }

        
        if (processed % 100 == 0) {
            fmt::println("Vertices removed: {} / {}", processed, totalInitVs-count);
        }
    }
    tri.mesh.compress();
    fmt::println("Simplification done.");
}
