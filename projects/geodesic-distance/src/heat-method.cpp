// Implement member functions HeatMethod class.
#include "heat-method.h"
#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    // TODO: Build Laplace and flow matrices.
    // Note: core/geometry.cpp has meanEdgeLength() function
    this->A = geometry->laplaceMatrix();

    double h = geometry->meanEdgeLength();
    this->F = geometry->massMatrix() + h * h * A;
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {

    FaceData<Vector3> X(*mesh, {0, 0, 0});
    for (Face f : mesh->faces()) {
        double area = geometry->faceArea(f);
        Vector3 N = geometry->faceNormal(f);
        Vector3 Du {0, 0, 0};

        Halfedge st = f.halfedge();
        Halfedge he = st;
        do {
            size_t idxV = he.next().tipVertex().getIndex();
            double ui = u[idxV];
            Vector3 ei = geometry->halfedgeVector(he);
            Du += ui * cross(N, ei);

            he = he.next();
        } while (he != st);

        Du /= 2 * area;
        X[f] = -Du / Du.norm();
    }

    return X;
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {

    Vector<double> divX = Vector<double>::Zero(mesh->nVertices());
    for (Vertex v : mesh->vertices()) {
        Halfedge st = v.halfedge();
        Halfedge he = st;

        double divXv = 0.0;

        do {
            Face f = he.face();
            if (!f.isBoundaryLoop()) {
                Vector3 Xj = X[f];
                Vector3 e1 = geometry->halfedgeVector(he);
                Vector3 e2 = geometry->halfedgeVector(he.next().next().twin());
                divXv += geometry->cotan(he) * dot(e1, Xj) + geometry->cotan(he.next().next()) * dot(e2, Xj);
            }
            he = he.twin().next();
        } while (he != st);

        divXv /= 2;
        divX[v.getIndex()] = divXv;
    }
    return divX;
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    SparseMatrix<double> _A = A;
    SparseMatrix<double> _F = F;

    Vector<double> divX = -computeDivergence(computeVectorField(solvePositiveDefinite(_F, delta)));
    Vector<double> phi = solvePositiveDefinite(_A, divX);

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}
