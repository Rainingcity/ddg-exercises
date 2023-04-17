// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

    if (!he.isInterior()) return 0.0;

    Vector3 v = halfedgeVector(he.next().next());
    Vector3 w = halfedgeVector(he.next().twin());
    return dot(v, w) / norm(cross(v, w));
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;

    double totalArea = 0.0;

    do {
        if (!he.face().isBoundaryLoop())
            totalArea += faceArea(he.face());
        he = he.twin().next();
    } while (he != st);

    return totalArea / 3.0;
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    Halfedge he = c.halfedge();
    while(he.next() != c.halfedge()) he = he.next();
    Vector3 v = halfedgeVector(c.halfedge());
    Vector3 w = halfedgeVector(he.twin());
    double cosine = dot(v, w) / norm(v) / norm(w);
    return acos(cosine);
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    Vector3 v = faceNormal(he.face());
    Vector3 w = faceNormal(he.twin().face());
    Vector3 e = halfedgeVector(he);
    return atan2(dot(e / e.norm(), cross(v, w)), dot(v, w));
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    int count = 0;
    Vector3 totalNormal = {0, 0, 0};
    do {
        if (!he.face().isBoundaryLoop()) {
            if (he.edge().halfedge() == he)
                totalNormal += faceNormal(he.face());
            else
                totalNormal -= faceNormal(he.face());
            count++;
        }

        he = he.twin().next();
    } while (he != st);
    return totalNormal / count;
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    Vector3 totalNormal = {0, 0, 0};
    do {
        if (!he.face().isBoundaryLoop()) {
            Corner c = he.corner();
            double currAngle = angle(c);
            Vector3 currNormal = faceNormal(he.face());
            totalNormal += currAngle * currNormal;
        }

        he = he.twin().next();
    } while (he != st);
    return totalNormal / totalNormal.norm();
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    Vector3 totalNormal = {0, 0, 0};
    do {
        if (!he.twin().face().isBoundaryLoop()) {
            Vector3 v = halfedgeVector(he);
            Vector3 w = halfedgeVector(he.twin().next());
            totalNormal += cross(w, v) / v.norm2() / w.norm2();
        }

        he = he.twin().next();
    } while (he != st);
    return totalNormal / totalNormal.norm();
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    Vector3 totalNormal = {0, 0, 0};
    do {
        Face f = he.face();
        if (!f.isBoundaryLoop()) {
            totalNormal += faceArea(f) * faceNormal(f);
        }

        he = he.twin().next();
    } while (he != st);
    return totalNormal / totalNormal.norm();
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    Vector3 totalNormal = {0, 0, 0};
    do {
        if (!he.twin().face().isBoundaryLoop() && !he.face().isBoundaryLoop()) {
            double theta = dihedralAngle(he);
            Vector3 e = halfedgeVector(he);
            totalNormal += theta * e / e.norm();
        }

        he = he.twin().next();
    } while (he != st);
    return totalNormal / totalNormal.norm();
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    Vector3 totalNormal = {0, 0, 0};
    do {
        if (!he.twin().face().isBoundaryLoop()) {
            Vector3 v = halfedgeVector(he);
            Vector3 w = halfedgeVector(he.twin().next());
            totalNormal += cotan(he.twin()) * v + cotan(he.twin().next()) * w;
        }

        he = he.twin().next();
    } while (he != st);
    return totalNormal / totalNormal.norm();
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    double defect = 2 * M_PI;
    do {
        if (!he.face().isBoundaryLoop()) {
            Corner c = he.corner();
            defect -= angle(c);
        }

        he = he.twin().next();
    } while (he != st);
    return defect;
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    double totalDefect = 0.0;
    for(Vertex v : mesh.vertices()) {
        totalDefect += angleDefect(v);
    }

    return totalDefect;
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    double meanCurvature = 0.0;
    do {
        if (!he.twin().face().isBoundaryLoop() && !he.face().isBoundaryLoop()) {
            double theta = dihedralAngle(he);
            Vector3 e = halfedgeVector(he);
            meanCurvature += theta * e.norm();
        }

        he = he.twin().next();
    } while (he != st);
    return meanCurvature / 2;
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    Halfedge st = v.halfedge();
    Halfedge he = st;
    double totalArea = 0.0;
    do {
        if (!he.twin().face().isBoundaryLoop()) {
            Vector3 v = halfedgeVector(he.twin());
            Vector3 w = halfedgeVector(he.twin().next());
            totalArea += cotan(he.twin()) * v.norm2() + cotan(he.twin().next()) * w.norm2();
        }

        he = he.twin().next();
    } while (he != st);
    return totalArea / 8;
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    double area = circumcentricDualArea(v);
    double K = angleDefect(v) / area;
    double H = scalarMeanCurvature(v) / area;
    double diff = sqrt(H * H - K);
    return std::make_pair(H - diff, H + diff);
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    typedef Eigen::Triplet<double> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for(Vertex v : mesh.vertices()) {
        size_t idxV = v.getIndex();
        Halfedge st = v.halfedge();
        Halfedge he = st;
        double sum = 0.0;
        do {
            double curr = 0.0;
            if (!he.face().isBoundaryLoop()) curr += cotan(he);
            if (!he.twin().face().isBoundaryLoop()) curr += cotan(he.twin());
            curr /= 2.0;
            sum += curr;

            tripletList.push_back(TripletEntry(idxV, he.tipVertex().getIndex(), -curr));

            he = he.twin().next();
        } while (he != st);

        tripletList.push_back(TripletEntry(idxV, idxV, sum + 1e-8));
    }

    // Construct sparse matrix from triplets
    SparseMatrix<double> newMatrix(mesh.nVertices(), mesh.nVertices());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    typedef Eigen::Triplet<double> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for(Vertex v : mesh.vertices()) {
        size_t idxV = v.getIndex();
        tripletList.push_back(TripletEntry(idxV, idxV, barycentricDualArea(v)));
    }

    // Construct sparse matrix from triplets
    SparseMatrix<double> newMatrix(mesh.nVertices(), mesh.nVertices());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    typedef Eigen::Triplet<std::complex<double>> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for(Vertex v : mesh.vertices()) {
        size_t idxV = v.getIndex();
        Halfedge st = v.halfedge();
        Halfedge he = st;
        double sum = 0.0;
        do {
            double curr = (cotan(he) + cotan(he.twin())) / 2.0;
            sum += curr;

            tripletList.push_back(TripletEntry(idxV, he.tipVertex().getIndex(), std::complex<double>(-curr, 0)));

            he = he.twin().next();
        } while (he != st);

        tripletList.push_back(TripletEntry(idxV, idxV, std::complex<double>(sum + 1e-8, 0)));
    }

    // Construct sparse matrix from triplets
    SparseMatrix<std::complex<double>> newMatrix(mesh.nVertices(), mesh.nVertices());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral
