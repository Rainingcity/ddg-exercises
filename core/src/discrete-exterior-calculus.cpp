// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    typedef Eigen::Triplet<double> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for (Vertex v : mesh.vertices()) {
        size_t idxV = v.getIndex();
        tripletList.push_back(TripletEntry(idxV, idxV, barycentricDualArea(v)));
    }

    // Construct sparse matrix from triplets
    SparseMatrix<double> newMatrix(mesh.nVertices(), mesh.nVertices());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    typedef Eigen::Triplet<double> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for (Edge e : mesh.edges()) {
        size_t idxE = e.getIndex();

        double sumCot = 0.0;
        if (!e.halfedge().face().isBoundaryLoop()) sumCot += cotan(e.halfedge());
        if (!e.halfedge().twin().face().isBoundaryLoop()) sumCot += cotan(e.halfedge().twin());

        tripletList.push_back(TripletEntry(idxE, idxE, sumCot / 2));
    }

    // Construct sparse matrix from triplets
    SparseMatrix<double> newMatrix(mesh.nEdges(), mesh.nEdges());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    typedef Eigen::Triplet<double> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for (Face f : mesh.faces()) {
        size_t idxF = f.getIndex();

        double s = 0.0;
        Halfedge st = f.halfedge();

        Halfedge he = st;
        do {
            s += edgeLength(he.edge());
            he = he.next();
        } while (he != st);
        s = s / 2;

        double area = s;
        he = st;
        do {
            area *= (s - edgeLength(he.edge()));
            he = he.next();
        } while (he != st);

        tripletList.push_back(TripletEntry(idxF, idxF, 1 / sqrt(area)));
    }

    // Construct sparse matrix from triplets
    SparseMatrix<double> newMatrix(mesh.nFaces(), mesh.nFaces());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    typedef Eigen::Triplet<double> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for (Edge e : mesh.edges()) {
        size_t idxE = e.getIndex();

        Vertex v1 = e.halfedge().vertex();
        size_t idxV1 = v1.getIndex();
        Vertex v2 = e.halfedge().twin().vertex();
        size_t idxV2 = v2.getIndex();

        tripletList.push_back(TripletEntry(idxE, idxV1, -1.0));
        tripletList.push_back(TripletEntry(idxE, idxV2, 1.0));
    }

    // Construct sparse matrix from triplets
    SparseMatrix<double> newMatrix(mesh.nEdges(), mesh.nVertices());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    typedef Eigen::Triplet<double> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for (Edge e : mesh.edges()) {
        size_t idxE = e.getIndex();

        Face f = e.halfedge().face();
        if (!f.isBoundaryLoop()) {
            size_t idxF = f.getIndex();
            tripletList.push_back(TripletEntry(idxF, idxE, 1.0));
        }

        f = e.halfedge().twin().face();
        if (!f.isBoundaryLoop()) {
            size_t idxF = f.getIndex();
            tripletList.push_back(TripletEntry(idxF, idxE, -1.0));
        }
    }

    // Construct sparse matrix from triplets
    SparseMatrix<double> newMatrix(mesh.nFaces(), mesh.nEdges());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

} // namespace surface
} // namespace geometrycentral
