// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {

    SparseMatrix<std::complex<double>> ED = geometry->complexLaplaceMatrix() / 2.0;

    typedef Eigen::Triplet<std::complex<double>> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for(BoundaryLoop bl : mesh->boundaryLoops()) {
        for(Halfedge he : bl.adjacentHalfedges()) {
            size_t idxV1 = he.tailVertex().getIndex();
            size_t idxV2 = he.tipVertex().getIndex();

            tripletList.push_back(TripletEntry(idxV1, idxV2, std::complex<double>(0.0,  .25)));
            tripletList.push_back(TripletEntry(idxV2, idxV1, std::complex<double>(0.0, -.25)));
        }
    }

    // Construct sparse matrix from triplets
    SparseMatrix<std::complex<double>> Area(mesh->nVertices(), mesh->nVertices());
    Area.setFromTriplets(tripletList.begin(), tripletList.end());

    return ED - Area;
}


/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {

    // TODO
    Vector<std::complex<double>> eigen = solveInversePowerMethod(buildConformalEnergy());

    VertexData<Vector2> dict(*mesh);

    for (Vertex v : mesh->vertices()) {
        size_t idxV = v.getIndex();
        Vector2 vec = Vector2::fromComplex(eigen[idxV]);
        dict[v] = vec;
    }

    return dict;
}
