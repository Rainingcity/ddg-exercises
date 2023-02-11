// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.

    typedef Eigen::Triplet<size_t> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for (Edge e : mesh->edges()) {
        // Get Edge index
        size_t idxE = e.getIndex();

        // Get Vertex index
        Vertex v1 = e.halfedge().vertex();
        size_t idxV1 = v1.getIndex();
        Vertex v2 = e.halfedge().twin().vertex();
        size_t idxV2 = v2.getIndex();

        // Add adjacent relationship
        tripletList.push_back(TripletEntry(idxE, idxV1, 1));
        tripletList.push_back(TripletEntry(idxE, idxV2, 1));
    }

    // Construct sparse matrix from triplets
    SparseMatrix<size_t> newMatrix(mesh->nEdges(), mesh->nVertices());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    typedef Eigen::Triplet<size_t> TripletEntry;
    std::vector<TripletEntry> tripletList{};

    for (Face f : mesh->faces()) {
        size_t idxF = f.getIndex();

        Halfedge start = f.halfedge();
        Halfedge he = start;

        do {
            Edge e = he.edge();
            size_t idxE = e.getIndex();

            tripletList.push_back(TripletEntry(idxF, idxE, 1));

            he = he.next();
        } while (he != start);
    }

    // Construct sparse matrix from triplets
    SparseMatrix<size_t> newMatrix(mesh->nFaces(), mesh->nEdges());
    newMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return newMatrix;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> subsetVertices = Vector<size_t>::Zero(mesh->nVertices());
    for (size_t vertex : subset.vertices) {
        subsetVertices[vertex] = 1;
    }
    return subsetVertices;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> subsetEdges = Vector<size_t>::Zero(mesh->nEdges());
    for (size_t edge : subset.edges) {
        subsetEdges[edge] = 1;
    }
    return subsetEdges;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> subsetFaces = Vector<size_t>::Zero(mesh->nFaces());
    for (size_t face : subset.faces) {
        subsetFaces[face] = 1;
    }
    return subsetFaces;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    MeshSubset star = subset.deepCopy();

    for (size_t vertex : subset.vertices) {
        Vertex v = mesh->vertex(vertex);
        Halfedge start = v.halfedge();
        Halfedge he = start;
        do {
            Edge e = he.edge();
            Face f = he.face();

            star.addEdge(e.getIndex());
            if (!f.isBoundaryLoop())
                star.addFace(f.getIndex());

            he = he.twin().next();
        } while (he != start);
    }

    for (size_t edge : subset.edges) {
        Edge e = mesh->edge(edge);

        if (!e.halfedge().face().isBoundaryLoop())
            star.addFace(e.halfedge().face().getIndex());

        if (!e.halfedge().twin().face().isBoundaryLoop())
            star.addFace(e.halfedge().twin().face().getIndex());
    }

    return star;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset closure = subset.deepCopy();

    for (size_t face : subset.faces) {
        Face f = mesh->face(face);
        Halfedge start = f.halfedge();
        Halfedge he = start;
        do {
            Edge e = he.edge();
            Vertex v = he.vertex();

            closure.addVertex(v.getIndex());
            closure.addEdge(e.getIndex());

            he = he.next();
        } while (he != start);
    }

    for (size_t edge : subset.edges) {
        Edge e = mesh->edge(edge);
        closure.addVertex(e.halfedge().vertex().getIndex());
        closure.addVertex(e.halfedge().twin().vertex().getIndex());
    }

    return closure;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    MeshSubset closure = SimplicialComplexOperators::closure(subset);
    MeshSubset stcl = SimplicialComplexOperators::star(closure);

    MeshSubset star = SimplicialComplexOperators::star(subset);
    MeshSubset clst = SimplicialComplexOperators::closure(star);

    clst.deleteSubset(stcl);

    return clst;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    MeshSubset closure = SimplicialComplexOperators::closure(subset);
    closure.deleteSubset(subset);
    return closure.vertices.empty() && closure.edges.empty() && closure.faces.empty();
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    if (!SimplicialComplexOperators::isComplex(subset)) return -1;
    MeshSubset pure;
    int deg = -1;

    if (!subset.faces.empty()) {
        pure.addFaces(subset.faces);
        deg = 2;
    } else if (!subset.edges.empty()) {
        pure.addEdges(subset.edges);
        deg = 1;
    } else if (!subset.vertices.empty()) {
        pure.addVertices(subset.vertices);
        deg = 0;
    } else {
        deg = 0;
    }
    MeshSubset cl_pure = SimplicialComplexOperators::closure(pure);

    MeshSubset copy;
    copy.addSubset(subset);
    copy.deleteSubset(cl_pure);

    return (copy.vertices.empty() && copy.edges.empty() && copy.faces.empty()) ? deg : -1;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    MeshSubset boundary;
    int deg = SimplicialComplexOperators::isPureComplex(subset);
    std::set<size_t> garbage;
    if (deg == 2) {
        for (size_t face : subset.faces) {
            Face f = mesh->face(face);

            Halfedge start = f.halfedge();
            Halfedge he = start;

            do {
                size_t idxE = he.edge().getIndex();

                if (garbage.find(idxE) == garbage.end()) {
                    if (boundary.edges.find(idxE) == boundary.edges.end()) {
                        boundary.addEdge(idxE);
                    } else {
                        boundary.deleteEdge(idxE);
                        garbage.insert(idxE);
                    }
                }

                he = he.next();
            } while (he != start);
        }
    } else if (deg == 1) {
        for (size_t edge : subset.edges) {
            Edge e = mesh->edge(edge);
            size_t idxV1 = e.halfedge().vertex().getIndex();
            size_t idxV2 = e.halfedge().twin().vertex().getIndex();

            if (garbage.find(idxV1) == garbage.end()) {
                if (boundary.vertices.find(idxV1) == boundary.vertices.end()) {
                    boundary.addVertex(idxV1);
                } else {
                    boundary.deleteVertex(idxV1);
                    garbage.insert(idxV1);
                }
            }

            if (garbage.find(idxV2) == garbage.end()) {
                if (boundary.vertices.find(idxV2) == boundary.vertices.end()) {
                    boundary.addVertex(idxV2);
                } else {
                    boundary.deleteVertex(idxV2);
                    garbage.insert(idxV2);
                }
            }

        }
    }

    return SimplicialComplexOperators::closure(boundary);
}
