// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    return M + h * geometry->laplaceMatrix();
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // Note: Update positions via geometry->inputVertexPositions
    SparseMatrix<double> Delta = buildFlowOperator(geometry->massMatrix(), h);
    SparseMatrix<double> M = geometry->massMatrix();
    Vector<double> x_rho(mesh->nVertices()), y_rho(mesh->nVertices()), z_rho(mesh->nVertices());
    for (Vertex v : mesh->vertices()) {
        size_t idxV = v.getIndex();
        x_rho[idxV] = geometry->inputVertexPositions[v].x;
        y_rho[idxV] = geometry->inputVertexPositions[v].y;
        z_rho[idxV] = geometry->inputVertexPositions[v].z;
    }

    x_rho = M * x_rho;
    y_rho = M * y_rho;
    z_rho = M * z_rho;

    Vector<double> xh_rho = solvePositiveDefinite(Delta, x_rho);
    Vector<double> yh_rho = solvePositiveDefinite(Delta, y_rho);
    Vector<double> zh_rho = solvePositiveDefinite(Delta, z_rho);

    for (Vertex v : mesh->vertices()) {
        size_t idxV = v.getIndex();
        geometry->inputVertexPositions[v].x = xh_rho[idxV];
        geometry->inputVertexPositions[v].y = yh_rho[idxV];
        geometry->inputVertexPositions[v].z = zh_rho[idxV];
    }
}
