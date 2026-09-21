#pragma once
#include "SimplicialDual3D.hpp"
#include "common/LinearSolveOptions.hpp"

namespace simplicial3d {
// Section 4.5: the missing wall circulation is Omega - C U. This small
// constrained reconstruction makes that integral participate in interpolation.
// The paper does not specify the underdetermined boundary-vector solve; we pick
// the minimum change to the supplied tangential trace, subject to its integrals.
class BoundaryCirculation {
public:
    BoundaryCirculation(const Mesh& mesh,const DualMesh& dual);
    static Point tangent(Point p,unsigned wallMask);
    LinearSolveResult reconstruct(const DualMesh& dual,std::span<const double> missing,
                                 std::vector<Point>& velocity,int maxIterations) const;
private:
    struct Coefficient {int vertex;Point value;};
    struct Row {int primalEdge;std::vector<Coefficient> coefficients;double diagonal=0;};
    std::vector<Row> rows_;
};
}
