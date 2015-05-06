/**
 * solidblock.d
 *
 * Base class for a block representing a solid.
 * Typically, we want to compute the heat transfer
 * through the solid and its effect on the adjoining
 * flow field.
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-22-04
 */

module solidblock;

import std.conv;
import luad.all;
import geom;
import solidfvcell;
import solidbc;
import solidprops;


class SolidBlock {
public:
    int id; // block identifier
    string label;
    bool active; // if true, block participates in time integration
    double energyResidual; // monitor this for steady state
    Vector3 energyResidualLoc; // location of worst case
    int hncell; // number of history cells
    SolidProps sp; // material properties for whole block
    LuaState udfSourceTerms;

    SolidFVCell[] activeCells; // collection of references to active cells in the domain
    SolidBoundaryCondition[] bc; // collection of references to boundary conditions

    override string toString() const { return "SolidBlock(id=" ~ to!string(id) ~ ")"; }
    abstract void assembleArrays();
    abstract void bindFacesAndVerticesToCells();
    abstract void readGrid(string filename);
    abstract void writeGrid(string filename, double sim_time);
    abstract void readSolution(string filename);
    abstract void writeSolution(string fileName, double simTime);
    abstract void computePrimaryCellGeometricData();
    abstract void computeSecondaryCellGeometricData();

    abstract void applyPreSpatialDerivAction(double t, int tLevel);
    abstract void applyPostFluxAction(double t, int tLevel);
    abstract void computeSpatialDerivatives(int ftl);
    abstract void computeFluxes();
}

