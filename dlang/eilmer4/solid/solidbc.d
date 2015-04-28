/**
 * solidbc.d
 *
 * Author: Rowan G. and Peter J.
 */

module solidbc;

import solid_boundary_interface_effect;
import solid_boundary_flux_effect;

class SolidBoundaryCondition {
public:
    int blkId;
    int whichBoundary;
    bool setsFluxDirectly;

    final void applyPreSpatialDerivAction(double t, int tLevel)
    {
	foreach ( sie; preSpatialDerivAction ) sie.apply(t, tLevel);
    }

    final void applyPostFluxAction(double t, int tLevel)
    {
	foreach ( sfe; postFluxAction ) sfe.apply(t, tLevel);
    }

private:
    SolidBoundaryInterfaceEffect[] preSpatialDerivAction;
    SolidBoundaryFluxEffect[] postFluxAction;

}
