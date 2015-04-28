module solid_boundary_interface_effect;

import geom;
import globaldata;
import solidfvinterface;

class SolidBoundaryInterfaceEffect {
public:
    int blkId;
    int whichBoundary;
    
    this(int id, int boundary) {
	blkId = id;
	whichBoundary = boundary;
    }

    void apply(double t, int tLevel) {}
}

class FixedTInterfaceEffect : SolidBoundaryInterfaceEffect {
public:
    this(int id, int boundary, double Twall)
    {
	super(id, boundary);
	_Twall = Twall;
    }

    override void apply(double t, int tLevel)
    {
	size_t i, j, k;
	SolidFVInterface IFace;
	auto blk = allSolidBlocks[blkId];

	final switch (whichBoundary) {
	case Face.north:
	    j = blk.jmax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.east:
	    i = blk.imax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.top:
	    throw new Error("[TODO] FixedT bc not implemented for TOP face.");
	case Face.bottom:
	    throw new Error("[TODO] FixedT bc not implemented for BOTTOM face.");

	}

    }

private:
    double _Twall;
}
