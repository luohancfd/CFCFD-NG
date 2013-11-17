/* \file baldwin_lomax.cxx
 * \ingroup eilmer3
 * \brief Baldwin-Lomax turbulence model for cns4u.c.
 *
 * \author PA Jacobs  ICASE 1991, UQ 2008
 * \author Richard Goozee UQ, 2002
 *
 * \version 16-Jan-07 : Changed the search for the maximum of F 
 *                      to the absolute maximum 
 * \version 26-Jun-97 : add iywall variable 
 *                      correction in step 2 iy --> iywall
 *                      speed up step 9 by doing sqrt only once
 *                      check for zero F_B in step 10
 * \version 09-Jul-97 : Compressibility factor for VanDriest damping
 * \version 12-Jul-97 : Improve robustness and update the data structures.
 * \version 10-Sep-97 : compute distance from wall more generally
 * \version 16-Oct-97 : Generalise Prandtl number and Cp
 * \version 20-Oct-97 : Fixed the logic for computing the cross-over point.
 * \version 19-Dec-97 : Search for the first maximum of F from the wall
 *                      (as suggested by Andrew M. and Klaus H.)
 * \version 30-Sep-2008: ported to Elmer3 for Rainer Kirchhartz
 * \version 28-Oct-2008: The Baldwin-Lomax model now allows turbulent 
 *                       viscous calculations on SOUTH boundaries as well. 
 *
 */

/*-----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "visc.hh"
#include "diffusion.hh"
#include "baldwin_lomax.hh"

//-------------------------------------------------------------------------
// Coefficients for the Baldwin-Lomax Eddy Viscosity Model 

// Overall characteristics of the boundary layer model.
const int COMPRESSIBILITY = 1;
const int TRANSITIONAL = 1;

// Parameters that are not often varied.
const double C_mutm = 14.0;
const double C_WK = 1.0;
const double K = 0.0168;

// A_plus should be varied for boundary layers in pressure gradients.
const double A_plus = 26.0;

// Baldwin and Lomax (1978) parameters for the flat plate boundary layer.
// Used successfully in Peter's models of a Mach 8 nozzle using air but
// worked badly in Richard Goozee's simulations of the entire Drummond tunnel. 
//
// See the turbulence model code in mbcns2 for more detail.
const double C_CP = 1.6;
const double C_KLEB = 0.3;
const double k_inner = 0.4;


//---------------------------------------------------------------------------
/// \brief Apply the Baldwin-Lomax turbulence model by computing the
///        turbulent transport coefficients.
///
/// \verbatim
/// On entry, this routine assumes that the laminar coefficients
/// (at the cell interface midpoints) and the velocity derivatives 
/// (at the cell vertices) have been previously computed.
///
/// The geometric assumptions are that the North boundary is a 
/// no-slip boundary, is located above the South boundary and
/// is aligned with the x-axis.  Of course, this is not true
/// generally but we hope that it is not to far off.
///
/// The coefficients (some modifiable and some not) are stored in cns_turb.hh
///
/// Reference:
/// ---------
/// B. S. Baldwin and H. Lomax
/// Thin layer approximation and algebraic model for separated
/// turbulent flows.
/// AIAA Paper 78-257
///
/// Application of this turbulence model is problem dependent
/// so this section of code must be altered to suit the problem
/// at hand.  Originally, we have coded it to compute the 
/// Mach 8 nozzle problem as defined in ICASE Report 91-01.
///
/// We assume that the viscous wall is the North boundary 
/// and that the j index moves us normal to the wall. 
/// We apply the model one i-station at a time.
/// \endverbatim
///
int baldwin_lomax_turbulence_model(global_data& gdata, Block& blk, size_t gtl)
{
    int i, j, j_limit, jwall;
    int ifirst, ilast, jfirst, jlast;
    int j_crossover;
    int keep_turbulent_viscosity;

    FV_Cell *cell;
    vector<double> vort, eta, y_plus, D, F, F_KLEB;
    vector<double> mu_t_inner, mu_t_outer, mu_t;
    double tau_max, mu_wall, rho_wall, y_wall, x_wall, dx, dy;
    double factor, temporary;
    double F_MAX, eta_MAX, U_MAX, U_MIN, U_DIFF, F_A, F_B, F_WAKE;
    double Pr_t = gdata.turbulence_prandtl;

    // Dimension of the work arrays, allowing for ghost cells.
    int ndim = blk.nnj + 4;
    vort.resize(ndim);
    eta.resize(ndim);
    y_plus.resize(ndim);
    D.resize(ndim);
    F.resize(ndim);
    F_KLEB.resize(ndim);
    mu_t_inner.resize(ndim);
    mu_t_outer.resize(ndim);
    mu_t.resize(ndim);

    Gas_model *gm = get_gas_model_ptr();
    int status_flag;

    if ( gdata.dimensions != 2 ) {
	cout << "baldwin_lomax_turbulence_model(): only works for 2D flow"
	     << " with the NORTH face being a solid wall." << endl;
	exit( NOT_IMPLEMENTED_ERROR );
    }


    // If NORTH boundary is the wall ..
    if ( blk.bcp[NORTH]->is_wall() && blk.bcp[NORTH]->type_code != SLIP_WALL ) {

        // Index limits.
        ifirst = blk.imin;
        ilast = blk.imax;
        jfirst = blk.jmax; // Start at the wall ...
        jwall = blk.jmax;
        jlast = blk.jmin;  // ...and move to the centre of the tube.
        // Set the search limit for determining inner/outer layer.
        j_limit = blk.jmax - (blk.nnj / 4);

        for ( i = ifirst; i <= ilast; ++i ) {
            /*
             * Before proceeding with the turbulence calculation,
             * check to see if the gas is moving with reasonable speed.
             * Use U_MAX later in step 9. Note that U_MAX is presently the
             * velocity squared.
             */
            U_MAX = 0.0;
            for ( j = jfirst; j >= jlast; --j ) {
                cell = blk.get_cell(i,j);
                temporary = cell->fs->vel.x * cell->fs->vel.x + cell->fs->vel.y * cell->fs->vel.y;
                U_MAX = max(U_MAX, temporary);
            }
            if (U_MAX < 10.0) {
	        for ( j = jfirst; j >= jlast; --j ) {
	            cell = blk.get_cell(i,j);
		    cell->fs->mu_t = 0.0;
		    cell->fs->k_t = 0.0;
	        }
	        continue; // skip to the next i-position.
	    }

            // Step 1: Compute the vorticity magnitude at the vertical interface
            //         midpoints. Compute the normal coordinate assuming that 
            //         the j index direction is aligned with the wall-normal direction.
            //
	    // Vertex indices as shown here.
	    // 3------2
	    // |      |
	    // | cell |
	    // |      |
	    // 0------1
	    cell = blk.get_cell(i,jwall);
            y_wall = 0.5 * (cell->vtx[3]->pos[gtl].y + cell->vtx[2]->pos[gtl].y);
            x_wall = 0.5 * (cell->vtx[3]->pos[gtl].x + cell->vtx[2]->pos[gtl].x);
            for ( j = jwall; j >= jlast; --j ) {
	        cell = blk.get_cell(i,j);
                temporary = 0.25*(cell->vtx[0]->dudy + cell->vtx[1]->dudy + 
		    	          cell->vtx[2]->dudy + cell->vtx[3]->dudy)
		    - 0.25*(cell->vtx[0]->dvdx + cell->vtx[1]->dvdx + 
			    cell->vtx[2]->dvdx + cell->vtx[3]->dvdx);
                vort[j] = fabs(temporary);
                dy = y_wall - cell->pos[gtl].y;
                dx = x_wall - cell->pos[gtl].x;
                eta[j] = sqrt(dx * dx + dy * dy);
            }    

            // Step 2: Compute the maximum shear stress --
            //         we will assume that it occurs at the wall.
	    cell = blk.get_cell(i,jwall);
            mu_wall = cell->fs->gas->mu;
            tau_max = mu_wall * vort[jwall];
            rho_wall = cell->fs->gas->rho;

            // Step 3: Compute y-plus across the domain.
            temporary = sqrt(rho_wall * tau_max) / mu_wall;
            for ( j = jfirst; j >= jlast; --j )
                y_plus[j] = temporary * eta[j];

            // Step 4: Compute the vanDriest damping factor
            for ( j = jfirst; j >= jlast; --j ) {
#               if COMPRESSIBILITY == 1
	        cell = blk.get_cell(i,j);
                factor = sqrt(cell->fs->gas->rho / rho_wall) * mu_wall / cell->fs->gas->mu;
#               else
                factor = 1.0;
#               endif
                D[j] = 1.0 - exp(-y_plus[j] * factor / A_plus);
            }

            // Step 5: Compute the inner-layer turbulent viscosity.
            for ( j = jfirst; j >= jlast; --j ) {
                temporary = k_inner * eta[j] * D[j];
                mu_t_inner[j] = blk.get_cell(i,j)->fs->gas->rho * temporary * temporary * vort[j];
            }

            // Step 6: Start computing the pieces for the outer-layer viscosity.
	    //         Compute F.
            for ( j = jfirst; j >= jlast; --j ) {
                F[j] = eta[j] * vort[j] * D[j];
            }

            // Step 7: Find the maximum of F and the point at which it occurs.
            F_MAX = F[jfirst];
            eta_MAX = eta[jfirst];
            for ( j = jfirst - 1; j >= jlast; --j ) {
                if ( F[j] > F_MAX ) {
                    F_MAX = F[j];
                    eta_MAX = eta[j];
	        }
            } // end for j
            // If F_MAX is zero, there isn't anything happening at this i station; 
            // skip to the next i station.
            if ( F_MAX < 1.0e-10 ) {
	        for ( j = jfirst; j >= jlast; --j ) {
		    cell = blk.get_cell(i,j);
		    cell->fs->mu_t = 0.0;
		    cell->fs->k_t = 0.0;
	        }
                continue; // ...with next i-station
	    }

            // Step 8: Compute the Klebanoff Intermittency factor.
            for ( j = jfirst; j >= jlast; --j ) {
                temporary = C_KLEB * eta[j] / eta_MAX;
                F_KLEB[j] = 1.0 / (1.0 + 5.5 * pow(temporary,6) );
            }

            // Step 9: Find the maximum and minimum total velocities.
            //         Currently set U_MIN = 0.
            U_MIN = 0.0;
            U_MAX = sqrt(U_MAX); // square of total velocity was computed before
            U_DIFF = U_MAX - U_MIN;

            // Step 10: Compute F_WAKE.
            F_A = eta_MAX * F_MAX;
            if ( F_MAX > 1.0e-10 )
                F_B = C_WK * eta_MAX * U_DIFF * U_DIFF / F_MAX;
            else
                F_B = F_A;
            F_WAKE = min(F_B, F_A);

            // Step 11: Compute the outer-layer turbulent viscosity.
            for ( j = jfirst; j >= jlast; --j ) {
                mu_t_outer[j] = K * C_CP * blk.get_cell(i,j)->fs->gas->rho * F_WAKE * F_KLEB[j];
            }

            // Step 12: Find the cross-over point for the inner- and
            //          outer-layers.  Assume that, near the wall, mu_t_inner
            //          is smaller than mu_t_outer and work out to find the first
            //          point for which that is not true.
            j_crossover = jfirst - 1; // set it safely, just in case...
            for ( j = jfirst - 1; j >= j_limit; --j ) {
                j_crossover = j;
                if ( mu_t_inner[j] > mu_t_outer[j] ) {
   		    // Add 1 to prevent an overshoot in the viscosity
		    if ( j_crossover < jfirst-1 ) j_crossover += 1;   
		    break;
	        }
            }

            // Step 13: Select the turbulent viscosity from the appropriate layer.
            for ( j = jfirst; j >= jlast; --j ) {
                if ( j >= j_crossover )
                    mu_t[j] = mu_t_inner[j];
                else
                    mu_t[j] = mu_t_outer[j];
            }

            // Step 14: Keep the turbulent transport coefficients, maybe.
            //
            // Transition extension to the model:
            // if the maximum mu_t in the profile from the wall for this x index is less 
            // than some specified transition level, discard the turbulent transport.
#           if TRANSITIONAL == 1
  	    double mu_t_max = 0.0;
	    for ( j = jfirst; j >= jlast; --j ) {
	        if ( mu_t[j] > mu_t_max ) mu_t_max = mu_t[j];
	    }
	    // comparing to the viscosity of the cell in the middle of the duct
	    // which is assumed to be at index jmin for "tube-type" calculations.
	    double mu_freestream = blk.get_cell(i,blk.jmin)->fs->gas->mu;
	    keep_turbulent_viscosity = (mu_t_max >= C_mutm * mu_freestream);
#           else
            keep_turbulent_viscosity = 1;
#           endif

            if ( keep_turbulent_viscosity ) {
	        for ( j = jfirst; j >= jlast; --j ) {
	  	    cell = blk.get_cell(i,j);
		    cell->fs->mu_t = mu_t[j];
		    cell->fs->k_t = gm->Cp(*(cell->fs->gas), status_flag) * cell->fs->mu_t / Pr_t;
	        }
	    } else {
	        for ( j = jfirst; j >= jlast; --j ) {
		    cell = blk.get_cell(i,j);
		    cell->fs->mu_t = 0.0;
		    cell->fs->k_t = 0.0;
	        }
 	    }
        } // end of i-loop
    } // end of "If NORTH boundary is the wall..."


    // Else if SOUTH boundary is the wall ..
    else if ( blk.bcp[SOUTH]->is_wall() && blk.bcp[SOUTH]->type_code != SLIP_WALL ) {
    
        // Index limits.
        ifirst = blk.imin;
        ilast = blk.imax;
        jfirst = blk.jmin; // Start at SOUTH wall ...
        jwall = blk.jmin;
        jlast = blk.jmax; // ...and move to the centre of the tube.
        // Set the search limit for determining inner/outer layer.
        j_limit = blk.jmin + (blk.nnj / 4);

        for ( i = ifirst; i <= ilast; ++i ) {
            /*
             * Before proceeding with the turbulence calculation,
             * check to see if the gas is moving with reasonable speed.
             * Use U_MAX later in step 9. Note that U_MAX is presently the
	     * velocity squared.
             */
            U_MAX = 0.0;
            for ( j = jfirst; j <= jlast; ++j ) {
                cell = blk.get_cell(i,j);
                temporary = cell->fs->vel.x * cell->fs->vel.x + cell->fs->vel.y * cell->fs->vel.y;
                U_MAX = max(U_MAX, temporary);
            }
            if (U_MAX < 10.0) {
	        for ( j = jfirst; j <= jlast; ++j ) {
		    cell = blk.get_cell(i,j);
		    cell->fs->mu_t = 0.0;
		    cell->fs->k_t = 0.0;
	        }
	        continue; // skip to the next i-position.
	    }

            // Step 1: Compute the vorticity magnitude at the vertical interface
            //         midpoints. Compute the normal coordinate assuming that 
            //         the j index direction is aligned with the wall-normal direction.
	    //
	    // Vertex indices as shown here.
	    // 3------2
   	    // |      |
	    // | cell |
	    // |      |
	    // 0------1
	    cell = blk.get_cell(i,jwall);
            y_wall = 0.5 * (cell->vtx[0]->pos[gtl].y + cell->vtx[1]->pos[gtl].y);
            x_wall = 0.5 * (cell->vtx[0]->pos[gtl].x + cell->vtx[1]->pos[gtl].x);
            for ( j = jwall; j <= jlast; ++j ) {
	        cell = blk.get_cell(i,j);
                temporary = 0.25*(cell->vtx[0]->dudy + cell->vtx[1]->dudy + 
	    		          cell->vtx[2]->dudy + cell->vtx[3]->dudy)
		    - 0.25*(cell->vtx[0]->dvdx + cell->vtx[1]->dvdx + 
		  	    cell->vtx[2]->dvdx + cell->vtx[3]->dvdx);
                vort[j] = fabs(temporary);
                dy = cell->pos[gtl].y - y_wall;
                dx = cell->pos[gtl].x - x_wall;
                eta[j] = sqrt(dx * dx + dy * dy);
            }

            // Step 2: Compute the maximum shear stress --
            //         we will assume that it occurs at the wall.
	    cell = blk.get_cell(i,jwall);
            mu_wall = cell->fs->gas->mu;
            tau_max = mu_wall * vort[jwall];
            rho_wall = cell->fs->gas->rho;

            // Step 3: Compute y-plus across the domain.
            temporary = sqrt(rho_wall * tau_max) / mu_wall;
            for ( j = jfirst; j <= jlast; ++j )
                y_plus[j] = temporary * eta[j];

            // Step 4: Compute the vanDriest damping factor
            for ( j = jfirst; j <= jlast; ++j ) {
#               if COMPRESSIBILITY == 1
	        cell = blk.get_cell(i,j);
                factor = sqrt(cell->fs->gas->rho / rho_wall) * mu_wall / cell->fs->gas->mu;
#               else
                factor = 1.0;
#               endif
                D[j] = 1.0 - exp(-y_plus[j] * factor / A_plus);
            }

            // Step 5: Compute the inner-layer turbulent viscosity.
            for ( j = jfirst; j <= jlast; ++j ) {
                temporary = k_inner * eta[j] * D[j];
                mu_t_inner[j] = blk.get_cell(i,j)->fs->gas->rho * temporary * temporary * vort[j];
            }

            // Step 6: Start computing the pieces for the outer-layer viscosity.
	    //         Compute F.
	    for ( j = jfirst; j <= jlast; ++j ) {
                F[j] = eta[j] * vort[j] * D[j];
            }

            // Step 7: Find the maximum of F and the point at which it occurs.
            F_MAX = F[jfirst];
            eta_MAX = eta[jfirst];
            for ( j = jfirst + 1; j <= jlast; ++j ) {
                if ( F[j] > F_MAX ) {
                    F_MAX = F[j];
                    eta_MAX = eta[j];
	        }
            } // end for j
            // If F_MAX is zero, there isn't anything happening at this i station; 
            // skip to the next i station.
            if ( F_MAX < 1.0e-10 ) {
                for ( j = jfirst; j <= jlast; ++j ) {
		    cell = blk.get_cell(i,j);
		    cell->fs->mu_t = 0.0;
		    cell->fs->k_t = 0.0;
	        }
                continue; // ...with next i-station
 	    }

            // Step 8: Compute the Klebanoff Intermittency factor.
            for ( j = jfirst; j <= jlast; ++j ) {
                temporary = C_KLEB * eta[j] / eta_MAX;
                F_KLEB[j] = 1.0 / (1.0 + 5.5 * pow(temporary,6) );
            }

            // Step 9: Find the maximum and minimum total velocities.
            //         Currently set U_MIN = 0.
            U_MIN = 0.0;
            U_MAX = sqrt(U_MAX); // square of total velocity was computed before
            U_DIFF = U_MAX - U_MIN;

            // Step 10: Compute F_WAKE.
            F_A = eta_MAX * F_MAX;
            if ( F_MAX > 1.0e-10 )
                F_B = C_WK * eta_MAX * U_DIFF * U_DIFF / F_MAX;
            else
                F_B = F_A;
            F_WAKE = min(F_B, F_A);

            // Step 11: Compute the outer-layer turbulent viscosity.
            for ( j = jfirst; j <= jlast; ++j ) {
                mu_t_outer[j] = K * C_CP * blk.get_cell(i,j)->fs->gas->rho * F_WAKE * F_KLEB[j];
            }

            // Step 12: Find the cross-over point for the inner- and
            //          outer-layers.  Assume that, near the wall, mu_t_inner
            //          is smaller than mu_t_outer and work out to find the first
            //          point for which that is not true.
            j_crossover = jfirst + 1; // set it safely, just in case...
            for ( j = jfirst + 1; j <= j_limit; ++j ) {
                j_crossover = j;
                if ( mu_t_inner[j] > mu_t_outer[j] ) {
                    // Substract 1 to prevent an overshoot in the viscosity
                    if ( j_crossover > jfirst+1 ) j_crossover -= 1;
		    break;
	        }
            }

            // Step 13: Select the turbulent viscosity from the appropriate layer.
            for ( j = jfirst; j <= jlast; ++j ) {
                if ( j <= j_crossover )
                    mu_t[j] = mu_t_inner[j];
                else
                    mu_t[j] = mu_t_outer[j];
            }

            // Step 14: Keep the turbulent transport coefficients, maybe.
            //
            // Transition extension to the model:
            // if the maximum mu_t in the profile from the wall for this x index is less 
            // than some specified transition level, discard the turbulent transport.
#           if TRANSITIONAL == 1
	    mu_t_max = 0.0;
            for ( j = jfirst; j <= jlast; ++j ) {
	        if ( mu_t[j] > mu_t_max ) mu_t_max = mu_t[j];
	    }
	    // comparing to the viscosity of the cell in the middle of the duct
	    // which is assumed to be at index jmax for "tube-type" calculations.
	    mu_freestream = blk.get_cell(i,blk.jmax)->fs->gas->mu;
	    keep_turbulent_viscosity = (mu_t_max >= C_mutm * mu_freestream);
#           else
            keep_turbulent_viscosity = 1;
#           endif

            if ( keep_turbulent_viscosity ) {
                for ( j = jfirst; j <= jlast; ++j ) {
		    cell = blk.get_cell(i,j);
		    cell->fs->mu_t = mu_t[j];
		    cell->fs->k_t = gm->Cp(*(cell->fs->gas), status_flag) * cell->fs->mu_t / Pr_t;
	        }
	    } else {
                for ( j = jfirst; j <= jlast; ++j ) {
		    cell = blk.get_cell(i,j);
		    cell->fs->mu_t = 0.0;
		    cell->fs->k_t = 0.0;
	        }
	    }
        } // end of i-loop
    } // end of "If SOUTH boundary is the wall..."


    // Else ..
    else {
    //    cout << "baldwin_lomax_turbulence_model(): "
    //         << "no acceptable face being chosen as solid wall\n"
    //         << "Currently, only NORTH or SOUTH face can be accepted "
    //         << "as solid walls." << endl;
    } // end of "else"

    
    vort.clear();
    eta.clear();
    y_plus.clear();
    D.clear();
    F.clear();
    F_KLEB.clear();
    mu_t_inner.clear();
    mu_t_outer.clear();
    mu_t.clear();
    return SUCCESS;
} // end of baldwin_lomax_turbulence_model()

