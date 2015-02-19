/**
 * This module provides routines to update a gas phase
 * chemistry system.
 *
 * Author: Rowan G.
 */

module kinetics.chemistry_update;

import std.algorithm;

immutable int MAX_SUBCYCLES = 10000; // maximum number of subcycles to perform over tInterval
immutable int MAX_STEP_ATTEMPTS = 3; // maximum number of attempts on any one step
immutable double DT_INCREASE_PERCENT = 10.0; // allowable percentage increase on succesful step
immutable double DT_DECREASE_PERCENT = 50.0; // allowable percentage decrease on succesful step
immutable double DT_REDUCTION_FACTOR = 2.0; // factor by which to reduce timestep
                                            // after a failed attempt
alias ChemistryStep = int(double[], double, ref double);
static bool working_memory_allocated = false;
static GasState Qinit;
static double[] conc;
static double[] concSave;
enum ResultOfStep { success, failure };

int update_chemistry(GasState Q, double tInterval, ref double dtSuggest,
		     in GasModel gmodel, in ReactionMechanism rmech, ChemistryStep cstep,
		     bool tightTempCoupling, int maxSubcycles=MAX_SUBCYCLE, int maxAttempts=MAX_STEP_ATTEMPTS)
{
    // 0. On entry take a copy of the GasState in case we bugger it up.
    if ( !working_memory_allocated ) {
	Qinit = new GasState(gmodel.n_species, gmodel.n_modes);
	conc.length = gmodel.n_species;
	concSave.length = gmodel.n_species;
	working_memory_allocated = true;
    }
    Qinit.copy_values_from(Q);
    massf2conc(Q.massf, gmodel.mol_masses, conc);
    // 1. Sort out the time step for possible subcycling.
    double t = 0.0;
    double h;
    if ( dtSuggest > tInterval )
	h = dtSuggest;
    else if ( dtSuggest <= 0.0 )
	h = estimateStepSize(conc);
    else
	h = dtSuggest;
    // 2. Now do the interesting stuff, increment species change
    // 2a. Evaluate the rate constants (kf, kb) at least initially
    rmech.eval_rate_constants(Q);
    // 2b. Begin cycling
    int cycle = 0;
    int attempt = 0;
    for ( ; cycle < MAX_SUBCYCLES; ++cycle ) {
	/* Copy the last good timestep before moving on.
	 * If we're near the end of the interval, we want
	 * the penultimate step. In other words, we don't
	 * want to store the fractional step of (tInterval - t)
	 * that is taken as the last step.
	 */
	dtSuggest = h;
	h = min(h, tInterval - t);
	attempt = 0;
	for ( ; attempt < MAX_STEP_ATTEMPTS; ++attempt ) {
	    if ( cstep(rmech, Q, conc, h, hSuggest) == success ) {
		/* We succesfully took a step of size h
		 * so increment the total time.
		 */
		t += h;
		/* We can now make some decision about how to
		 * increase the timestep. We will take some
		 * guidance from the ODE step, but also check
		 * that it's sensible. For example, we won't
		 * increase by more than 10% (or INCREASE_PERCENT).
		 * We'll also balk if the ODE step wants to reduce
		 * the stepsize by 50% (or DECREASE_PERCENT).
		 * In that case, we'll set the new
		 * timestep to 50% of what it was and keep
		 * on going. Our reasoning being that we were
		 * successful on the previous step and things
		 * shouldn't have changed by that much.
		 */
		double hMax = h*(1.0 + DT_INCREASE_PERCENT/100.0);
		double hMin = h*(1.0 - DT_DECREASE_PERCENT/100.0);
		h = min(hSuggest, hMax);
		h = max(h, hMin);
		break;
	    }
	    else { // in the case of failure...
		/* We now need to make some decision about 
		 * what timestep to attempt next. We follow
		 * David Mott's suggestion in his thesis (on p. 51)
		 * and reduce the timestep by a factor of 2 or 3.
		 * (The actual value is set as DT_REDUCTION_FACTOR).
		 */
		h *= DT_REDUCTION_FACTOR; 
	    }
	} // end attempts at single step.
	if ( attempt == MAX_STEP_ATTEMPTS ) {
	    // We did poorly. Let the outside world know by returning -1.
	    return -1;
	}
	/* Otherwise, we've done well.
	 * If tight temperature coupling is requested, we can reevaluate
	 * the temperature at this point. With regards to tight coupling,
	 * we follow Oran and Boris's advice on pp. 140-141.
	 * To paraphrase: solving a separate differential equation for
	 * temperature is computationally expensive, however, it usually
	 * suffices to reevaluate the temperature assuming that total internal
	 * energy of the system has not changed but has been redistributed
	 * among chemical states. Since the chemistry has not altered much, 
	 * the iteration for temperature should converge in one or two
	 * iterations.
	 *
	 * My own additional argument for avoiding a temperature differential
	 * equation is that it does not play nicely with the special
	 * ODE methods for chemistry that exploit structure in the species
	 * differential equations.
	 */
	if ( tightTempCoupling ) {
	    gmodel.conc2massf(conc, Q);
	    gmodel.eval_thermo_from_rhoe(Q);
	    rmech.eval_rate_constants(Q);
	}
    }
    if ( cycle == MAX_SUBCYCLES ) {
	// We stopped because things are taking too long.
	return -1;
    }
}
		     
