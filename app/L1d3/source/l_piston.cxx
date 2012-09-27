// piston.cxx
// Refactored from l1d code 26-Sep-2012

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "l_piston.hh"
#include "l1d.hh"


PistonData::PistonData(int indx, double dt_init, std::string config_file_name, int echo_input)
{
    ConfigParser dict = ConfigParser(config_file_name);
    std::stringstream tag;
    tag << indx;
    std::string section = "piston-" + tag.str();
    if (echo_input == 1) cout << "Reading piston " << indx << " parameters..." << endl;

    dict.parse_double(section, "mass", mass, 1.0);
    dict.parse_double(section, "diameter", diam, 1.0);
    dict.parse_double(section, "length", length, 1.0);
    const double myPI = 4.0*atan(1.0);
    area = myPI * 0.25 * diam * diam;
    dict.parse_double(section, "front_seal_f", front_seal_f, 0.0);
    dict.parse_double(section, "front_seal_area", front_seal_area, 0.0);
    dict.parse_double(section, "back_seal_f", back_seal_f, 0.0);
    dict.parse_double(section, "back_seal_area", back_seal_area, 0.0);
    if (echo_input == 1) {
	cout << "    mass = " << mass << endl;
	cout << "    diameter = " << diam << endl;
	cout << "    length = " << length << endl;
	cout << "    area = " << area << endl;
	cout << "    front_seal_f = " << front_seal_f << endl;
	cout << "    front_seal_area = " << front_seal_area << endl;
	cout << "    back_seal_f = " << back_seal_f << endl;
	cout << "    back_seal_area = " << back_seal_area << endl;
    }
    dict.parse_double(section, "p_restrain", p_restrain, 1.0);
    dict.parse_int(section, "is_restrain", is_restrain, 0);
    dict.parse_double(section, "x_buffer", x_buffer, 0.0);
    dict.parse_int(section, "hit_buffer", on_buffer, 0);
    hit_buffer_count = 0;  /* initialize counter for buffer strikes */
    V_buffer = 0.0;
    dict.parse_int(section, "with_brakes", with_brakes, 0);
    dict.parse_int(section, "brakes_on", brakes_on, 0);
    if (echo_input == 1) {
	cout << "    p_restrain = " << p_restrain << endl;
	cout << "    is_restrain = " << is_restrain << endl;
	cout << "    x_buffer = " << x_buffer << endl;
	cout << "    on_buffer = " << on_buffer << endl;
	cout << "    hit_buffer_count = " << hit_buffer_count << endl;
	cout << "    with_brakes = " << with_brakes << endl;
	cout << "    brakes_on = " << brakes_on << endl;
    }
    /* By default, the piston is free. */
    left_slug_id = -1;
    left_slug_end_id = -1;
    dict.parse_int(section, "left-slug-id", left_slug_id, -1);
    std::string label;
    dict.parse_string(section, "left-slug-end-id", label, "");
    if (label[0] == 'L' || label[0] == 'l' || label[0] == '0')
        left_slug_end_id = LEFT;
    if (label[0] == 'R' || label[0] == 'r' || label[0] == '1')
        left_slug_end_id = RIGHT;
    if (echo_input == 1) {
	cout << "    left-slug-id = " << left_slug_id << endl;
	cout << "    left-slug-end-id = " << left_slug_end_id << endl;
    }
    right_slug_id = -1;
    right_slug_end_id = -1;
    dict.parse_int(section, "right-slug-id", right_slug_id, -1);
    dict.parse_string(section, "right-slug-end-id", label, "");
    if (label[0] == 'L' || label[0] == 'l' || label[0] == '0')
        right_slug_end_id = LEFT;
    if (label[0] == 'R' || label[0] == 'r' || label[0] == '1')
        right_slug_end_id = RIGHT;
    if (echo_input == 1) {
	cout << "    right-slug-id = " << right_slug_id << endl;
	cout << "    right-slug-end-id = " << right_slug_end_id << endl;
    }
    // Initial position and velocity.
    dt = dt_init;
    dict.parse_double(section, "x0", x0, 0.0);
    dict.parse_double(section, "v0", V0, 0.0);
    x = x0;
    V = V0;
    if (echo_input == 1) {
	cout << "    x0 = " << x0 << endl;
	cout << "    v0 = " << V0 << endl;
    }
    // Mass decay
    dict.parse_double(section, "f_decay", f_decay, 0.0);
    dict.parse_double(section, "mass_limit", mass_limit, 0.0);
    if (echo_input == 1) {
	cout << "    f_decay = " << f_decay << endl;
	cout << "    mass_limit = " << mass_limit << endl;
    }
    if (f_decay != 0.0) 
	apply_decay = 1;
    else 
	apply_decay = 0;
}


PistonData::PistonData(const PistonData& pd)
{
    mass = pd.mass;
    mass_old = pd.mass_old;
    diam = pd.diam;
    area = pd.area;
    length = pd.length;
    front_seal_f = pd.front_seal_f;
    front_seal_area = pd.front_seal_area;
    back_seal_f = pd.back_seal_f;
    back_seal_area = pd.back_seal_area;
    with_brakes = pd.with_brakes;
    p_restrain = pd.p_restrain;
    x_buffer = pd.x_buffer;
    V_buffer = pd.V_buffer;
    f_decay = pd.f_decay;
    mass_limit = pd.mass_limit;
    apply_decay = pd.apply_decay;
    brakes_on = pd.brakes_on;
    brakes_on_old = pd.brakes_on_old;
    is_restrain = pd.is_restrain;
    is_restrain_old = pd.is_restrain_old;
    hit_buffer_count = pd.hit_buffer_count;
    hit_buffer_count_old = pd.hit_buffer_count_old;
    on_buffer = pd.on_buffer;
    on_buffer_old = pd.on_buffer_old;
    x0 = pd.x0;
    V0 = pd.V0;
    sim_time = pd.sim_time;
    x = pd.x;
    V = pd.V;
    x_old = pd.x_old;
    V_old = pd.V_old;
    Pf = pd.Pf;
    Pb = pd.Pb;
    Friction = pd.Friction;
    dt = pd.dt;
    for ( int i = 0; i < NL; ++i ) {
	DxDt[i] = pd.DxDt[i];
	DVDt[i] = pd.DVDt[i];
    }
    left_slug_id = pd.left_slug_id;
    right_slug_id = pd.right_slug_id;
    left_slug_end_id = pd.left_slug_end_id;
    right_slug_end_id = pd.right_slug_end_id;
} // end PistonData copy constructor


PistonData::~PistonData()
{
    // nothing to do
}


int PistonData::read_state(FILE* infile)
// Read the piston state from an already open file.
{
#   define NCHAR 320
    char line[NCHAR];
    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty solution file.\n");
        return FAILURE;
    }
    sscanf(line, "%lf", &sim_time);
    // Position, Velocity and acceleration.
    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty solution file.\n");
        return FAILURE;
    }
    sscanf(line, "%lf %lf %lf %lf", &x, &V, &(DVDt[0]), &mass );
    // Flags
    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty solution file.\n");
        return FAILURE;
    }
    sscanf(line, "%d %d %d", &is_restrain, &on_buffer, &brakes_on);
#   undef NCHAR
    return SUCCESS;
}


int PistonData::write_state(FILE* outfile)
// Write the piston state to an already open file.
{
    fprintf(outfile, "%e  # begin piston data: sim_time\n", sim_time);
    fprintf(outfile, "%e %e %e %e  # x, V, accel, mass\n", x, V, DVDt[0], mass);
    fprintf(outfile, "%d %d %d  # is_restrain, on_buffer, brakes_on\n",
	    is_restrain, on_buffer, brakes_on);
    fflush(outfile);
    return SUCCESS;
}


int PistonData::record_state(void)
// Record the piston position and velocity before attempting a time-step.
{
    x_old = x;
    V_old = V;
    brakes_on_old = brakes_on;
    hit_buffer_count_old = hit_buffer_count;
    on_buffer_old = on_buffer;
    is_restrain_old = is_restrain;
    mass_old = mass;
    return SUCCESS;
}


int PistonData::restore_state(void)
// Restore the piston position and velocity to that 
// which existed before attempting the time-step.
{
    x = x_old;
    V = V_old;
    brakes_on = brakes_on_old;
    hit_buffer_count = hit_buffer_count_old;
    on_buffer = on_buffer_old;
    is_restrain = is_restrain_old;
    mass = mass_old;
    return SUCCESS;
}


int PistonData::time_derivatives(int time_level, double sim_time)
// Compute the time derivatives for the piston dynamics
//
// Input...
// time_level : specifies where the computed derivatives
//              are to be stored.
// sim_time   : current simulation time
{
    double pressure_force, friction_force;
    // The tolerance below which the velocity is assumed zero.
    const double V_TOL=1.0e-6;

    // If we have a restrained piston, check for release.
    if ( is_restrain == 1 ) {
        if ( Pb > p_restrain ) is_restrain = 0;
    }
    // If the piston is still restrained, we are going nowhere.
    if ( is_restrain == 1 ) {
        /* Set speed and acceleration to zero. */
        V = 0.0;
        DxDt[time_level] = 0.0;
        DVDt[time_level] = 0.0;
        return 0;
    }

    // Check for buffer strike only while we're moving forward.
    if ( x > x_buffer && V > V_TOL && on_buffer == 0 ) {
        hit_buffer_count += 1;
        on_buffer = 1;
        // Record striking speed.
        V_buffer = V;
        printf("Buffer strike: speed = %e\n", V_buffer);
    }
    // If we were on the buffer and the velocity is now negative,
    // we must have moved off the buffer.
    if ( on_buffer == 1 && V < 0.0 ) {
        on_buffer = 0;
    }
    // Consequence of being on the buffer is that we have zero speed.
    // There may still be pressure and friction forces to consider.
    if ( on_buffer == 1 ) {
        // Set speed to zero for this step.
        V = 0.0;
    }

    // Check for the application of brakes.
    if ( (with_brakes == 1) && (V < 0.0) ) brakes_on = 1;
    // If the brakes are on, we are going nowhere.
    if ( brakes_on == 1 ) {
        // Set speed and acceleration to zero.
        V = 0.0;
        DxDt[time_level] = 0.0;
        DVDt[time_level] = 0.0;
        return SUCCESS;
    }

    // ***********************************************
    // Now continue on with a free-piston calculation.
    // ***********************************************
    // The (signed) pressure force.
    pressure_force = area * (Pb - Pf);
    // The magnitude of the friction force.
    friction_force = front_seal_f*front_seal_area*Pf + back_seal_f*back_seal_area*Pb;

    // Update the state vector for the piston dynamics.
    //
    // This is the velocity.
    DxDt[time_level] = V;
    // Now, do the acceleration.
    if ( V > V_TOL ) {
	// Moving forward, apply full friction in reverse 
	DVDt[time_level] = (pressure_force - friction_force) / mass;
    } else if ( V < -(V_TOL) ) {
	// Moving backward, apply full friction forward 
	DVDt[time_level] = (pressure_force + friction_force) / mass;
    } else {
	// We are effectively stationary.
	//
	// If the pressure force is larger than the friction
	// force then apply the difference else remain stationary.
	if ( fabs(pressure_force) > friction_force ) {
	    // Pressure force dominates, apply the remainder 
	    if ( pressure_force > 0.0 ) {
		// Accelerate forward.
		DVDt[time_level] = (pressure_force - friction_force) / mass;
	    } else {
		// Accelerate backwards
		DVDt[time_level] = (pressure_force + friction_force) / mass;
	    }
	} else {
	    // Friction force dominates; let's remain stationary. 
	    DVDt[time_level] = 0.0;
	} // end if sufficient pressure to accelerate
    } // end if ... moving or stationary

    return SUCCESS;
} // end time_derivatives


int PistonData::predictor_step(void)
// Use the time derivatives to advance the piston dynamics forward by time step dt.
{
    if ( is_restrain || brakes_on ) {
        x = x_old;
        V = 0.0;
    } else if ( on_buffer == 1 ) {
        x = x_old;
        V = 0.0 + dt * DVDt[0];
    } else {
        x = x_old + dt * DxDt[0];
        V = V_old + dt * DVDt[0];
    }
    // Apply mass decay for flagged and unrestrained pistons only.
    if ( is_restrain == 0 && apply_decay == 1 ) {
        // Limit the mass decay to a user-defined level
        if ( mass > mass_limit ) {
            mass = mass_old * ( 1.0 - f_decay * dt );
        } else {
            mass = mass_limit;
        }
    }
    return SUCCESS;
} // end predictor_step


int PistonData::corrector_step(void)
// Use the time derivatives to advance the piston dynamics forward by time step dt.
{
    if ( is_restrain || brakes_on ) {
        x = x_old;
        V = 0.0;
    } else if ( on_buffer == 1 ) {
        x = x_old;
        V = 0.0 + dt * 0.5 * (DVDt[1] + DVDt[0]);
    } else {
        x = x_old + dt * 0.5 * (DxDt[1] + DxDt[0]);
        V = V_old + dt * 0.5 * (DVDt[1] + DVDt[0]);
    }
    return SUCCESS;
} // end corrector_step

