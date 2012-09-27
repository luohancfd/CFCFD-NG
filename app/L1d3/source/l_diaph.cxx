// l_diaph.cxx
// Refactored from l1d code 26-Sep-2012

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "l_diaph.hh"
#include "l1d.hh"

DiaphragmData::DiaphragmData(int indx, std::string config_file_name, int echo_input)
{
    ConfigParser dict = ConfigParser(config_file_name);
    std::stringstream tag;
    tag << indx;
    std::string section = "diaphragm-" + tag.str();
    if (echo_input == 1) cout << "Reading diaphragm " 
			      << indx << " parameters..." << endl;

    dict.parse_int(section, "is_burst", is_burst, 0);
    dict.parse_double(section, "p_burst", P_burst, 0.0);
    dict.parse_double(section, "dt_hold", hold_period, 0.0);
    dict.parse_double(section, "dt_blend", blend_delay, 0.0);
    dict.parse_double(section, "dx_blend", blend_dx, 0.0);
    if (echo_input == 1) {
	cout << "    is_burst = " << is_burst << endl;
	cout << "    p_burst = " << P_burst << endl;
	cout << "    dt_hold = " << hold_period << endl;
	cout << "    dt_blend = " << blend_delay << endl;
	cout << "    dx_blend = " << blend_dx << endl;
    }
    // Initially set the trigger_time to a negative number.
    trigger_time = -1.0;
    already_blended = 0;
    // By default, the diaphragm is unconnected.
    left_slug_id = -1;
    left_slug_end_id = -1;
    dict.parse_int(section, "left-slug-id", left_slug_id, -1);
    std::string label;
    dict.parse_string(section, "left-slug-end-id", label, "");
    if (label[0] == 'L' || label[0] == 'l' || label[0] == '0')
        left_slug_end_id = LEFT;
    if (label[0] == 'R' || label[0] == 'r' || label[0] == '1')
        left_slug_end_id = RIGHT;
    dict.parse_double(section, "dxL", left_slug_dx, 0.0);
    if (echo_input == 1) {
	cout << "    left-slug-id = " << left_slug_id << endl;
	cout << "    left-slug-end-id = " << left_slug_end_id << endl;
	cout << "    dxL = " << left_slug_dx << endl;
    }
    right_slug_id = -1;
    right_slug_end_id = -1;
    dict.parse_int(section, "right-slug-id", right_slug_id, -1);
    dict.parse_string(section, "right-slug-end-id", label, "");
    if (label[0] == 'L' || label[0] == 'l' || label[0] == '0')
        right_slug_end_id = LEFT;
    if (label[0] == 'R' || label[0] == 'r' || label[0] == '1')
        right_slug_end_id = RIGHT;
    dict.parse_double(section, "dxR", right_slug_dx, 0.0);
    if (echo_input == 1) {
	cout << "    right-slug-id = " << right_slug_id << endl;
	cout << "    right-slug-end-id = " << right_slug_end_id << endl;
	cout << "    dxR = " << right_slug_dx << endl;
    }
}


DiaphragmData::DiaphragmData(const DiaphragmData& dd)
{
    sim_time = dd.sim_time;
    is_burst = dd.is_burst;
    P_burst = dd.P_burst;
    hold_period = dd.hold_period;
    trigger_time = dd.trigger_time;
    already_blended = dd.already_blended;
    blend_dx = dd.blend_dx;
    blend_delay = dd.blend_delay;
    left_slug_id = dd.left_slug_id;
    right_slug_id = dd.right_slug_id;
    left_slug_end_id = dd.left_slug_end_id;
    right_slug_end_id = dd.right_slug_end_id;
    left_slug_dx = dd.left_slug_dx;
    right_slug_dx = dd.right_slug_dx;
} // end DiaphragmData copy constructor


DiaphragmData::~DiaphragmData()
{
    // nothing to do
}


int DiaphragmData::read_state(FILE* infile)
// Read the diaphragm state from an already open file.
{
#   define NCHAR 320
    char line[NCHAR];
    int nread;
    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty solution file.\n");
        return FAILURE;
    }
    nread = sscanf(line, "%lf", &sim_time);
    if ( nread != 1 ) {
	printf( "read_diaphragm_solution(): didn't correctly read sim_time\n" );
        printf( "from line:\n%s\n", line );
	return FAILURE ;
    }
    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty solution file.\n");
        return FAILURE;
    }
    nread = sscanf(line, "%d %d %lf", &is_burst, &already_blended, &trigger_time);
    if ( nread != 3 ) {
	printf( "read_diaphragm_solution(): " );
	printf( "didn't correctly read is_burst, already_blended, trigger_time\n" );
        printf( "from line:\n%s\n", line );
	return FAILURE;
    }
    return SUCCESS;
#   undef NCHAR
}


int DiaphragmData::write_state(FILE* outfile)
// Write the diaphragm state to an already open file.
{
    fprintf(outfile, "%e  # begin diaphragm data: sim_time\n", sim_time);
    fprintf(outfile, "%d %d %e  # is_burst, already_blended, trigger_time\n", 
	    is_burst, already_blended, trigger_time);
    fflush(outfile);
    return SUCCESS;
}

