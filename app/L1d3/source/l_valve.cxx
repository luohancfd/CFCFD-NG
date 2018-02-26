// l_valve.cxx
// Refactored from l1d code 26-Sep-2012

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "l_valve.hh"
#include "l1d.hh"

ValveData::ValveData(int indx, std::string config_file_name, int echo_input)
{
    ConfigParser dict = ConfigParser(config_file_name);
    std::stringstream tag;
    tag << indx;
    std::string section = "valve-" + tag.str();
    if (echo_input == 2) cout << "Reading valve " 
			      << indx << " parameters..." << endl;

    dict.parse_int(section, "is_open", is_open, 0);
    //dict.parse_double(section, "p_open", P_open, 0.0);
    dict.parse_double(section, "open_period", open_period, 0.0);
    dict.parse_double(section, "open_time", open_time, 0.0);
    //dict.parse_double(section, "dx_blend", blend_dx, 0.0);
    if (echo_input == 2) {
	cout << "    is_open = " << is_open << endl;
	cout << "    open_period = " << open_period << endl;
	cout << "    open_time = " << open_time << endl;
    }
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
    if (echo_input == 2) {
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
    if (echo_input == 2) {
	cout << "    right-slug-id = " << right_slug_id << endl;
	cout << "    right-slug-end-id = " << right_slug_end_id << endl;
	cout << "    dxR = " << right_slug_dx << endl;
    }
}

ValveData::ValveData(const ValveData& vd)
{
    sim_time = vd.sim_time;
    is_open = vd.is_open;
    //P_open = vd.P_open;
    open_period = vd.open_period;
    open_time = vd.open_time;
    //already_blended = vd.already_blended;
    //blend_dx = vd.blend_dx;
    //blend_delay = vd.blend_delay;
    left_slug_id = vd.left_slug_id;
    right_slug_id = vd.right_slug_id;
    left_slug_end_id = vd.left_slug_end_id;
    right_slug_end_id = vd.right_slug_end_id;
    left_slug_dx = vd.left_slug_dx;
    right_slug_dx = vd.right_slug_dx;
} // end ValveData copy constructor


ValveData::~ValveData()
{
    // nothing to do
}


int ValveData::read_state(FILE* infile)
// Read the valve state from an already open file.
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
	printf( "read_valve_solution(): didn't correctly read sim_time\n" );
        printf( "from line:\n%s\n", line );
	return FAILURE ;
    }
    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty solution file.\n");
        return FAILURE;
    }
    nread = sscanf(line, "%d %lf %lf", &is_open, &open_period, &open_time);
    if ( nread != 3 ) {
	printf( "read_valve_solution(): " );
	printf( "didn't correctly read is_open, open_period, open_time\n" );
        printf( "from line:\n%s\n", line );
	return FAILURE;
    }
    return SUCCESS;
#   undef NCHAR
}


int ValveData::write_state(FILE* outfile)
// Write the diaphragm state to an already open file.
{
    fprintf(outfile, "%e  # begin valve data: sim_time\n", sim_time);
    fprintf(outfile, "%d %e %e  # is_open, open_period, open_time\n", 
	    is_open, open_period, open_time);
    fflush(outfile);
    return SUCCESS;
}

