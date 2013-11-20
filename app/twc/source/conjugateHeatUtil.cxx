#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "../../../lib/util/source/useful.h"
#include "init.hh"
#include "conjugateHeatUtil.hh"

void wall_printer(valarray<double> &wall, list_of_inputs &inputs) {

/*	---------------------------------------------------------------------
	---------------------------------------------------------------------

	This is a debugging function which prints out any Valarray describing
	nodes on the wall in the shape of the wall.  It gives a more intuitive
	view of the wall, making it easier to check that correct parameters
	are assigned to nodes.

	---------------------------------------------------------------------

	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------
	---------------------------------------------------------------------

*/

	int counter = 0;
	int width = 10;
	cout << endl;
	for (int i = 0; i<inputs.N; i++) {
		cout << "[ ";
		for (int j = 0; j < inputs.M; j++) {
			cout << setprecision(5) << setw(width) << wall[counter] << " ";
			counter++;
		}
		cout << "]" << endl;
	}
	cout << endl;

}

bool file_exists(const char *filename) {
/*

	---------------------------------------------------------------------
	---------------------------------------------------------------------
	This is just a utility function with checks if a file exists.

	Returns 1 if the file exists.  
	
	It's taken directly from a c++ forum:
	http://www.cplusplus.com/forum/general/1796/


	---------------------------------------------------------------------

	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------
	---------------------------------------------------------------------

*/	
	ifstream ifile(filename);
	return ifile;
}

int write_row(list_of_vars &vars, list_of_inputs &inputs, int row) {

	/*
	--------------------------------------------------------------------------
	--------------------------------------------------------------------------
	This function writes a particular row out to a text file.  The format of the 
	text file is just x and T separated by a space, allowing it to be plotted
	simply using GNU Plot.

	24/10/2013 - At the moment just use GNUPlot from the command window


	---------------------------------------------------------------------------
	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------------
	---------------------------------------------------------------------------
	*/

	const char* filename = "../user/row_temp.txt";

	ofstream output;
	output.open(filename);

    for(int i = inputs.M*(row-1); i < inputs.M*row; i++) {
        output << vars.x_vals[i] << " " << vars.temps[i] << endl;
	}

	output.close();

	return 0;

}

int write_column(list_of_vars &vars, list_of_inputs &inputs, int column) {

	/*
	--------------------------------------------------------------------------
	--------------------------------------------------------------------------
	This function writes a particular column out to a text file.  The format of the 
	text file is just x and T separated by a space, allowing it to be plotted
	simply using GNU Plot.

	24/10/2013 - At the moment just use GNUPlot from the command window


	---------------------------------------------------------------------------
	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------------
	---------------------------------------------------------------------------
	*/

	const char* filename = "../user/column_temp.txt";

	ofstream output;
	output.open(filename);

    for(int i = (column-1); i < inputs.M*inputs.N; i+=inputs.M) {
        output << vars.y_vals[i] << " " << vars.temps[i] << endl;
	}

	output.close();

	return 0;

}


int write_soln(list_of_vars &vars, list_of_inputs &inputs, double time_elapsed, int print_number) {

    /*
      --------------------------------------------------------------------------
      --------------------------------------------------------------------------
      Writes the solution out in the same way that the Eilmer3 output is produced.
      
      This also creates and runs a GNUPlot script which plots the temperature 
      profile and stores it in the same directory as the printed results text
      file.
      
      
      ---------------------------------------------------------------------------
      \author Jared Clifford
      \version November 2013
      
      ---------------------------------------------------------------------------
      ---------------------------------------------------------------------------
    */

    char buffer[100];
    int flag;
    sprintf(buffer, "mkdir -p user/temp/t%04d", print_number);
    flag = system(buffer); // Runs the command string "buffer" in terminal
    
    sprintf(buffer, "user/temp/t%04d/conjugateHeat.temp.b0000.t%04d.txt", print_number, print_number);
    ofstream output;
    output.open(buffer);
    
    output << time_elapsed << endl;
    output << "''pos.x'' ''pos.y'' ''temp''" << endl;
    output << inputs.M << " " << inputs.N << endl;
    for (int i = 0; i< inputs.M*inputs.N; i++) {
	output << vars.x_vals [i] << " " << vars.y_vals[i] << " " << vars.temps[i] << endl;
    }
    output.close();
    
    // Next part of the code writes the GNUPlot script for graphing each separate temp profile.
    
    ofstream op;
    sprintf(buffer, "user/temp_pos_plot.p");
    op.open(buffer);
    op << "set term png" << endl;
    sprintf(buffer, " set output \"user/temp/t%04d/pos_plot.png\"", print_number);
    op << buffer << endl;
    sprintf(buffer, "set title \"Temperature versus Position after %.3f seconds\"", time_elapsed);
    op << buffer << endl;
    op << "set xlabel \"Position (m)\"" << endl;
    op << "set ylabel \"Temperature (degrees C)\"" << endl;
    op << "set key bottom left" << endl;
    //op << "set xrange [0:0.025]" << endl;
    //op << "set yrange [0:0.025]" << endl;
    op << "set cbrange [295:355]" << endl;
    sprintf(buffer, "plot \"user/temp/t%04d/conjugateHeat.temp.b0000.t%04d.txt\" using 1:2:3 with image", print_number, print_number);
    op << buffer << endl;
    op.close();
    
    sprintf(buffer, "gnuplot \"user/temp_pos_plot.p\"");
    flag = system(buffer);
    if ( flag != 0 ) return FAILURE;
    return SUCCESS;
}

int write_hist(list_of_vars &vars, list_of_inputs &inputs, double time_elapsed) {

	/*
	--------------------------------------------------------------------------
	--------------------------------------------------------------------------
	Writes the temperature out at the selected node to the history file.



	---------------------------------------------------------------------------
	\author Jared Clifford
	\version November 2013

	---------------------------------------------------------------------------
	---------------------------------------------------------------------------
	*/

	char filename[50] = "../user/hist/hist.txt";

	ofstream output;
	output.open(filename, ofstream::app);
	output << time_elapsed << " " << vars.temps[inputs.nn_hist]-273. << endl;
	output.close();
	return SUCCESS;
}


int old_file_delete_catcher(list_of_inputs &inputs) {

	/*
	--------------------------------------------------------------------------
	--------------------------------------------------------------------------
	This function is designed to try and catch myself before overwriting
	old files.  It checks to see if the temp and hist directories exist 
	(indicating that an old solution is still in memory).  If the directory/ies
	exist, then it prompts the user as to whether they would like to archive the
	existing directory.  Archiving is performed via a shell script which stores 
	the temp and hist directories in a directory named after the current date
	and time.


	---------------------------------------------------------------------------
	\author Jared Clifford
	\version November 2013

	---------------------------------------------------------------------------
	---------------------------------------------------------------------------
	*/
	
	const char* temp_dir = "../user/temp"; const char* hist_dir = "../user/hist";
	bool temp_dir_exists = file_exists(temp_dir);
	bool hist_dir_exists = file_exists(hist_dir);
	int flag = 0;

	string yn;

	if (temp_dir_exists == 1 || hist_dir_exists == 1) {
	    cout << "Existing solution files are present.  Erase (selecting no means files are archived)? (y/n)" << endl;
	    cin >> yn;
	    if (yn == "y" || yn == "Y") {
		if (temp_dir_exists == 1) {
		    flag = system("rm -r ../user/temp");
		}
		if (hist_dir_exists == 1) {
		    flag = system("rm -r ../user/hist");
		}
	    }
	    else if (yn == "n" || yn == "N") {
		flag = system("sh arch.sh");
	    }
	    else {
		cout << "Inappropriate letter given" << endl;
		exit(MISMATCHED_DIMENSIONS);
	    }
	}
	
	flag = system("mkdir ../user/temp");
	if (inputs.flag_hist == 1) {
	    flag = system("mkdir ../user/hist");
	}
	if ( flag != 0 ) return FAILURE;
	return SUCCESS;
}


double linear_interp(valarray<double> in, valarray<double> out, double searchval) {

	if (in.size() != out.size()) {
		cerr << "Valarrays aren't the same length- no interpolation possible" << endl;
		exit(MISMATCHED_DIMENSIONS);
	}

	double ans = 0; double frac = 0;
	size_t in_maxindex = 0;
	size_t in_minindex = 0;
	
	for ( size_t i = 0; i < in.size(); i++ ) {

		if (in[i] == in.max()) {
			in_maxindex = i;
		}
		if (in[i] == in.min()) {
			in_minindex = i;
		}
	}
	for ( size_t i = 0; i < in.size(); i++ ) {
		if ((searchval > in[i] && searchval < in[i+1]) || (searchval < in[i] && searchval > in[i+1])) {
			frac = (searchval-in[i])/(in[i+1]-in[i]);
			ans = out[i] + frac*(out[i+1]-out[i]);
		}
		else if (searchval < in.min()) {
			ans = out[in_minindex];
		}
		else if (searchval > in.max()) {
			ans = out[in_maxindex];
		}

	}
	
	return ans;

}

