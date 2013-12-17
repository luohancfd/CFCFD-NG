/** \file fit_spline.cxx
 *
 *  \author Rowan J. Gollan
 *  \version 15-Jan-2007
 *
 **/
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>

#include "../../nm/source/fobject.hh"
#include "../../nm/source/nelmin.hh"
#include "gpath.hh"
#include "best_fit_spline.hh"

using namespace std;

void printUsage()
{
    cout << "fit_spline.x\n";
    cout << "Usage:\n";
    cout << "> fit_spline.x input.data nsp python|text\n";
    cout << " where...\n";
    cout << " input.data       : rows of data, x-y pairs on each row\n";
    cout << " nsp              : integer number of points to fit\n";
    cout << " python or text   : if python, a python representation is given\n";
    cout << "                  : if text, a plain list of spline points is given.\n";
    exit(1);
}

int main(int argc, char *argv[])
{

    if(argc != 4)
	printUsage();

    string input_file(argv[1]);
    int nsp = atoi(argv[2]);
    string output_type(argv[3]);

    if( nsp < 2 ) {
	cout << "No. spline points should be greater than 2.\n";
	cout << "Bailing Out!\n";
	exit(1);
    }
    
    MultivariateFunction *best_spline = new Best_fit_spline(input_file, 1.0, 0.005, 0.0, 0.0);

    double f_min = 0;
    int n_fe = 0;
    int n_restart = 0;
    
    // Set up initial guess.
    double x, y;
    ifstream infile;
    infile.open(input_file.c_str());
    if( infile.fail() ) {
	cout << "Could not open " << input_file << "; Bailing out!\n";
	exit(1);
    }

    vector<Vector3> points;
    while( ! infile.eof() ) {
	infile >> x;
	infile >> y;
	if( infile.eof() )
	    break;
	if( x > 0.0) continue;
	points.push_back(Vector3(x, y));
    }

    infile.close();
    
    vector<int> dn;
    int n = points.size();
    int nx = 0;
    for(int i = 0; i < (nsp-1); ++i ) {
	nx = n / (nsp-i-1);
	dn.push_back(nx);
	n = n - nx;
    }

    vector<double> xin;
    xin.push_back(points[0].x);
    int index = 0;
    for( size_t i = 0; i < (dn.size()-1); ++i ) {
	index += dn[i];
	xin.push_back(points[index].x);
	xin.push_back(points[index].y);
    }
    xin.push_back(points[points.size()-1].y);

    minimize(best_spline, xin, &f_min, &n_fe, &n_restart);

    // But we'd like the control points equidistant.
    vector<Vector3> p;
    p.push_back(Vector3(xin[0], 0.0));
    for(size_t i = 1; i < xin.size()-1; i += 2) {
	p.push_back(Vector3(xin[i],xin[i+1]));
    }
    p.push_back(Vector3(0.0, xin[xin.size()-1]));

    Path* myPath = new Spline(p);

    double length = myPath->length();
    double dL = length / (nsp-1);

    vector<Vector3> new_p;
    new_p.push_back(p[0]);
    double dummy_t = 0.0;
    for(int i = 1; i < (nsp-1); ++i) {
	new_p.push_back(myPath->point_from_length(i*dL, dummy_t));
    }
    new_p.push_back(p[p.size()-1]);

    cout << setprecision(6) << showpoint;

    if(output_type == "python" ) {
	cout << "Spline([\n";
	for( size_t i = 0; i < new_p.size(); ++i ) {
	    cout << "    Vector(" << setw(12) << new_p[i].x << "," << setw(12) << new_p[i].y << "),\n";
	}
	cout << "])\n";
    }
    else {
	for( size_t i = 0; i < new_p.size(); ++i) {
	    cout << setw(12) << new_p[i].x << " " << setw(12) << new_p[i].y << endl;
	}
    }

    return 0;

}



