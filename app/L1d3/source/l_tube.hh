// l_tube.hh
// Tube Area Specification

#ifndef L_TUBE_HH
#define L_TUBE_HH

#include <vector>
#include <string>
#include <stdio.h>

class TubeModel {
public:
    // This first part of the data is used internally. 
    int n;                        // Number of x-stations
    double x1;                    // Starting point
    double dx;                    // x-spacing between stations
    double x2;                     // End point
    std::vector<double> diam;     // Effective diameter
    std::vector<double> area;     // Area at each station
    std::vector<double> T_Wall;   // Wall temperature
    std::vector<double> K_over_L; // Loss coefficient at each
                                  // station (div by length)

    // This second part of the data is for the new user input.
    // Each segment of tube i=1..nseg is between xb[i-1] and xb[i]
    // and has area varying from diameters Diamb[i-1] to Diamb[i].
    // If linear[i-1] == 1 then the variation in diameter is linear.
    // Otherwise a cubic polynomial is used to give ddiam/dx == 0
    // at the end points.
    int nseg;
    std::vector<int>linear;
    std::vector<double> xb;
    std::vector<double> Diamb;
    
   int nv;
   std::vector<int>n_points;
   std::vector<double> x_loc;
   std::vector<double> d_max;

    // Loss coefficients are selectively applied to patches of the tube.
    int nKL;
    std::vector<double> xbeginK;
    std::vector<double> xendK;
    std::vector<double> K;

    // Wall temperature patches.
    int nT;
    std::vector<double> xbeginT;
    std::vector<double> xendT;
    std::vector<double> Tlocal;
    double Tnominal;

    TubeModel(std::string config_file_name, int echo_input=0);
    ~TubeModel();
    int read_area(std::string file_name);
    int write_area(std::string file_name);
    int write_dump_file(std::string file_name);
}; // end of class tube_data

#endif
