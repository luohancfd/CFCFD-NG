//This will read in a gzip compressed file with vetex locations produced from e3prep.py
/*format of read in is:
  41 21 1  # ni nj nk
  0.000000000000e+00   0.000000000000e+00   0.000000000000e+00
  1.000000000000e-03   0.000000000000e+00   0.000000000000e+00
  2.000000000000e-03   0.000000000000e+00   0.000000000000e+00
  ...
*/

#include <iostream>
#include <stdlib.h>
#include <fstream> //class to read and write to files
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <cstdlib>   


#include "output.hh"

using namespace std;
//

int create_directory(string path){
    mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);//path, mode_t read write permissions
    
    // check directory exists
    struct stat sb;
    if (!(stat(path.c_str(), &sb)) == 0 && !(S_ISDIR(sb.st_mode))){
	cout << "Failed to create directory temperature" << endl;
	exit(1);
    };
    return 0;
}

int write_grid_to_file(SolidBlock &blk)
{
    cout <<"Writing to  grid to file...";
    int nni = blk.nni;
    int nnj= blk.nnj; //get dimensions
    
    ofstream grid_file;
    
    string grid = "solid/grid/cell_grid.dat";
    
    
    grid_file.open(grid.c_str());
    if(!(grid_file.is_open())){
	cout << "File did not open" <<endl;
	exit(1);
    };
    
    
    for (int j=0; j<(nnj-1); j++){
        for (int i=0; i<(nni-1); i++){
	    
            grid_file << blk.block_cells[i][j].pos[0] << "\t" << blk.block_cells[i][j].pos[1] << "\n";
	    
        }
    }
    
    grid_file.close ();
    cout << "Finished writing to file" << endl;
    return 0;
}

int write_solid_times(double time, int print_number, int first_flag){
    if( first_flag == 0){
	string fname = "solid/solid.times";

	ofstream t_file;
	t_file.open(fname.c_str());
	if (!(t_file.is_open())) {
	    cout << ".time file did not open" <<endl;
	    exit(1);
	};
	t_file << "# tindx sim_time\n";
	t_file.close();
    };

    string fname = "solid/solid.times";

    ofstream t_file;
    t_file.open(fname.c_str(), ios::app);
    if (!(t_file.is_open())) {
	cout << ".times file did not open" <<endl;
	exit(1);
    };
    t_file << setfill('0') << setw(4) << print_number << " ";
    t_file << std::scientific << std::setprecision(16) << time << "\n"; //current time
    t_file.close();
    return 0;
}


int write_temperatures_to_file(SolidBlock &blk, double time, int print_number)
{
    cout <<"Writing solid solution\ttime " << time << "\t\ttindx " << print_number;

    //Time of write
    stringstream block_id;
    block_id << setfill('0') << setw(4) << blk.id; //pad id with zeros
    stringstream print_num;
    print_num << setfill('0') << setw(4) << print_number ;

    create_directory("solid/temperature/t" + print_num.str() );
    
    int nni = blk.nni;
    int nnj= blk.nnj; //get number of nodes
    int nnk(2); //2D simulation
    double pos_z = 0.0;
    
    ofstream grid_file;
    ofstream T_file;
        
    string temp = "solid/temperature/t" + print_num.str();
    
    string temperature = temp + "/b" +  block_id.str() + ".t"+ print_num.str();
    
    T_file.open(temperature.c_str());
    
    if (!(T_file.is_open())) {
	cout << "Temperature file did not open" <<endl;
	exit(1);
	};

    
    //Write stuff to file
    T_file << std::scientific << std::setprecision(16) << time << "\n"; //current time
    T_file << "\"pos.x\" \"pos.y\" \"pos.z\" \"volume\" \"rho\" \"cp\" \"k\" \"source\"  \"e\" \"T\"\n"; //Cell properties
    T_file << nni-1 << " " << nnj-1 << " " << nnk-1 << " \n";//grid size number of cells
    for (int j=0; j<(nnj-1); j++){
        for (int i=0; i<(nni-1); i++){
            //cout << i << "," << j<<endl;
	    //
            T_file  << blk.block_cells[i][j].pos[0] << " "
		    << blk.block_cells[i][j].pos[1] << " "
		    << pos_z << " " //2D simulation
		    << blk.block_cells[i][j].volume << " "
		    << blk.block_cells[i][j].rho << " "
		    << blk.block_cells[i][j].cp << " "
		    << blk.block_cells[i][j].k << " "
		    << blk.block_cells[i][j].source << " "
		    << blk.block_cells[i][j].e << " "
		    << blk.block_cells[i][j].T << "\n";
	    
        }
    }
    T_file.close ();
    if ( print_number == 0  ) {
	write_solid_times(time, print_number, print_number);}
    else {
	write_solid_times(time, print_number, 1);}
    cout << "." << endl;
    return 0;
}

