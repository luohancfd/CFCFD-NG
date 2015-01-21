#ifndef SOLID_CELL_HH
#define SOLID_CELL_HH

#include <vector>


//fluxes between cell interfaces at boundaries
class Solid_FV_Interface {
public: //public access permission

    //Geometry
    std::vector<double> pos; //mid position of interface based on vertex pos
    double length; //side length
    std::vector<double> normal; //normal vector of interface
    std::vector<double> tangent;

    std::vector<double> q; // heat flux
    double T; //Temperature

    //Constructors
    Solid_FV_Interface();

    //member functions
    void print_pos();
    void print_normal();
    void print_tangent();
    void print_T();
    void print_q();
    void print_length();


};

//holds data for calculation of gradients
class Solid_FV_Vertex {
public:

    //Geomerty
    std::vector<double> pos; //vertex position (grid node position)


    //member fuctions
    void print_pos();
    void print_gradT();

    //Secondary cell information
    double volume; // secondary cell volume for grad_q calc
    std::vector<double> grad_T; //flux derivatives
    std::vector<Solid_FV_Interface *> iface; // pointer to secondary cell normals

    //Consturctor
    Solid_FV_Vertex();
};

//cell-averaged data
class Solid_FV_Cell {
public:

    //Geometry
    std::vector<double> pos;//cell centre position approx
    double volume; //volume of cell

    //Properties
    //Maybe assign only homogenious properties to block. Anisotropy is different but can also be block based. create connecting blocks to allow non-homgenity
    double rho,k,cp;
    double k11,k12,k22;

    //Energy
    double e; //energy perunit volume
    double e0;
    double T; //temperature
    double T1; //temperature
    double deondt; // time derivatives for update stage
    double de0ondt;
    double source; //Source term

    //Boundaryflags
    int is_on_boundary;
    int which_boundary;

    //connections
    std::vector<Solid_FV_Vertex *> vrtx; //Pointers to defining vertex
    std::vector<Solid_FV_Interface *> iface; //pointers to defining interfaces

    Solid_FV_Cell();
    void print_pos();
    void print_T();
    void print_source();
    void print_deondt();
    void print_V();
    void one_stage_update(double &dt);
    void two_stage_update(double &dt);
    
};

#endif // SOLID_CELL_HH
