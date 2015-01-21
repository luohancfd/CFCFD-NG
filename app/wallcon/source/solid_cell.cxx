#include <iostream>

#include "solid_cell.hh"
#include <iomanip>


Solid_FV_Interface::Solid_FV_Interface() //default constructor
    :pos(2,0.0), length(0.0), normal(2,0.0), tangent(2,0.0), q(2,0.0), T(0.0) //initialsise. vector(2 elements, initial value)
{}

void Solid_FV_Interface::print_pos(){
    std::cout << "(" << pos[0] << "," << pos[1] << ")" << std::endl;
}

void Solid_FV_Interface::print_normal(){
    std::cout << "(" << normal[0] << "," << normal[1] << ")" << std::endl;
}

void Solid_FV_Interface::print_tangent(){
    std::cout << "(" << tangent[0] << "," << tangent[1] << ")" << std::endl;
}

void Solid_FV_Interface::print_T(){
    std::cout << T << std::endl;
}

void Solid_FV_Interface::print_length(){
    std::cout << length << std::endl;
}

void Solid_FV_Interface::print_q(){
    std::cout << "(" << q[0] << "," << q[1] << ")" << std::endl;
}

Solid_FV_Vertex::Solid_FV_Vertex()//default constructor
    :pos(2,0.0),volume(0.0), grad_T(2,0.0)
{
    iface.push_back(new Solid_FV_Interface() );
    iface.push_back(new Solid_FV_Interface() );
    iface.push_back(new Solid_FV_Interface() );
    iface.push_back(new Solid_FV_Interface() );
}

void Solid_FV_Vertex::print_pos(){
    std::cout << "(" << pos[0] << "," << pos[1] << ")" << std::endl;
}

void Solid_FV_Vertex::print_gradT(){
    std::cout << "(" << grad_T[0] << "," << grad_T[1] << ")" << std::endl;
}


Solid_FV_Cell::Solid_FV_Cell()
    :pos(2,0.0), volume(0.0), rho(0.0), k(0.0), cp(0.0), k11(0.0), k12(0.0), k22(0.0),
     e(0.0), e0(0.0), T(0.0), T1(0.0), deondt(0), de0ondt(0.0), source(0)
{
    //This doesn't like to init above. Put here instead.
    vrtx.push_back(new Solid_FV_Vertex() );
    vrtx.push_back(new Solid_FV_Vertex() );
    vrtx.push_back(new Solid_FV_Vertex() );
    vrtx.push_back(new Solid_FV_Vertex() );
    iface.push_back(new Solid_FV_Interface() );
    iface.push_back(new Solid_FV_Interface() );
    iface.push_back(new Solid_FV_Interface() );
    iface.push_back(new Solid_FV_Interface() );
}


void Solid_FV_Cell::print_pos(){
    std::cout << "(" << pos[0] << "," << pos[1] << ")" << std::endl;
}

void Solid_FV_Cell::print_source(){
    std::cout << "Source: " << source << std::endl;
}

void Solid_FV_Cell::print_T(){
    std::cout << "Temperature: " << T << std::endl;
}

void Solid_FV_Cell::print_V(){
    std::cout << "Volume: " << volume << std::endl;
}

void Solid_FV_Cell::print_deondt(){
    std::cout << "Energy derivative " <<  deondt << std::endl;
}

void Solid_FV_Cell::one_stage_update(double &dt){
    e = e + dt * deondt;
    T = e/ ( rho * cp);
}

void Solid_FV_Cell::two_stage_update(double &dt){
    e = e0 + 0.5*dt*( de0ondt + deondt);
    //Checking differences
    // std::cout << "de0ondt "<< de0ondt << " deondt " << deondt << std::endl;
    // std::cout << std::setprecision(12) << "T0 " << T;
    T = e/ ( rho * cp);
    // std::cout << " T " << T << std::endl;
}
