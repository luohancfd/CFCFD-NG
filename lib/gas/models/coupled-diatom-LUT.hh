// Author: Daniel F. Potter
// Version: 09-Apr-2012

#ifndef COUPLED_DIATOM_LUT_HH
#define COUPLED_DIATOM_LUT_HH

#include <vector>
#include <string>

class NoneqCoupledDiatomicLUT {
public:
    NoneqCoupledDiatomicLUT( std::string fname, int icol );
    ~NoneqCoupledDiatomicLUT();
    
    std::string str();

    double eval( double T_el, double T_vib, double T_rot );

private:
    void find_bounding_temperature_indices( double T, int &i_l, int &i_u );
    
private:
    std::string fname_;
    int icol_;
    std::vector<double> T_list_;
    double *** data_;
    int ndim_;
};

class EqCoupledDiatomicLUT {
public:
    EqCoupledDiatomicLUT( std::string fname, int icol );
    ~EqCoupledDiatomicLUT();
    
    std::string str();

    double eval( double T );

private:
    void find_bounding_temperature_indices( double T, int &i_l, int &i_u );
    
private:
    std::string fname_;
    int icol_;
    std::vector<double> T_list_;
    std::vector<double> data_;
    int ndim_;
};

#endif
