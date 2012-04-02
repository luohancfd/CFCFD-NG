/** \file gas_data.hh
 *  \ingroup gas
 *
 *  \author Rowan J Gollan
 *  \version 24-May-2008 -- initial coding, refactored
 *                          from gas_models2/source/gas.hh
 **/

#ifndef GAS_DATA_HH
#define GAS_DATA_HH

#include <vector>
#include <string>

typedef std::vector<std::vector<double> > matrix;

double get_matrix_element(const matrix &M, int i, int j );
int copy_matrix_elements(const matrix &srcM, matrix &destM );
int average_matrix_elements(const matrix &srcM0, const matrix &srcM1, matrix &destM);

class Gas_model;

class Gas_data {
public:
    // -- Thermodynamic properties
    double rho;                /**< \brief density, kg/m**3               */
    double p;                  /**< \brief [total] pressure, Pa           */
    double p_e;                /**< \brief electron pressure, Pa          */
    double a;                  /**< \brief speed of sound, m/s            */
    std::vector<double> e;     /**< \brief specific internal energy, J/kg */
    std::vector<double> T;     /**< \brief array of temperatures, K       */
    // -- Transport properties
    double mu;                 /**< \brief viscosity, Pa.s                 */
    std::vector<double> k;     /**< \brief array of heat flux coefficients */
    matrix D_AB;               /**< \brief binary diffusion coefficients.  */
    // -- composition
    std::vector<double> massf; /**< \brief species mass fractions          */
    //
    Gas_data(Gas_model *gm);
    Gas_data(const Gas_data &Q);
    Gas_data(int nsp, int nmodes);
    ~Gas_data();
    Gas_data & operator= (const Gas_data &Q);
    void print_values(bool print_transport_data=true) const;
    void write_values(std::string fname, bool print_transport_data) const;
    void copy_values_from(const Gas_data &src);
    void average_values_from(const Gas_data &src0, const Gas_data &src1, bool with_diff_coeff);
    int accumulate_values_from(const Gas_data &src, double alpha);
    double* copy_values_to_buffer(double *buf) const;
    double* copy_values_from_buffer(double *buf);
    int check_values(bool print_message=true ) const;
};

double mass_average(const std::vector<double> &massf, const std::vector<double> &vec);
double mass_average_inverse(const std::vector<double> &massf, const std::vector<double> &vec);
double mole_average(const std::vector<double> &molef, const std::vector<double> &vec);
int scale_mass_fractions(std::vector<double> &massf, double tolerance=0.0);
int copy_mass_fractions(const std::vector<double> &f_src, std::vector<double> &f_target);

#endif
