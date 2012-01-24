#ifndef RR_COEFFS_HH
#define RR_COEFFS_HH

#include <string>
#include "../../util/source/config_parser.hh"
#include "../models/gas.hh"
#include "../models/molecule.hh"
#include "../models/particle_data.hh"

class RateCoeffModel {
public:
    RateCoeffModel();
    
    /// \brief Normal constructor
    RateCoeffModel( const std::string model, const std::string suffix );
    
    /// \brief Copy constructor
    RateCoeffModel( const RateCoeffModel &r );

    /// \brief Default destructor
    virtual ~RateCoeffModel();

    /// \brief clone function
    ///
    /// \author Rowan J Gollan
    /// \version 24-Feb-2006
    ///
    ///
    virtual RateCoeffModel* clone();

    /// \brief string representation for RateCoeffModel
    virtual std::string str() const;

    virtual std::string latex_str(int nc) const
    { return "No appropriate LaTeX string implemented."; }

    double k() { return k_; };
    void set_k( double k ) { k_ = k; }

    /// \brief Compute the reaction rate coefficient.
    ///
    /// \author Rowan J Gollan
    /// \version 24-Feb-2006
    ///
    virtual void eval( const Gas_data &Q );

    std::string model() const { return model_; }
protected:
    std::string model_;
    std::string suffix_;
    double k_;

};

class GeneralisedArrhenius : public RateCoeffModel {
public:
    /// \brief Default constructor
    GeneralisedArrhenius();
    
    /// \brief Normal constructor
    GeneralisedArrhenius( const std::string suffix, double A, double n, double E_a );
    
    /// \brief Read from a config object
    GeneralisedArrhenius( ConfigParser &cfg, const std::string section, const std::string suffix );

    /// \brief Copy constructor
    GeneralisedArrhenius( const GeneralisedArrhenius &g );

    /// \brief Assignment operator
    GeneralisedArrhenius& operator=(const GeneralisedArrhenius &g);

    /// \brief Default destructor
    virtual ~GeneralisedArrhenius();

    /// \brief clone() function
    GeneralisedArrhenius* clone();

    /// \brief string representation
    std::string str() const;

    std::string latex_str(int nc) const;

    /// \brief Compute the reaction rate coefficient
    virtual void eval( const Gas_data &Q );

    void set_A(double A) { A_ = A; }
    void set_n(double n) { n_ = n; }
    void set_E_a(double E_a) { E_a_ = E_a; }

    double get_A() const { return A_; }
    double get_n() const { return n_; }
    double get_E_a() const { return E_a_; }

protected:
    double A_;
    double n_;
    double E_a_;

};

class HeavisideZeroActivationEnergy : public RateCoeffModel {
public:
    /// \brief Normal constructor
    HeavisideZeroActivationEnergy( const std::string suffix, double alpha, double T_i );
    
    /// \brief Read from a config object
    HeavisideZeroActivationEnergy( ConfigParser &cfg, const std::string section, const std::string suffix );

    /// \brief Copy constructor
    HeavisideZeroActivationEnergy( const HeavisideZeroActivationEnergy &g );

    /// \brief Default destructor
    virtual ~HeavisideZeroActivationEnergy();

    /// \brief clone() function
    HeavisideZeroActivationEnergy* clone();

    /// \brief string representation
    std::string str() const;

    /// \brief Compute the reaction rate coefficient
    void eval( const Gas_data &Q );

private:
    double alpha_;
    double T_i_;

};

class Nonequilibrium_dissociation : public GeneralisedArrhenius {
public:
    /// \brief Normal constructor
    Nonequilibrium_dissociation(const std::string suffix, double A, double n, double E_a,
				int ivib, double alpha, double U);

    /// \brief Read from a config object
    Nonequilibrium_dissociation(ConfigParser &cfg, const std::string section, const std::string suffix);

    /// \brief Copy constructor
    Nonequilibrium_dissociation(const Nonequilibrium_dissociation &n);

    /// \brief Default destructor
    ~Nonequilibrium_dissociation();

    /// \brief clone() function
    Nonequilibrium_dissociation* clone();

    /// \brief string representation
    std::string str() const;

    /// \brief Compute the reaction rate coefficient
    void eval(const Gas_data &Q);
private:
    int ivib_;
    double alpha_;
    double U_;
    double D_;
};


class Nonequilibrium_exchange : public GeneralisedArrhenius {
public:
    /// \brief Normal constructor
    Nonequilibrium_exchange(const std::string suffix, double A, double n, double E_a,
			    int ivib, double alpha, double U, double D);

    /// \brief Read from a config object
    Nonequilibrium_exchange(ConfigParser &cfg, const std::string section, const std::string suffix);

    /// \brief Copy constructor
    Nonequilibrium_exchange(const Nonequilibrium_exchange &n);

    /// \brief Default destructor
    ~Nonequilibrium_exchange();

    /// \brief clone() function
    Nonequilibrium_exchange* clone();

    /// \brief string representation
    std::string str() const;

    /// \brief Compute the reaction rate coefficient
    void eval(const Gas_data &Q);
private:
    int ivib_;
    double alpha_;
    double U_;
    double Aj_;
    double D_;
};

double nonequilibrium_factor(const Gas_data &Q, int ivib, double theta_v,
			     double alpha, double A, double U, double D);

double polyatomic_nonequilibrium_factor(const Gas_data &Q, int ivib,
                                        std::valarray<double> theta_vs,
                                        double alpha, double A, double U, double D);

double vib_partition_function(double theta_v, double Y, double T);

double polyatomic_vib_partition_function(std::valarray<double> theta_vs, double Y, double T);

class ParkNonequilibrium : public RateCoeffModel {
public:
    /// \brief Normal constructor
    ParkNonequilibrium( const std::string suffix, double A, double n, double E_a,
                        int ivib, double s );
    
    /// \brief Read from a config object
    ParkNonequilibrium( ConfigParser &cfg, const std::string section, const std::string suffix );

    /// \brief Copy constructor
    ParkNonequilibrium( const ParkNonequilibrium &g );

    /// \brief Default destructor
    virtual ~ParkNonequilibrium();

    /// \brief clone() function
    ParkNonequilibrium* clone();

    /// \brief string representation
    std::string str() const;

    /// \brief Compute the reaction rate coefficient
    virtual void eval( const Gas_data &Q );

protected:
    double A_;
    double n_;
    double E_a_;
    int ivib_;
    double s_;

};

class RadiativeDecay : public RateCoeffModel {
public:
    /// \brief Normal constructor
    RadiativeDecay( const std::string suffix, double lambda, double tau );
    
    /// \brief Read from a config object
    RadiativeDecay( ConfigParser &cfg, const std::string section, const std::string suffix );

    /// \brief Copy constructor
    RadiativeDecay( const RadiativeDecay &g );

    /// \brief Default destructor
    virtual ~RadiativeDecay();

    /// \brief clone() function
    RadiativeDecay* clone();

    /// \brief string representation
    std::string str() const;

    /// \brief Compute the reaction rate coefficient
    virtual void eval( const Gas_data &Q );

protected:
    double lambda_;
    double tau_;
};

#endif
