/** \file stark_width_models.hh
 *  \ingroup radiation
 *
 *  \brief Classes for describing Stark broadening of lines
 *
 *  \author Daniel F. Potter
 *  \version 16-Sept-13: initial implementation
 *
 **/

#ifndef STARK_WIDTH_MODELS_HH
#define STARK_WIDTH_MODELS_HH

class StarkWidthModel {
public:
    /// \brief Constructor
    StarkWidthModel( std::string name, double n );

    /// \brief Deconstructor
    virtual ~StarkWidthModel();

public:
    /// \brief Evaluate the Stark width for the given conditions
    virtual double eval( double T_e, double N_e  ) = 0;

protected:
    /// \brief Model name
    std::string name;

    /// \brief Exponential factor
    double n;
};

class ApproxStarkWidth : public StarkWidthModel {
public:
    /// \brief Constructor
    ApproxStarkWidth( double n, double constA, double constB, double I, double E_u );

    /// \brief Deconstructor
    virtual ~ApproxStarkWidth();

public:
    /// \brief Evaluate the Stark width for the given conditions
    double eval( double T_e, double N_e  );

private:
    /// \brief Parameter for determining the reference Stark width
    double constA;

    /// \brief Parameter for determining the reference Stark width
    double constB;

    /// \brief Ionisation energy in cm-1
    double I_icm;

    /// \brief Upper level energy in cm-1
    double E_u_icm;
};

class GriemStarkWidth : public StarkWidthModel {
public:
    /// \brief Constructor
    GriemStarkWidth( double n, double gamma_s0 );

    /// \brief Deconstructor
    virtual ~GriemStarkWidth();

public:
    /// \brief Evaluate the Stark width for the given conditions
    double eval( double T_e, double N_e  );

private:
    /// \brief Reference Stark width in Hz
    double gamma_S0;
};

#endif
