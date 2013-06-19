/** \file diatomic_band.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Mar-08: initial implementation
 *           11-Aug-09: Photaura version
 *  \brief Declarations for DiatomicBand class
 *
 **/

#ifndef DIATOMIC_BAND_HH
#define DIATOMIC_BAND_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "spectra_pieces.hh"

#define DIATOMIC_VOIGT_PROFILE_METHOD   1      /* [0] Accurate Whiting expression, [1] approx. Whiting expression    */
#define DIATOMIC_STARK_WIDTH            0      /* [0] Johnston 2006 curve fit, [1] Cowley 1971, [2] Arnold 1979      */
#define DIATOMIC_COLLISION_WIDTH_METHOD 1      /* [0] Johnston 2006, [1] Spradian07                                  */
#define DIATOMIC_NAUTRAL_WIDTH_METHOD   0      /* [0] constant value of Spradian07, [1] classical expression         */
#define DIATOMIC_LINE_TYPE              0      /* Select Voigt [0], Lorentz [1] or Doppler [2]                       */
#define DIATOMIC_LINE_EXTENT            10     /* One-sided line extent in Voigt half-width units                    */
#define DIATOMIC_LINE_POINTS            4      /* Points-per-VHW describing an atomic line                           */
#define DIATOMIC_LIMITED_LINE_EXTENT    1      /* Line extents are unlimited [0] or limited [1]                      */
#define BAND_AVERAGE_FREQUENCY_METHOD   2      /* [0] Exact averaging, [1] Hartung's approximate expression, [2] 0-0 */
#define HUND_AB_DOUBLET_HLF_METHOD      2      /* [0] Whiting's expressions for 2Sigma - 2Pi, [1] Kovacs generic or [2] Hund's case (a) */

#define F_ROT_LINE_LIMIT 1.0e-3

/* Declaration of Diatomic_line class */

class DiatomicLine {
public:
    /// \brief Constructor
    DiatomicLine( double nu_00, double E_u, double m_w, double I, double sigma_nm );
    
    /// \brief Deconstructor
    ~DiatomicLine();
    
public:
    /* Representative frequency for calculating line-widths */
    double nu_00;                        /**< \brief zero-zero transition frequency         */
    
    /* Level data for calculating line-shape */
    double E_u;
    
    /* Species data for calculating line-shape */
    double m_w;
    double I;
    double sigma_nm;
    
    /* Temporary storage of line properties */
    double nu_ul;			 /**< \brief temp centreline transition frequency   */
    double j_ul;			 /**< \brief temp integrated emission coefficient   */
    double kappa_lu;			 /**< \brief temp integrated absorption coefficient */
    
    /* Further data for describing the line-shape */
    double gamma_L;                      /**< \brief Combined Lorentz (half) half-width     */
    double gamma_D;                      /**< \brief Doppler (half) half-width              */
    double gamma_V;                      /**< \brief Voigt (half) half-width                */
    
public:
    /* Initialization functions */
    void set_A_ul();
    void initialize( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av );
    double calculate_Lorentz_width( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av );
    double calculate_collision_width( double T, double p, double N_hvy, double mw_av );
    double calculate_Stark_width( double Te, double N_elecs );
    double calculate_natural_width();
    double calculate_Doppler_width( double T );
    double calculate_Voigt_width();
    
    /* Line-shape functions */
    double get_lorentz_point( double delta_nu );
    double get_doppler_point( double delta_nu );
    double get_voigt_point(double delta_nu );
    
    /* Spectrum generation */
    void calculate_spectrum( CoeffSpectra &X );
    
    /* Access functions */
    std::string line_width_string( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av );
    
};

// Forward declaration of DiatomicElecLev
class DiatomicElecLev;

class DiatomicBand {
public: 
    /// \brief Constructor
    DiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, double Re_vib );
    
    /// \brief Deconstructor
    virtual ~DiatomicBand();
    
public:
    /* Initialization functions */
    virtual void initialize( double T, double Te, double p, double N_hvy = 0.0, double N_elecs = 0.0, double mw_av = 0.0 ) = 0;
    
    virtual std::string line_width_string( double T, double Te, double p, double N_hvy = 0.0, double N_elecs = 0.0, double mw_av = 0.0 ) = 0;
    
    /* Analysis functions */
    virtual void write_lines_to_file( std::string fname ) = 0;
    
    double calculate_average_frequency();
    
    /* Energy functions */
    double calculate_j_ul( double N_u );
    
    virtual void calculate_spectrum( CoeffSpectra &X, double Tv, double Tr ) = 0;
    
public:
    int Vu;
    int Vl;
    
    // Pointers to electronic levels
    DiatomicElecLev * elev_u;
    DiatomicElecLev * elev_l;
    
    double nu_bh;
    double nu_00;
    double nu_av;
    
    double Re_vib;
    double A_ul_av;
    
    int transition_type;
    int tS;
    
    int tJu_max;
    int tJl_max;
    
    int lambda_u;
    int lambda_l;
    
    bool reverse_band;
};

class LBLDiatomicBand : public DiatomicBand {
public:
    /// \brief Constructor
    LBLDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		     double Re_vib, double m_w, double I, double sigma_nm );
    
    /// \brief Deconstructor
    ~LBLDiatomicBand();
    
public:
    void initialize( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av );
    
    std::string line_width_string( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av );
    
    virtual void calculate_spectrum( CoeffSpectra &X, double Tv, double Tr ) = 0;
    
public:
    DiatomicLine * line;
};

class HundALBLDiatomicBand : public LBLDiatomicBand {
public:
    /// \brief Constructor
    HundALBLDiatomicBand( int Vu, int Vl,
    			    DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		            double Re_vib, double m_w, double I, double sigma_nm );
    
    /// \brief Deconstructor
    ~HundALBLDiatomicBand() {}
    
public:
    void write_lines_to_file( std::string fname );
    
    void calculate_spectrum( CoeffSpectra &X, double Tv, double Tr );
    
private:
    void calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, double &j_ul_00 );
    
    double get_HLF( int tJu, int tJl);
};

class SigmaTripletLBLDiatomicBand : public LBLDiatomicBand {
public:
    /// \brief Constructor
    SigmaTripletLBLDiatomicBand( int Vu, int Vl,
    				  DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		                  double Re_vib, double m_w, double I, double sigma_nm );
    
    /// \brief Deconstructor
    ~SigmaTripletLBLDiatomicBand() {}
    
public:
    void write_lines_to_file( std::string fname );
    
    void calculate_spectrum( CoeffSpectra &X, double Tv, double Tr );
    
private:
    void calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, double &j_ul_00 );
    
    double get_HLF( int tJu, int tJl  );
};

class HundBDoubletLBLDiatomicBand : public LBLDiatomicBand {
public:
    /// \brief Constructor
    HundBDoubletLBLDiatomicBand( int Vu, int Vl,
    				     DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		                     double Re_vib, double m_w, double I, double sigma_nm );
    
    /// \brief Deconstructor
    ~HundBDoubletLBLDiatomicBand() {}
    
public:
    void write_lines_to_file( std::string fname );
    
    void calculate_spectrum( CoeffSpectra &X, double Tv, double Tr );
    
private:
    void calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, int tSig_u, int tSig_l, double &j_ul_00 );
    
    double get_HLF( int tJu, int tJl, int tSig_u, int tSig_l );
};

class HundABDoubletLBLDiatomicBand : public LBLDiatomicBand {
public:
    /// \brief Constructor
    HundABDoubletLBLDiatomicBand( int Vu, int Vl,
    				  DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		                  double Re_vib, double m_w, double I, double sigma_nm );
    
    /// \brief Deconstructor
    ~HundABDoubletLBLDiatomicBand() {}
    
public:
    void write_lines_to_file( std::string fname );
    
    void calculate_spectrum( CoeffSpectra &X, double Tv, double Tr );
    
private:
    void calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, int tSig_u, int tSig_l, double &j_ul_00 );
    
    double get_HLF( int tJu, int tJl, int tSig_u, int tSig_l );
    
    double u_plus( double lambda, double A, double B_v, double J );
    
    double u_minus( double lambda, double A, double B_v, double J );
    
    double C_plus( double lambda, double A, double B_v, double J );
    
    double C_minus( double lambda, double A, double B_v, double J );
};

class HundABTripletLBLDiatomicBand : public LBLDiatomicBand {
public:
    /// \brief Constructor
    HundABTripletLBLDiatomicBand( int Vu, int Vl,
    				  DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		                  double Re_vib, double m_w, double I, double sigma_nm );
    
    /// \brief Deconstructor
    ~HundABTripletLBLDiatomicBand() {}
    
public:
    void write_lines_to_file( std::string fname );
    
    void calculate_spectrum( CoeffSpectra &X, double Tv, double Tr );
    
private:
    void calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, int tSig_u, int tSig_l, double &j_ul_00 );
    
    double get_HLF( int tJu, int tJl, int tSig_u, int tSig_l );
    
    double u_1_plus( double lambda, double A, double B_v, double J );
    
    double u_1_minus( double lambda, double A, double B_v, double J );
    
    double u_3_plus( double lambda, double A, double B_v, double J );
    
    double u_3_minus( double lambda, double A, double B_v, double J );
    
    double C_1( double lambda, double A, double B_v, double J );
    
    double C_2( double lambda, double A, double B_v, double J );
    
    double C_3( double lambda, double A, double B_v, double J );
};

class SRBDiatomicBand : public DiatomicBand {
public:
    /// \brief Constructor
    SRBDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		     double Re_vib );
    
    /// \brief Deconstructor
    ~SRBDiatomicBand();
    
public:
    void initialize( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av );
    
    std::string line_width_string( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av);
    
    void write_lines_to_file( std::string fname );
    
    void calculate_spectrum( CoeffSpectra &X, double Tv, double Tr );
    
public:
    double SB_constA;
    double SB_constB;
    double SB_constC;
};

#endif

