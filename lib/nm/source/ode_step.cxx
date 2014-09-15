/*  \file ode_step.cxx
 *  \ingroup nm
 *  \brief Definitions for the OdeStep class and derivatives.
 *  \author Rowan J Gollan
 *  \version 20-Feb-2006
 **/

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "ode_system.hh"
#include "ode_step.hh"
using namespace std;

/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
/// \param name: name of  the ODE stepping algorithm
///
OdeStep::OdeStep( const string name, int ndim )
    : name_( name ), ndim_( ndim ) {}

/// \brief Copy constructor
/// 
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
OdeStep::OdeStep( const OdeStep &o )
    : name_( o.name() ), ndim_( o.ndim() ) {}

/// \brief Default destuctor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
OdeStep::~OdeStep() {}

/// \brief string representation of the step
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
string OdeStep::str() const
{
    ostringstream ost;
    ost << "OdeStep(\n"
	<< "   name= " << name_ << endl
	<< ")\n";
    return ost.str();
}

/// \brief Overloaded ostream operator
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
ostream& operator<<( ostream &os, const OdeStep &o )
{
    os << o.str();
    return os;
}

//-----------------------------------------------------------------------//
// Derived OdeSteps: actual ODE stepping algorithms                      //
//-----------------------------------------------------------------------//

/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
/// \param name: name of  the ODE stepping algorithm
///
EulerStep::EulerStep( const string name, int ndim )
    : OdeStep( name, ndim )
{
    y_dot_.resize(ndim);
}

/// \brief Copy constructor
/// 
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
EulerStep::EulerStep( const EulerStep &o )
    : OdeStep( o.name(), o.ndim() )
{
    y_dot_.resize(o.ndim());
}

/// \brief Default destuctor
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
EulerStep::~EulerStep() {}

/// \brief Clone of Euler step
/// 
/// \author Rowan J Gollan
/// \version 20-Feb-2006
///
EulerStep* EulerStep::clone()
{
    return new EulerStep( *this );
}

/// \brief Takes an Euler step
///
/// \author Rowan J Gollan
/// \version 20-Feb-2006
/// 
/// \param ode  : an OdeSystem object providing the eval function
/// \param yin  : the input vector
/// \param yout : the output vector
/// \param h    : the stepsize - this will be unchanged as the
///               Euler method provides no means of predicting
///               a new stepsize
/// \returns    : true if successful, false otherwise.
///
/// This function implements the time-honoured Euler algorithm
/// for numerically advancing a first-order OdeSystem.
///
bool EulerStep::advance( OdeSystem &ode, const vector<double> &yin,
			 vector<double> &yout, double *h )
{
    int ndim = yin.size();

    ode.eval( yin, y_dot_ );
    ode.called_at_least_once = true;


    for( int i = 0; i < ndim; ++i ) {
	yout[i] = yin[i] + (*h) * y_dot_[i];
    }

    return true;
}

/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
ModEulerStep::ModEulerStep( const string name, int ndim )
    : OdeStep( name, ndim )
{
    y_nplus1_.resize(ndim);
    dy_n_.resize(ndim);
    dy_nplus1_.resize(ndim);
}

/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
ModEulerStep::ModEulerStep( const ModEulerStep &m )
    : OdeStep( m.name_, m.ndim_ )
{

    y_nplus1_.resize(ndim_);
    dy_n_.resize(ndim_);
    dy_nplus1_.resize(ndim_);

}

/// \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
ModEulerStep::~ModEulerStep() {}

/// \brief Clone of the ModEuler object
/// 
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
ModEulerStep* ModEulerStep::clone()
{
    return new ModEulerStep( *this );
}


/// \brief Takes a modified Euler step
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
/// \param ode  : an OdeSystem object providing the eval function
/// \param yin  : the input vector
/// \param yout : the output vector
/// \param h    : the stepsize - this will be unchanged as the
///               modified Euler method provides no means of predicting
///               a new stepsize
/// \returns    : true if successful, false otherwise.
///
/// This function implements the modified Euler algorithm
/// for numerically advancing a first-order OdeSystem.
///
bool ModEulerStep::advance( OdeSystem &ode, const vector<double> &yin, 
			vector<double> &yout, double *h )
{
    int ndim = yin.size();

    ode.eval( yin, dy_n_ );
    ode.called_at_least_once = true;

    // 1. Evaluate y_n+1
    for( int i = 0; i < ndim; ++i ) {
	y_nplus1_[i] = yin[i] + (*h) * dy_n_[i];
    }
    
    // 2. Now find the estimate for y'(n+1)
    ode.eval( y_nplus1_, dy_nplus1_ );
    
    // 3. finally evaluate yout
    for( int i = 0; i < ndim; ++i ) {
	yout[i] = yin[i] + (*h) * 0.5 * ( dy_n_[i] + dy_nplus1_[i] );
    }

    return true;
}

/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
RKFStep::RKFStep( const string name, int ndim, double tol )
    : OdeStep( name, ndim ), tol_( tol )
{
    if( tol_ <= 0.0 ) {
	cout << "RKFStep::RKFStep() WARNING: input tolerance is less than 0.0\n";
	cout << "                   _tol has been set to 1.0e-8.\n";
	tol_ = 1.0e-8;
    }
		
    tmp_.resize(ndim);
    yerr_.resize(ndim);
    k1_.resize(ndim);
    k2_.resize(ndim);
    k3_.resize(ndim);
    k4_.resize(ndim);
    k5_.resize(ndim);
    k6_.resize(ndim);

    set_constants();
    
}

/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
RKFStep::RKFStep( const RKFStep &r )
    : OdeStep( r.name_, r.ndim_ ), tol_( r.tol_ )
{
    tmp_.resize(ndim_);
    yerr_.resize(ndim_);
    k1_.resize(ndim_);
    k2_.resize(ndim_);
    k3_.resize(ndim_);
    k4_.resize(ndim_);
    k5_.resize(ndim_);
    k6_.resize(ndim_);

    set_constants();

}

void
RKFStep::
set_constants()
{
    // Constants from
    // Cash and Karp (1990)
    // A Variable Order Runge-Kutta Method for
    // Initial Value Problems woth Rapidly Varying
    // Right-Hand Sides
    // ACM Transactions on Mathematical Software, 16:3 pp 201--222
    
    a21_ = 1.0/5.0;
    a31_ = 3.0/40.0;
    a32_ = 9.0/40.0;
    a41_ = 3.0/10.0;
    a42_ = -9.0/10.0;
    a43_ = 6.0/5.0;
    a51_ = -11.0/54.0;
    a52_ = 5.0/2.0;
    a53_ = -70.0/27.0;
    a54_ = 35.0/27.0;
    a61_ = 1631.0/55296.0;
    a62_ = 175.0/512.0;
    a63_ = 575.0/13824.0;
    a64_ = 44275.0/110592.0;
    a65_ = 253.0/4096.0;

    // For autonomous systems (which we treat here)
    // these are not needed.
//     const double c2_ = 1.0/5.0;
//     const double c3_ = 3.0/10.0;
//     const double c4_ = 3.0/5.0;
//     const double c5_ = 1.0;
//     const double c6_ = 7.0/8.0;

    b51_ = 37.0/378.0;
    b53_ = 250.0/621.0;
    b54_ = 125.0/594.0;
    b56_ = 512.0/1771.0;
    
    b41_ = 2825.0/27648.0;
    b43_ = 18575.0/48384.0;
    b44_ = 13525.0/55296.0;
    b45_ = 277.0/14336.0;
    b46_ = 1.0/4.0;

}

/// \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
RKFStep::~RKFStep() {}

/// \brief Clone of the RKFStep object
/// 
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
RKFStep* RKFStep::clone()
{
    return new RKFStep( *this );
}


/// \brief Takes a step using the Runge-Kutta-Fehlberg method
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
/// \param ode  : an OdeSystem object providing the eval function
/// \param yin  : the input vector
/// \param yout : the output vector
/// \param h    : the stepsize - this will be changed at the end of the
///               step as a suggestion for further stepping
/// \returns    : true if successful, false otherwise.
///
/// This function implements the Runge-Kutta-Fehlberg algorithm
/// for numerically advancing a first-order OdeSystem.
///
/// The constants for integration come from Cash and Karp (1990)
/// and the adaptive step size control is based on Press et al. (2007)
///
/// References:
/// 
/// Cash, J. R. and Karp, A. H. (1990)
/// A Variable Order Runge-Kutta Method for Initial Value
/// Problems with Rapidly Varying Right-Hand Sides
/// ACM Transactions on Mathematical Software, 16:3, pp. 201--222
///
/// Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007)
/// Numerical Recipes: The Art of Scientific Computing, Third Edition
/// Cambridge University Press, New York, USA
///
bool RKFStep::advance( OdeSystem &ode, const vector<double> &yin, 
		       vector<double> &yout, double *h )
{
		
    ode.eval(yin, k1_);
    ode.called_at_least_once = true;

    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a21_*k1_[i]);
    }
    ode.eval(tmp_, k2_);

    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a31_*k1_[i] + a32_*k2_[i]);
    }
    ode.eval(tmp_, k3_);

    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a41_*k1_[i] + a42_*k2_[i] + a43_*k3_[i]);
    }
    ode.eval(tmp_, k4_);

    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a51_*k1_[i] + a52_*k2_[i] + a53_*k3_[i] + 
				 a54_*k4_[i]);
    }
    ode.eval(tmp_, k5_);

    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a61_*k1_[i] + a62_*k2_[i] + a63_*k3_[i] +
				 a64_*k4_[i] + a65_*k5_[i]);
    }
    ode.eval(tmp_, k6_);

    for ( int i = 0; i < ndim_; ++i ) {
	yout[i] = yin[i] + (*h)*(b51_*k1_[i] + b53_*k3_[i] + b54_*k4_[i] +
				 b56_*k6_[i]);
    }

    for ( int i = 0; i < ndim_; ++i ) {
	// 5th-order estimate minus 4th-order estimate
	yerr_[i] = yout[i] - (yin[i] + (*h)*(b41_*k1_[i] + b43_*k3_[i] + b44_*k4_[i] +
					     b45_*k5_[i] + b46_*k6_[i]));
    }

    // Compute error using tol as atol and rtol as suggested in
    // Press et al. (2007)
    double err = 0.0;
    double sk = 0.0;
    double atol = tol_;
    double rtol = tol_;

    for ( int i = 0; i < ndim_; ++i ) {
	sk = atol + rtol*max(fabs(yin[i]), fabs(yout[i]));
	err += (yerr_[i]/sk)*(yerr_[i]/sk);
    }
    err = sqrt(err/ndim_);

    // Now use error as an estimate for new step size
    double scale = 0.0;
    const double maxscale = 10.0;
    const double minscale = 0.2;
    const double safe = 0.9;
    const double alpha = 0.2;
    if ( err <= 1.0 ) {
	// A successful step
	if( err == 0.0 ) {
	    scale = maxscale;
	}
	else {
	    // We are NOT using the PI control version as
	    // given by Press et al. as I do not want to store
	    // the old error from previous integration steps.
	    scale = safe * pow(err, -alpha);
	    if ( scale < minscale )
		scale = minscale;
	    if ( scale > maxscale )
		scale = maxscale;
	}

	*h *= scale;
	return true;
    }
    // else, failed step
    scale = max(safe*pow(err, -alpha), minscale);
    *h *= scale;
    
    return false;
}

/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 17-Apr-2010
///
DP853Step::DP853Step(const string name, int ndim, double tol)
    : OdeStep( name, ndim ), tol_( tol )
{
    if( tol_ <= 0.0 ) {
	cout << "DP853Step::DP853Step() WARNING: input tolerance is less than 0.0\n";
	cout << "                                tol_ has been set to 1.0e-8.\n";
	tol_ = 1.0e-8;
    }
		
    tmp_.resize(ndim);
    yerr_.resize(ndim);
    yerr2_.resize(ndim);
    k1_.resize(ndim);
    k2_.resize(ndim);
    k3_.resize(ndim);
    k4_.resize(ndim);
    k5_.resize(ndim);
    k6_.resize(ndim);
    k7_.resize(ndim);
    k8_.resize(ndim);
    k9_.resize(ndim);
    k10_.resize(ndim);

    set_constants();
    
}

/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 17-Apr-2010
///
DP853Step::DP853Step(const DP853Step &d)
    : OdeStep(d.name_, d.ndim_), tol_(d.tol_)
{
    tmp_.resize(ndim_);
    yerr_.resize(ndim_);
    yerr2_.resize(ndim_);
    k1_.resize(ndim_);
    k2_.resize(ndim_);
    k3_.resize(ndim_);
    k4_.resize(ndim_);
    k5_.resize(ndim_);
    k6_.resize(ndim_);
    k7_.resize(ndim_);
    k8_.resize(ndim_);
    k9_.resize(ndim_);
    k10_.resize(ndim_);

    set_constants();
    
}

/// \brief Set constants for Dormand-Prince 8(53) method
///
/// \author Rowan J Gollan
/// \version 17-Apr-2010
///
void
DP853Step::
set_constants()
{
    // Copy-n-pasted from the Numerical Recipes in C
    // PDF file available on web for this method.
    b1_ = 5.42937341165687622380535766363e-2;
    b6_ = 4.45031289275240888144113950566e0;
    b7_ = 1.89151789931450038304281599044e0;
    b8_ = -5.8012039600105847814672114227e0;
    b9_ =   3.1116436695781989440891606237e-1;
    b10_ = -1.52160949662516078556178806805e-1;
    b11_ = 2.01365400804030348374776537501e-1;
    b12_ = 4.47106157277725905176885569043e-2;

    bhh1_ = 0.244094488188976377952755905512e+00;
    bhh2_ = 0.733846688281611857341361741547e+00;
    bhh3_ = 0.220588235294117647058823529412e-01;

    er1_  =  0.1312004499419488073250102996e-01;
    er6_  = -0.1225156446376204440720569753e+01;
    er7_  = -0.4957589496572501915214079952e+00;
    er8_  =  0.1664377182454986536961530415e+01;
    er9_  = -0.3503288487499736816886487290e+00;
    er10_ =  0.3341791187130174790297318841e+00;
    er11_ =  0.8192320648511571246570742613e-01;
    er12_ = -0.2235530786388629525884427845e-01;

    a21_ =  5.26001519587677318785587544488e-2;
    a31_ =  1.97250569845378994544595329183e-2;
    a32_ =  5.91751709536136983633785987549e-2;
    a41_ =  2.95875854768068491816892993775e-2;
    a43_ =  8.87627564304205475450678981324e-2;
    a51_ =  2.41365134159266685502369798665e-1;
    a53_ = -8.84549479328286085344864962717e-1;
    a54_ =  9.24834003261792003115737966543e-1;
    a61_ =  3.7037037037037037037037037037e-2;
    a64_ =  1.70828608729473871279604482173e-1;
    a65_ =  1.25467687566822425016691814123e-1;
    a71_ =  3.7109375e-2;
    a74_ =  1.70252211019544039314978060272e-1;
    a75_ =  6.02165389804559606850219397283e-2;
    a76_ = -1.7578125e-2;
    a81_ =   3.70920001185047927108779319836e-2;
    a84_ =   1.70383925712239993810214054705e-1;
    a85_ =   1.07262030446373284651809199168e-1;
    a86_ =  -1.53194377486244017527936158236e-2;
    a87_ =   8.27378916381402288758473766002e-3;
    a91_ =   6.24110958716075717114429577812e-1;
    a94_ =  -3.36089262944694129406857109825e0;
    a95_ =  -8.68219346841726006818189891453e-1;
    a96_ =   2.75920996994467083049415600797e1;
    a97_ =   2.01540675504778934086186788979e1;
    a98_ =  -4.34898841810699588477366255144e1;
    a101_ =  4.77662536438264365890433908527e-1;
    a104_ = -2.48811461997166764192642586468e0;
    a105_ = -5.90290826836842996371446475743e-1;
    a106_ =  2.12300514481811942347288949897e1;
    a107_ =  1.52792336328824235832596922938e1;
    a108_ = -3.32882109689848629194453265587e1;
    a109_ = -2.03312017085086261358222928593e-2;
    a111_ =  -9.3714243008598732571704021658e-1;
    a114_ =   5.18637242884406370830023853209e0;
    a115_ =   1.09143734899672957818500254654e0;
    a116_ =  -8.14978701074692612513997267357e0;
    a117_ =  -1.85200656599969598641566180701e1;
    a118_ =   2.27394870993505042818970056734e1;
    a119_ =   2.49360555267965238987089396762e0;
    a1110_ = -3.0467644718982195003823669022e0;
    a121_ =   2.27331014751653820792359768449e0;
    a124_ =  -1.05344954667372501984066689879e1;
    a125_ =  -2.00087205822486249909675718444e0;
    a126_ =  -1.79589318631187989172765950534e1;
    a127_ =   2.79488845294199600508499808837e1;
    a128_ =  -2.85899827713502369474065508674e0;
    a129_ =  -8.87285693353062954433549289258e0;
    a1210_ =  1.23605671757943030647266201528e1;
    a1211_ =  6.43392746015763530355970484046e-1;
    a141_ = 5.61675022830479523392909219681e-2;
    a147_ = 2.53500210216624811088794765333e-1;
    a148_ = -2.46239037470802489917441475441e-1;
    a149_ = -1.24191423263816360469010140626e-1;
    a1410_ = 1.5329179827876569731206322685e-1;
    a1411_ = 8.20105229563468988491666602057e-3;
    a1412_ = 7.56789766054569976138603589584e-3;
    a1413_ = -8.298e-3;
    a151_ = 3.18346481635021405060768473261e-2;
    a156_ = 2.83009096723667755288322961402e-2;
    a157_ = 5.35419883074385676223797384372e-2;
    a158_ = -5.49237485713909884646569340306e-2;
    a1511_ = -1.08347328697249322858509316994e-4;
    a1512_ = 3.82571090835658412954920192323e-4;
    a1513_ = -3.40465008687404560802977114492e-4;
    a1514_ = 1.41312443674632500278074618366e-1;
    a161_ = -4.28896301583791923408573538692e-1;
    a166_ = -4.69762141536116384314449447206e0;
    a167_ = 7.68342119606259904184240953878e0;
    a168_ = 4.06898981839711007970213554331e0;
    a169_ = 3.56727187455281109270669543021e-1;
    a1613_ = -1.39902416515901462129418009734e-3;
    a1614_ = 2.9475147891527723389556272149e0;
    a1615_ = -9.15095847217987001081870187138e0;

    d41_  = -0.84289382761090128651353491142e+01;
    d46_  =  0.56671495351937776962531783590e+00;
    d47_  = -0.30689499459498916912797304727e+01;
    d48_  =  0.23846676565120698287728149680e+01;
    d49_  =  0.21170345824450282767155149946e+01;
    d410_ = -0.87139158377797299206789907490e+00;
    d411_ =  0.22404374302607882758541771650e+01;
    d412_ =  0.63157877876946881815570249290e+00;
    d413_ = -0.88990336451333310820698117400e-01;
    d414_ =  0.18148505520854727256656404962e+02;
    d415_ = -0.91946323924783554000451984436e+01;
    d416_ = -0.44360363875948939664310572000e+01;
    d51_  =  0.10427508642579134603413151009e+02;
    d56_  =  0.24228349177525818288430175319e+03;
    d57_  =  0.16520045171727028198505394887e+03;
    d58_  = -0.37454675472269020279518312152e+03;
    d59_  = -0.22113666853125306036270938578e+02;
    d510_ =  0.77334326684722638389603898808e+01;
    d511_ = -0.30674084731089398182061213626e+02;
    d512_ = -0.93321305264302278729567221706e+01;
    d513_ =  0.15697238121770843886131091075e+02;
    d514_ = -0.31139403219565177677282850411e+02;
    d515_ = -0.93529243588444783865713862664e+01;
    d516_ =  0.35816841486394083752465898540e+02;
    d61_ = 0.19985053242002433820987653617e+02;
    d66_ = -0.38703730874935176555105901742e+03;
    d67_ = -0.18917813819516756882830838328e+03;
    d68_ = 0.52780815920542364900561016686e+03;
    d69_ = -0.11573902539959630126141871134e+02;
    d610_ = 0.68812326946963000169666922661e+01;
    d611_ = -0.10006050966910838403183860980e+01;
    d612_ = 0.77771377980534432092869265740e+00;
    d613_ = -0.27782057523535084065932004339e+01;
    d614_ = -0.60196695231264120758267380846e+02;
    d615_ = 0.84320405506677161018159903784e+02;
    d616_ = 0.11992291136182789328035130030e+02;
    d71_  = -0.25693933462703749003312586129e+02;
    d76_  = -0.15418974869023643374053993627e+03;
    d77_  = -0.23152937917604549567536039109e+03;
    d78_  =  0.35763911791061412378285349910e+03;
    d79_  =  0.93405324183624310003907691704e+02;
    d710_ = -0.37458323136451633156875139351e+02;
    d711_ =  0.10409964950896230045147246184e+03;
    d712_ =  0.29840293426660503123344363579e+02;
    d713_ = -0.43533456590011143754432175058e+02;
    d714_ =  0.96324553959188282948394950600e+02;
    d715_ = -0.39177261675615439165231486172e+02;
    d716_ = -0.14972683625798562581422125276e+03;

}

/// \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 17-Apr-2010
///
DP853Step::~DP853Step() {}

/// \brief Clone of the DP853Step object
/// 
/// \author Rowan J Gollan
/// \version 17-Apr-2010
///
DP853Step* DP853Step::clone()
{
    return new DP853Step( *this );
}

bool
DP853Step::advance(OdeSystem &ode, const vector<double> &yin, 
		   vector<double> &yout, double *h )
{
    ode.eval(yin, k1_);
    ode.called_at_least_once = true;
    
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a21_*k1_[i]);
    }

    ode.eval(tmp_, k2_);
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a31_*k1_[i]+a32_*k2_[i]);
    }

    ode.eval(tmp_, k3_);
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a41_*k1_[i]+a43_*k3_[i]);
    }

    ode.eval(tmp_, k4_);
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a51_*k1_[i]+a53_*k3_[i]+a54_*k4_[i]);
    }

    ode.eval(tmp_, k5_);
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a61_*k1_[i]+a64_*k4_[i]+a65_*k5_[i]);
    }

    ode.eval(tmp_, k6_);
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a71_*k1_[i]+a74_*k4_[i]+a75_*k5_[i]+a76_*k6_[i]);
    }

    ode.eval(tmp_, k7_);
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a81_*k1_[i]+a84_*k4_[i]+a85_*k5_[i]+a86_*k6_[i]+a87_*k7_[i]);
    }

    ode.eval(tmp_, k8_);
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a91_*k1_[i]+a94_*k4_[i]+a95_*k5_[i]+
				 a96_*k6_[i]+a97_*k7_[i]+a98_*k8_[i]);
    }

    ode.eval(tmp_, k9_);
    for ( int i = 0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a101_*k1_[i]+a104_*k4_[i]+a105_*k5_[i]+
				 a106_*k6_[i]+a107_*k7_[i]+a108_*k8_[i]+a109_*k9_[i]);
    }

    ode.eval(tmp_, k10_);

    for(int i=0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a111_*k1_[i]+a114_*k4_[i]+a115_*k5_[i]+a116_*k6_[i]+
				 a117_*k7_[i]+a118_*k8_[i]+a119_*k9_[i]+a1110_*k10_[i]);
    }

    ode.eval(tmp_, k2_);
    for( int i=0; i < ndim_; ++i ) {
	tmp_[i] = yin[i] + (*h)*(a121_*k1_[i]+a124_*k4_[i]+a125_*k5_[i]+a126_*k6_[i]+
				 a127_*k7_[i]+a128_*k8_[i]+a129_*k9_[i]+a1210_*k10_[i]+a1211_*k2_[i]);
    }

    ode.eval(tmp_, k3_);
    for ( int i = 0; i < ndim_; ++i ) {
	k4_[i] = b1_*k1_[i]+b6_*k6_[i]+b7_*k7_[i]+b8_*k8_[i]+b9_*k9_[i]+b10_*k10_[i]+b11_*k2_[i]+b12_*k3_[i];
	yout[i] = yin[i] + (*h)*k4_[i];
    }

    for ( int i = 0; i < ndim_; ++i ) {
	yerr_[i] = k4_[i]-bhh1_*k1_[i]-bhh2_*k9_[i]-bhh3_*k3_[i];
	yerr2_[i] = er1_*k1_[i]+er6_*k6_[i]+er7_*k7_[i]+er8_*k8_[i]+er9_*k9_[i]+er10_*k10_[i]+er11_*k2_[i]+er12_*k3_[i];
    }

    double err = 0.0;
    double err2 = 0.0;
    double sk, denom;

    double atol = tol_;
    double rtol = tol_;

    for ( int i = 0; i < ndim_; ++i ) {
	sk = atol + rtol*max(fabs(yin[i]),fabs(yout[i]));
	err2 += (yerr_[i]/sk)*(yerr_[i]/sk);
	err += (yerr2_[i]/sk)*(yerr2_[i]/sk);
    }

    denom = err + 0.01*err2;
    if ( denom <= 0.0 ) {
	denom = 1.0;
    }
   
    double error = fabs(*h)*err*sqrt(1.0/(ndim_*denom));

    // And use error to adapt stepsize
    const double alpha = 1.0/8.0;
    const double safe = 0.9;
    const double minscale = 0.333;
    const double maxscale = 6.0;
    double scale;

    if ( error <= 1.0 ) {
	if ( error == 0.0 ) {
	    scale = maxscale;
	}
	else { 
	    scale = safe*pow(error, -alpha);
	    if ( scale < minscale ) {
		scale = minscale;
	    }
	    if ( scale > maxscale ) {
		scale = maxscale;
	    }
	}

	*h *= scale;
	return true;
    }

    // else in case of failure
    scale = max(safe*pow(error, -alpha), minscale);
    *h *= scale;
    return false;
}


/// \brief Normal constructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
QssStep::QssStep( const string name, int ndim, int max_correctors,
		  double qss_eps1, double c, double delta )

    : OdeStep( name, ndim ), max_correctors_( max_correctors ),
      qss_eps1_( qss_eps1 ), c_( c ), delta_( delta )
{
    if( (qss_eps1_ <= 0.0) || (c_ < 1.0) ) {
	cout << "QssStep::QssStep() WARNING: one of the tolerance input parameters was invalid.\n";
	cout << "                   qss_eps1_= " << qss_eps1_ << " c_= " << c_ << endl;
	cout << "                   The values have been set to\n";
	cout << "                   qss_eps1_ = 1.1e-5, c_ = 1.1 --> qss_eps2_ = 1.0e-5\n";
	qss_eps1_ = 1.1e-5;
	c_ = 1.1;
	qss_eps2_ = 1.0e-5;

    }
    else
	qss_eps2_ = qss_eps1_ / c_;

    p0_.resize(ndim);
    q0_.resize(ndim);
    L0_.resize(ndim);
    L_p_.resize(ndim);
    a0_.resize(ndim);
    y_p_.resize(ndim);
    y_p1_.resize(ndim);
    y_c_.resize(ndim);
    p_p_.resize(ndim);
    p_bar_.resize(ndim);
    a_bar_.resize(ndim);
    q_p_.resize(ndim);
    q_til_.resize(ndim);
}

/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
QssStep::QssStep( const QssStep &q )
    : OdeStep( q.name_, q.ndim_ ), max_correctors_( q.max_correctors_ ),
      qss_eps1_( q.qss_eps1_ ), c_( q.c_ ), qss_eps2_( q.qss_eps2_ ), delta_(q.delta_)
{
    p0_.resize(ndim_);
    q0_.resize(ndim_);
    L0_.resize(ndim_);
    L_p_.resize(ndim_);
    a0_.resize(ndim_);
    y_p_.resize(ndim_);
    y_p1_.resize(ndim_);
    y_c_.resize(ndim_);
    p_p_.resize(ndim_);
    p_bar_.resize(ndim_);
    a_bar_.resize(ndim_);
    q_p_.resize(ndim_);
    q_til_.resize(ndim_);

}

/// \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
QssStep::~QssStep() {}

/// \brief Clone of the QssStep object
/// 
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
QssStep* QssStep::clone()
{
    return new QssStep( *this );
}


const double ZERO_EPS = 1.0e-50;


/// \brief Takes a step using the \f$\alpha\f%-QSS method
///
/// \author Rowan J Gollan
/// \version 21-Feb-2006
///
/// PLEASE PUT SOME COMMENTS HERE!
///
bool QssStep::advance( OdeSystem &ode, const vector<double> &yin, 
		       vector<double> &yout, double *h )
{
    
    ode.eval_split( yin, q0_, L0_ );
    ode.called_at_least_once = true;

    // --- Predictor Step --- //

    p_on_y(L0_, yin, p0_);
    alpha( *h, p0_, a0_);

    for (int i = 0; i < ndim_; ++i) {
	y_p_[i] = yin[i] + ( ((*h) * (q0_[i] - L0_[i])) / (1.0 + a0_[i] * (*h) * p0_[i]));
	// Save initial predictor for use in error criterion
	y_p1_[i] = y_p_[i];
    }

    // --- Corrector step(s) --- //
    for (int corr = 0; corr <= max_correctors_; ++corr) {
	ode.eval_split( y_p_, q_p_, L_p_ );
	p_on_y( L_p_, y_p_, p_p_);
	p_bar( p0_, p_p_, p_bar_);
	alpha( *h, p_bar_, a_bar_);
	q_tilde(q0_, q_p_, a_bar_, q_til_);

	/* Actual corrector */
	for( int i = 0; i < ndim_; ++i) {
	    y_c_[i] = yin[i] + (( (*h) * ( q_til_[i] - p_bar_[i] * yin[i] ) ) / ( 1.0 + a_bar_[i] * (*h) * p_bar_[i] ));
	}

	bool converged = test_converged(y_c_, y_p1_);

	if (converged) {
	    copy_vector(y_c_, yout);
	    *h = step_suggest( *h, y_c_, y_p1_ );
	    return true;
	}
	copy_vector(y_c_, y_p_);
    }
    
    *h = step_suggest( *h, y_c_, y_p1_ );
    
    return false;

}

void QssStep::p_on_y( const vector<double> &p, const vector<double> &y,
		      vector<double> &p_y )

{

    for( int i = 0; i < ndim_; ++i ) {
	p_y[i] = p[i] / (y[i] + ZERO_EPS);
    }
    return;
}

void QssStep::alpha( double h, const vector<double> &p, vector<double> &a )
{
    for( int i = 0; i < ndim_; ++i ) {
	double r = 1.0 / ( p[i]*h + ZERO_EPS );
	double numer = 180.0 * r * r * r + 60.0 * r * r + 11.0 * r + 1;
	double denom = 360.0 * r * r * r + 60.0 * r * r + 12.0 * r + 1;
	a[i] = numer/denom;
    }
    return;
}

void QssStep::p_bar( const vector<double> &p0, const vector<double> &p_p, 
		     vector<double> &p_bar )
{
    for( int i = 0; i < ndim_; ++i ) {
	p_bar[i] = 0.5 * ( p0[i] + p_p[i] );
    }
    return;
}

void QssStep::q_tilde( const vector<double> &q0, const vector<double> &q_p,
		       const vector<double> &a_bar, vector<double> &q_til )
{
    for( int i = 0; i < ndim_; ++i ) {
	q_til[i] = a_bar[i]*q_p[i] + (1.0 - a_bar[i])*q0[i];

    }
    return;
}

bool QssStep::test_converged( const vector<double> &y_c, const vector<double> &y_p )
{
    int flag = 0;
    double test = 0.0;
    for( int i = 0; i < ndim_; ++i ) {
	if( y_c[i] < ZERO_EPS )
	    continue;
	test = fabs(y_c[i] - y_p[i]);
	// +delta from Qureshi and Prosser (2007)
	if( test > (qss_eps1_ * y_c[i] + delta_) ) {
	    ++flag;
	}
    }
    if( flag == 0 )
	return true;
    else
	return false;
}


double QssStep::step_suggest( double h, const vector<double> &y_c,
			      const vector<double> &y_p )
{
    double test = 0.0;
    double sigma = 0.0;
    double h_new = 0.0;

    for( int i = 0; i < ndim_; ++i ) {
	if( y_c[i] < ZERO_EPS )
	    continue;
	test = fabs( y_c[i] - y_p[i]) / (qss_eps2_ * y_c[i]);
	if( test > sigma )
	    sigma = test;
    }
    
    if( sigma <= 0.0 )
	h_new = h;
    else {
	// Estimate the sqrt of sigma using three Newton iterations
	// as Mott suggests.
	double x0 = sigma;
	double x1 = 0.0;
	for( int j = 0; j < 3; ++j) {
	    x1 = x0 - (x0*x0 - sigma)/(2*x0);
	    x0 = x1;
	}
	h_new = h * ( (1.0 / x1) + 0.005 );
    }

    return h_new;
}
