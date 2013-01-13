/** \file librad2.i
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 02-July-2009
 **/

%define RAD2_DOCSTRING
"Python interface to the radiation 2 library."
%enddef
%module(docstring=RAD2_DOCSTRING) librad2

%include "std_string.i"
%include "std_vector.i"

// The following magic allows us to pass a list of numbers
// to be collected as a C++ vector.
#ifndef STD_VECTOR_TEMPLATES_ALREADY_DEFINED
%template(vectori) std::vector<int>;
%template(vectord) std::vector<double>;
#define STD_VECTOR_TEMPLATES_ALREADY_DEFINED
#endif

// to enable use of cpointers
%include "cpointer.i"
/* Create some functions for working with "int *" and "double *" */
#ifndef POINTER_FUNCTIONS_ALREADY_DEFINED
%pointer_functions(int, intp);
%pointer_functions(double, doublep);
#define POINTER_FUNCTIONS_ALREADY_DEFINED
#endif

%{
#include "spectral_model.hh"
#include "atomic_line.hh"
#include "spectra_pieces.hh"
#include "radiation_constants.hh"
#include "LOS_pieces.hh"
#include "photaura.hh"
#include "radiator.hh"
#include "atomic_radiator.hh"
#include "diatomic_radiator.hh"
#include "diatomic_system.hh"
#include "cr_rr_coeffs.hh"
%}

%include "spectral_model.hh"
%include "atomic_line.hh"
%include "spectra_pieces.hh"
%include "radiation_constants.hh"
%include "LOS_pieces.hh"
%include "photaura.hh"
%include "radiator.hh"
%include "atomic_radiator.hh"
%include "diatomic_radiator.hh"
%include "diatomic_system.hh"
%include "cr_rr_coeffs.hh"

%template(vectorSB) std::vector<SpectralBin*>;

%pythoncode %{
class Bar:
    def __init__(self):
	self.len = 0
	self.chars = (' ', '+')
	self.wrap = ('[', ']')
	self.filledc = 0
	self.blank = 0
    
    def fill(self, i):
	assert not (i > 100) or (i < 0)
	self._setP(i)
    
    def _setP(self, p):
	self.filledc = int(round(float(self.len*p)/100))
    
    def show(self):
	out = []
	out.append(self.wrap[0])
	out.append(self.filledc*self.chars[1])
	out.append((self.len-self.filledc)*self.chars[0])
	out.append(self.wrap[1])
	return "".join(out)

# convert wavenumber in 1/cm to K
def nu2T( nu_cm ):
    return RC_h_SI * RC_c * nu_cm / RC_k_SI

# convert CEA heat of formation (J/mol) to cfcfd heat of formation (J/kg)
# note that the both h_f and m_w should be the values straight out thermo.inp
# in those units
def cea_Hf( h_f, m_w ):
    return h_f / ( m_w*1.0e-3 )

# CEA heat of formation (J/mol) to TAU heat of formation (K)
def cea_Hf2T( h_f ):
    return h_f / RC_Na / RC_k_SI 
    
# convert energy in Rydberg's to energy in Kelvin
def Ry2T( Ry ):
    return Ry * RC_Ry * RC_h_SI * RC_c_SI / RC_k_SI
    
# convert eV to nm
def eV2nm( eV ):
    return RC_h_SI * RC_c_SI / ( eV * RC_e_SI ) *1.0e9 

# convert eV to K
def eV2K( eV ):
    return eV * RC_e_SI / RC_k_SI 
    
# convert J to 1/cm
def J2nu( J ):
    return J / RC_h_SI / RC_c
    
# convert K to 1/cm
def T2nu( T ):
    return T * RC_k_SI / RC_h_SI / RC_c
%}
