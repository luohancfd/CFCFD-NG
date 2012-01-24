// Author: Rowan J. Gollan
// Date: 22-Jul-2008

#ifndef PHYSICAL_CONSTANTS_HH
#define PHYSICAL_CONSTANTS_HH

// Universal gas constant
//#define PC_R_u 8.314472             // libgas2 value
#define PC_R_u 8.31451              // J/(mol.K) -- Tipler (1991)
#define PC_R_u_kmol (PC_R_u*1000.0) // J/(kmol.K) 
#define PC_R_cal 1.9872065    	    // cal/(mol.K)

// Avogadro's number
#define PC_Avogadro 6.02214e23      // molecules/g-mol -- B,S,L (2001)

// Boltzmann's constant
//#define PC_k_SI 1.3806505e-23     // libgas2 value
#define PC_k_SI 1.380658e-23        // J/K -- Tipler (1991)

// Boltzmann's constant (in C.G.S)
#define PC_k_CGS 1.3806505e-16      // erg/K

// Reference temperature (as per CEA program)
#define PC_T_ref 298.15             // K

// One atmosphere, in Pascals
#define PC_P_atm 101.325e3          // Pa

// Planck's constant (in S.I.)
#define PC_h_SI 6.626076e-34        // J.s

// Speed of light (CGS units)
#define PC_c 2.99792458e10          // cm/s

// Stephan-Boltzmann constant
#define PC_sigma_SI 5.6704e-8       // W/(m2.K4)

// Charge on an electron (in C.G.S)
#define PC_e_CGS 4.803e-10          // esu

// Charge on an electron (in S.I)
#define PC_e_SI 1.60218e-19         // Coulombs

// Mass of an electron (in C.G.S)
#define PC_m_CGS 9.109390e-28       // g

// Mass of an electron (in S.I)
#define PC_m_SI 9.109390e-31        // kg


#endif

// References:
//
// Tipler (1991)
// Physics for Scientists and Engineers, Volume 2
// Third Editition, Worth Publishers, New York
//
// Bird, Stewart and Lightfoot (2001)
// Transport Phenomena, 2nd edition
// John Wiley & Sons, New York
//
// Mohr, Taylor and Newell (2008)
// CODATA Recommended Values of the Fundamental Physical Constants: 2006
// Rev. Mod. Phys. 80: 633Ð730. doi:10.1103/RevModPhys.80.633.
