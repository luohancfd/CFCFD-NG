/** \file zero_system.cxx
 * \brief  A class for creating a system for a zero finding technique.
 *
 * This collectes the required functions fo the abstract ZeroSystem class.
 *
 * \author Rowan J Gollan
 * \date 24-Apr-2006
 *
 **/

#include "no_fuss_linear_algebra.hh"
#include "zero_system.hh"

using namespace std;

ZeroSystem::ZeroSystem() {}
ZeroSystem::ZeroSystem( const ZeroSystem &z ) {}
ZeroSystem::~ZeroSystem() {}

ZeroFunction::ZeroFunction() {}
ZeroFunction::ZeroFunction( const ZeroFunction &z ) {}
ZeroFunction::~ZeroFunction() {}
