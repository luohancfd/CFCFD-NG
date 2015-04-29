/**
 * solidfvcore.d
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-22-04
 */

module solidfvcore;

import std.string;
import std.format;

enum SolidDomainUpdate {
    euler,
    pc
}

string solidDomainUpdateSchemeName(SolidDomainUpdate sdupdate)
{
    final switch ( sdupdate ) {
    case SolidDomainUpdate.euler: return "euler";
    case SolidDomainUpdate.pc: return "predictor-corrector";
    }
}

int numberOfStagesForUpdateScheme(SolidDomainUpdate sdupdate)
{
    final switch ( sdupdate ) {
    case SolidDomainUpdate.euler: return 1;
    case SolidDomainUpdate.pc: return 2;
    }
}

SolidDomainUpdate solidDomainUpdateSchemeFromName(string name)
{
    switch (name) {
    case "euler": return SolidDomainUpdate.euler;
    case "pc": return SolidDomainUpdate.pc;
    case "predictor-corrector": return SolidDomainUpdate.pc;
    case "predictor_corrector": return SolidDomainUpdate.pc;
    default:
	string errMsg = format("Invalid solid domain updated scheme name: %s", name);
	throw new Error(errMsg);
    }
}
