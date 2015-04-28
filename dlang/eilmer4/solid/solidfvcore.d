/**
 * solidfvcore.d
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-22-04
 */

module solidfvcore;

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
