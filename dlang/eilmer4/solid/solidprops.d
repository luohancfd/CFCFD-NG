/**
 * solidprops.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04
 */

module solidprops;

class SolidProps {
public:
    double rho;
    double k;
    double Cp;
    double k11;
    double k12;
    double k22;

    this(double rho_, double k_, double Cp_,
	 double k11_=0.0, double k12_=0.0, double k22_=0.0)
    {
	rho = rho_;
	k = k_;
	Cp = Cp_;
	k11 = k11_;
	k12 = k12_;
	k22 = k22_;
    }
}

double updateEnergy(SolidProps sp, double T)
{
    return sp.rho*sp.Cp*T;
}

double updateTemperature(SolidProps sp, double e)
{
    return e/(sp.rho*sp.Cp);
}
