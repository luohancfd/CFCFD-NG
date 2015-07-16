/** \file cns_diffusion.hh
 *  \ingroup mbcns2
 *  \brief Header file for the diffusion model class and functions.
 **/

#ifndef DIFFUSION_HH
#define DIFFUSION_HH

#include <string>

#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

using namespace std;

class DiffusionModel {
public:
    DiffusionModel(const string name, int nsp);
    
    DiffusionModel(const DiffusionModel &d);

    virtual ~DiffusionModel();
    
    virtual string str() const = 0;
    
    void fill_in_x(double rho, const vector<double> &massf);
    
    void fill_in_DAV_im(const matrix &D_AB);
    
    virtual void calculate_diffusion_fluxes(const Gas_data &Q,
					    double D_t,
					    const std::vector<double> &dfdx, 
					    const std::vector<double> &dfdy,
					    const std::vector<double> &dfdz,
					    std::vector<double> &jx, 
					    std::vector<double> &jy,
					    std::vector<double> &jz) = 0;
protected:
    string name_;
    int nsp_;
    int e_index_;
    double min_molef_;
    double small_DAB_;
    vector<double> DAV_im_;
    vector<double> x_;
    vector<double> M_;
    vector<double> Z_;
};

class NoDiffusionModel : public DiffusionModel {
public:
    NoDiffusionModel(const string name, int nsp);
    
    NoDiffusionModel(const NoDiffusionModel &d);

    virtual ~NoDiffusionModel();

    string str() const;

    void calculate_diffusion_fluxes(const Gas_data &Q,
				    double D_t,
				    const std::vector<double> &dfdx, 
				    const std::vector<double> &dfdy,
				    const std::vector<double> &dfdz,
				    std::vector<double> &jx, 
				    std::vector<double> &jy,
				    std::vector<double> &jz);
};

class FicksFirstLaw : public DiffusionModel {
public:
    FicksFirstLaw(const string name, int nsp);

    FicksFirstLaw(const FicksFirstLaw &f);

    virtual ~FicksFirstLaw();

    string str() const;

    void calculate_diffusion_fluxes(const Gas_data &Q,
				    double D_t,
				    const std::vector<double> &dfdx, 
				    const std::vector<double> &dfdy,
				    const std::vector<double> &dfdz,
				    std::vector<double> &jx, 
				    std::vector<double> &jy,
				    std::vector<double> &jz);
};

class StefanMaxwellModel : public DiffusionModel {
public:
    StefanMaxwellModel(const string name, int nsp, int no_iters);
    
    StefanMaxwellModel(const StefanMaxwellModel &s);

    virtual ~StefanMaxwellModel();

    string str() const;

    void calculate_diffusion_fluxes(const Gas_data &Q,
				    double D_t, 
				    const std::vector<double> &dfdx, 
				    const std::vector<double> &dfdy,
				    const std::vector<double> &dfdz,
				    std::vector<double> &jx, 
				    std::vector<double> &jy,
				    std::vector<double> &jz);

private:
    int no_iters_;
    std::vector<double> j0_;
    std::vector<double> j1_;

};

class RamshawChangModel : public DiffusionModel {
public:
    RamshawChangModel( const string name, int nsp );
    
    RamshawChangModel( const RamshawChangModel &rcm );

    virtual ~RamshawChangModel();

    string str() const;

    void calculate_diffusion_fluxes(const Gas_data &Q,
				    double D_t, 
				    const std::vector<double> &dfdx, 
				    const std::vector<double> &dfdy,
				    const std::vector<double> &dfdz,
				    std::vector<double> &jx, 
				    std::vector<double> &jy,
				    std::vector<double> &jz);
    
};

class ConstantLewisNumber : public DiffusionModel {
public:
    ConstantLewisNumber(const string name, int nsp, double Le);

    ConstantLewisNumber(const ConstantLewisNumber &c);

    virtual ~ConstantLewisNumber();

    string str() const;

    void calculate_diffusion_fluxes(const Gas_data &Q,
				    double D_t,
				    const std::vector<double> &dfdx, 
				    const std::vector<double> &dfdy,
				    const std::vector<double> &dfdz,
				    std::vector<double> &jx, 
				    std::vector<double> &jy,
				    std::vector<double> &jz);
protected:
    double Le_;
};

class ConstantSchmidtNumber : public DiffusionModel {
public:
    ConstantSchmidtNumber(const string name, int nsp, double Sc);

    ConstantSchmidtNumber(const ConstantSchmidtNumber &c);

    virtual ~ConstantSchmidtNumber();

    string str() const;

    void calculate_diffusion_fluxes(const Gas_data &Q,
                                    double D_t,
                                    const std::vector<double> &dfdx,
                                    const std::vector<double> &dfdy,
                                    const std::vector<double> &dfdz,
                                    std::vector<double> &jx,
                                    std::vector<double> &jy,
                                    std::vector<double> &jz);
protected:
    double Sc_;
};

class ConstantLewisNumber_DRM19 : public DiffusionModel {
public:
    ConstantLewisNumber_DRM19(const string name, int nsp );

    ConstantLewisNumber_DRM19(const ConstantLewisNumber_DRM19 &c);

    virtual ~ConstantLewisNumber_DRM19();

    string str() const;

    void calculate_diffusion_fluxes(const Gas_data &Q,
				    double D_t,
				    const std::vector<double> &dfdx, 
				    const std::vector<double> &dfdy,
				    const std::vector<double> &dfdz,
				    std::vector<double> &jx, 
				    std::vector<double> &jy,
				    std::vector<double> &jz);
};

int set_diffusion_model( const string diffusion_model="Stefan-Maxwell" );

void calculate_diffusion_fluxes(const Gas_data &Q,
				double D_t,
				const std::vector<double> &dfdx, 
				const std::vector<double> &dfdy,
				const std::vector<double> &dfdz,
				std::vector<double> &jx, 
				std::vector<double> &jy,
				std::vector<double> &jz);

#endif
