// bc_catalytic.hh
// FIX-ME Haven't put in default-constructors and copy-assignment operators. PJ 10-Mar-2013

const size_t MAX_EQ_WC_TABLE_ENTRIES = 100;

//-----------------------------------------------------------------
// A class to apply the catalytic-wall BC
class CatalyticWallBC {
public:
    CatalyticWallBC( int type_code );
    CatalyticWallBC( CatalyticWallBC &cw );
    virtual ~CatalyticWallBC();
    
    virtual int apply( Gas_data &Q, std::vector<double> &massf ) = 0;
    
public:
    int type_code;
};

class SuperCatalyticWallBC : public CatalyticWallBC {
public:
    SuperCatalyticWallBC( std::vector<double> massf_wall );
    SuperCatalyticWallBC( SuperCatalyticWallBC &cw );
    ~SuperCatalyticWallBC();
    
    int apply( Gas_data &Q, std::vector<double> &massf );
    
private:
    std::vector<double> massf_wall;
};

class EquilibriumCatalyticWallBC : public CatalyticWallBC {
public:
    EquilibriumCatalyticWallBC( std::string fname );
    EquilibriumCatalyticWallBC( EquilibriumCatalyticWallBC &cw );
    ~EquilibriumCatalyticWallBC();
    
    int apply( Gas_data &Q, std::vector<double> &massf );
    
private:
    std::vector<double> fC[MAX_EQ_WC_TABLE_ENTRIES];
    double lpmin, dlp;
    size_t ipmax;
};
