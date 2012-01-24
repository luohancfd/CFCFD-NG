// Author: Rowan J. Gollan
// Date: 17-Oct-2008
// Place: Hampton, Virginia, USA

#include <string>

#include "../../util/source/lua_service.hh"
#include "ode_setup.hh"

using namespace std;

OdeSolver* create_ode_solver(lua_State *L, int nsp, string name)
{
    string step_routine(get_string(L, -1, "step_routine"));
    int max_step_attempts = get_positive_int(L, -1, "max_step_attempts");
    double max_increase_factor = get_positive_number(L, -1, "max_increase_factor");
    double max_decrease_factor = get_positive_number(L, -1, "max_decrease_factor");
    double decrease_factor = get_positive_number(L, -1, "decrease_factor");

    return new OdeSolver(name,
			 nsp, step_routine, max_step_attempts,
			 max_increase_factor, max_decrease_factor,
			 decrease_factor);
}
