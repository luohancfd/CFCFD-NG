#include "dbc_assert.hh"

using namespace std;

int main()
{
    try
    {
        ASSERT(true);
        ASSERT(false);
        ASSERT(false);
        ASSERT<std::runtime_error>(false, "this is the error message");
    }
    catch (std::exception& e)
    {
        cout << e.what();
    }
}
