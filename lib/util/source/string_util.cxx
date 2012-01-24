// Author: Rowan J. Gollan
// Date: 26-May-2010
// Place: NASA Langley, Virginia, Hampton, USA
//
// History:
//   26-May-2010
//   Initial coding
//

#include <sstream>
#include <stdexcept>

#include "string_util.hh"

using namespace std;

int tokenize(const string &line, vector<string> &tokens)
{
    tokens.clear();
    stringstream parser(line);
    string temp;
    while ( parser >> temp ) {
	tokens.push_back(temp);
    }
    return tokens.size();

}

int str2int(const string &s)
{
    istringstream i(s);
    int x;
    if (!(i >> x))
	throw runtime_error("str2int(\"" + s + "\")");
    return x;
}

double str2dbl(const string &s)
{
    // From C++ FAQ Lite, 39.2
    istringstream i(s);
    double x;
    if (!(i >> x))
	throw runtime_error("str2dbl(\"" + s + "\")");
    return x;
}
