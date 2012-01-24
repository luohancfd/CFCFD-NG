// Author: Rowan J. Gollan
// Date: 26-May-2010
// Place: NASA Langley, Virginia, Hampton, USA
//
// History:
//   26-May-2010
//   Initial coding
//

#ifndef STRING_UTIL_HH
#define STRING_UTIL_HH

#include <string>
#include <vector>

int tokenize(const std::string &line, std::vector<std::string> &tokens);
int str2int(const std::string &s);
double str2dbl(const std::string &s);

#endif
