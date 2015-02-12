/**
 * msg_service.d
 * Home to some commonly used messages in assert statements
 * and error messages.
 *
 * Author: Rowan J. Gollan
 */

module util.msg_service;

import std.string;

pure string brokenPreCondition(string variable, int lineNo, string fileName)
{
    return format("Pre-condition contract broken for %s on line %d in file %s\n",
		  variable, lineNo, fileName);
}

pure string brokenPostCondition(string variable, int lineNo, string fileName)
{
    return format("Post-condition contract broken for %s on line %d in file %s\n",
		  variable, lineNo, fileName);
}

string failedUnitTest(int lineNo, string fileName)
{
    return format("Unit test failure on line %d in file %s\n",
		  lineNo, fileName);
}
