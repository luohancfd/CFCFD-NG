#include <iostream>

#include "no_fuss_linear_algebra.hh"

using namespace std;

int main()
{

    //    cout << "Initialise the NFVecs of length 100.\n";
    NFVec a(100, 10.0);
    NFVec b(100, 5.0);
    NFVec c(100, 15.0);

    //    cout << "Examine size in memory.. " << sizeof a << " bytes.\n";
    //    cout << "Data container size... " << a.get_data_size() << " bytes.\n";

    //    cout << "Now loop a million times and perform the following operations.\n";
    //    cout << "Add two NFVecs.\n";
    //  cout << "Subtract two NFVecs.\n";
    //    cout << "Scale an NFVec by a constant.\n";
    //    cout << "Scale an NFVec by a constant and store the result in another NFVec.\n";
    //    cout << "Copy one NFVec to another.\n";

    for(int i = 0; i < 1000000; ++i) {
	add_NFVecs(c, a, b);
	subtract_NFVecs(a, c, b);
	//scale_NFVec(a, 1.0);
	scale_NFVec2Vec(a, 3.0, c);
	copy_NFVec(a, c);
    }

    //    cout << "Done.\n";
    return 0;


}
