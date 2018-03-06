#include <iostream>
#include <string>

using namespace std;

/*
The goal of this dummy file is to force libtool to use g++ when a 
C++ library like Angpow is linked to CCL
*/
void angpow_cpp_placeholder()
{
  cout << " --- this is a dummy function to ensure that libtool uses the c++ driver to link / create libraries when c++ libraries are used with CCL"<<endl;
  return;
}
