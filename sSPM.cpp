#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::cout;
using std::endl;

#include "include.h"

/**
 * constructor, makes vector of dimension L*L, the sSPM is completely diagonal in momentum space, both in k_x and k_y.
 */
sSPM::sSPM(){

   sspm = new double [Tools::gL()];

}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
sSPM::sSPM(const sSPM &sspm_copy)  {

   sspm = new double [Tools::gL()];

   for(int i = 0;i < Tools::gL();++i)
      sspm[i] = sspm_copy[i];

}

/**
 * destructor
 */
sSPM::~sSPM(){

   delete [] sspm;

}

double &sSPM::operator[](int i){

   return sspm[i];

}

double sSPM::operator[](int i) const{

   return sspm[i];

}

ostream &operator<<(ostream &output,const sSPM &sspm_p){

   for(int a = 0;a < Tools::gL();++a)
      output << a << "\t" << sspm_p[a] << endl;

   return output;

}

void sSPM::transform(const SPM &spm){

   for(int x = 0;x < Tools::gL();++x){

      sspm[x] = 0;

      for(int y = 0;y < Tools::gL();++y)
         sspm[x] += spm[Hamiltonian::gxy_a(x,y)];

      sspm[x] /= (double)Tools::gL();

   }

}
