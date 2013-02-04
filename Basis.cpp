#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * constructor 
 */
Basis::Basis(){

   basis = new TPM * [TPTPM::gn()];

   int B,I,J;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      int S = TPM::gblock_char(B,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      basis[i] = new TPM();

      *basis[i] = 0.0;

      if(I == J)
         (*basis[i])(B,I,I) = 1.0/std::sqrt(2.0*S + 1.0);
      else
         (*basis[i])(B,I,J) = (*basis[i])(B,J,I) = 1.0/std::sqrt(2.0*(2.0*S+1.0));

   }

}

/**
 * copy constructor 
 */
Basis::Basis(const Basis &basis_c){

   basis = new TPM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i)
      basis[i] = new TPM(basis_c[i]);

}

/**
 * Destructor
 */
Basis::~Basis(){

   for(int i = 0;i < TPTPM::gn();++i)
      delete basis[i];

   delete [] basis;

}

TPM &Basis::operator[](int i){

   return *basis[i];

}

const TPM &Basis::operator[](int i) const{

   return *basis[i];

}

ostream &operator<<(ostream &output,const Basis &basis_p){

   int B,I,J;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);
      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      output << endl;
      output << "basismatrix\t" << i << "\t|\t" << B << "\t" << I << "\t" << J << endl;
      output << endl;

   }

   return output;

}
