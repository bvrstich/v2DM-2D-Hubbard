#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::cout;
using std::endl;

#include "include.h"

int SPM::L;
int SPM::N;
int SPM::M;

/**
 * static function that initializes the static variables
 */
void SPM::init(int L_in,int N_in){

   L = L_in;
   N = N_in;

   M = L*L*2;

}

/**
 * constructor, makes vector of dimension L*L, the SPM is completely diagonal in momentum space, both in k_x and k_y.
 */
SPM::SPM() : Vector(L*L) {

}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : Vector(spm_copy) {

}

/**
 * TPM constructor: Creates a SPM initialized on the "bar" of the TPM.
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be initiated.
 */
SPM::SPM(double scale,const TPM &tpm) : Vector(L*L) {

   this->bar(scale,tpm);

}

/**
 * PHM constructor: Creates a SPM initialized on the "bar" of the PHM.
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be initiated.
 */
SPM::SPM(double scale,const PHM &phm) : Vector(L*L) {

   this->bar(scale,phm);

}

/**
 * destructor
 */
SPM::~SPM(){

}

/**
 * @return nr of particles
 */
int SPM::gN() const{

   return N;

}

/**
 * @return dimension of sp space
 */
int SPM::gM() const{

   return M;

}

/**
 * @return dimension of the lattice
 */
int SPM::gL() const{

   return L;

}

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int a = 0;a < spm_p.gn();++a)
      output << a << "\t" << spm_p[a] << endl;

   return output;

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   for(int a = 0;a < L*L;++a){

      (*this)[a] = 0.0;

      //S = 0

      //b < a
      for(int b = 0;b < a;++b)
         (*this)[a] += tpm(0,(Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L,(Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L,a,b,a,b);

      //b == a
      (*this)[a] += 2.0*tpm(0,(Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(a,0))%L,(Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(a,1))%L,a,a,a,a);

      //b > a
      for(int b = a + 1;b < L*L;++b)
         (*this)[a] += tpm(0,(Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L,(Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L,a,b,a,b);

      //S = 1

      for(int b = 0;b < L*L;++b)
         (*this)[a] += 3.0*tpm(1,(Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L,(Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L,a,b,a,b);

      (*this)[a] *= 0.5 * scale;

   }

}

/**
 * Trace out a set of indices to create the "bar" matrix of a PHM, slight difference from the bar(TPM) function (normalization of the tp basisset).
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be filled
 */
void SPM::bar(double scale,const PHM &phm){

   int K_x,K_y;

   for(int a = 0;a < L*L;++a){

      (*this)[a] = 0.0;

      for(int b = 0;b < L*L;++b){

         K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L;
         K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L;

         (*this)[a] += phm(0,K_x,K_y,a,b,a,b) + 3.0*phm(1,K_x,K_y,a,b,a,b);

      }

      (*this)[a] *= 0.5 * scale;

   }

}
