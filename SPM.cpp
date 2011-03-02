#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes vector of dimension M/2, the SPM is completely diagonal in momentum space.
 * @param M dimension of single particle space
 * @param N nr of particles
 */
SPM::SPM(int M,int N) : Vector(M/2) {

   this->M = M;
   this->N = N;

}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : Vector(spm_copy) {

   this->M = spm_copy.gM();
   this->N = spm_copy.gN();

}

/**
 * TPM constructor: Creates a SPM initialized on the "bar" of the TPM.
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be initiated.
 */
SPM::SPM(double scale,const TPM &tpm) : Vector(tpm.gM()/2) {

   this->M = tpm.gM();
   this->N = tpm.gN();

   this->bar(scale,tpm);

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

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int k = 0;k < spm_p.gn();++k)
      output << k << "\t" << spm_p[k] << endl;

   return output;

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   for(int k = 0;k < M/2;++k){

      (*this)[k] = 0.0;

      for(int B = 0;B < tpm.gnr();++B){

         //p < k
         for(int p = 0;p < k;++p)
            (*this)[k] += tpm.gdeg(B) * tpm(B,k,p,k,p);

         //p == k: factor 2 for norm basis
         (*this)[k] += 2.0 * tpm.gdeg(B) * tpm(B,k,k,k,k);

         for(int p = k + 1;p < M/2;++p)
            (*this)[k] += tpm.gdeg(B) * tpm(B,k,p,k,p);

      }

      (*this)[k] *= 0.5 * scale;

   }

}
