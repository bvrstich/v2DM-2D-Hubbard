#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::cout;
using std::endl;

#include "include.h"

/**
 * constructor, makes vector of dimension Tools::gL()*Tools::gL(), the SPM is completely diagonal in momentum space, both in k_x and k_y.
 */
SPM::SPM() : Vector(Tools::gL()*Tools::gL()) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : Vector(spm_copy) { }

/**
 * destructor
 */
SPM::~SPM(){ }

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

   for(int a = 0;a < Tools::gL()*Tools::gL();++a){

      (*this)[a] = 0.0;

      //S = 0
      //b < a
      for(int b = 0;b < a;++b)
         (*this)[a] += tpm(0,a,b,a,b);

      //b == a
      (*this)[a] += 2.0*tpm(0,a,a,a,a);

      //b > a
      for(int b = a + 1;b < Tools::gL()*Tools::gL();++b)
         (*this)[a] += tpm(0,a,b,a,b);

      //S = 1
      for(int b = 0;b < Tools::gL()*Tools::gL();++b)
         (*this)[a] += 3.0*tpm(1,a,b,a,b);

      (*this)[a] *= 0.5 * scale;

   }

}

/**
 * Trace out a set of indices to create the "bar" matrix of a PHM, slight difference from the bar(TPM) function (normalization of the tp basisset).
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be filled
 */
void SPM::bar(double scale,const PHM &phm){

   for(int a = 0;a < Tools::gL()*Tools::gL();++a){

      (*this)[a] = 0.0;

      for(int b = 0;b < Tools::gL()*Tools::gL();++b)
         (*this)[a] += phm(0,a,b,a,b) + 3.0*phm(1,a,b,a,b);

      (*this)[a] *= 0.5 * scale;

   }

}

/** 
 * This bar function maps a PPHM object directly onto a SPM object, scaling it with a factor scale
 * @param scale the scalefactor
 * @param pphm Input PPHM object
 */
void SPM::bar(double scale,const PPHM &pphm){

   for(int c = 0;c < Tools::gL()*Tools::gL();++c){

      (*this)[c] = 0.0;

      //first S = 1/2 part
      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < Tools::gL()*Tools::gL();++a){

            for(int b = 0;b < a;++b)//b < a
               (*this)[c] += pphm.pph(0,S_ab,a,b,c,S_ab,a,b,c);

            (*this)[c] += 2.0 * pphm.pph(0,S_ab,a,a,c,S_ab,a,a,c);

            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b)//b > a
               (*this)[c] += pphm.pph(0,S_ab,a,b,c,S_ab,a,b,c);

         }

      }

      //then S = 3/2 part:
      for(int a = 0;a < Tools::gL()*Tools::gL();++a)
         for(int b = 0;b < Tools::gL()*Tools::gL();++b)
            (*this)[c] += 2.0 * pphm.pph(1,1,a,b,c,1,a,b,c);

      //scaling
      (*this)[c] *= scale;

   }

}
