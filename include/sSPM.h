#ifndef sSPM_H
#define sSPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"
#include "PPHM.h"

/**
 * @author Brecht Verstichel
 * @date 11-05-2010\n\n
 * This class sSPM was written for single particle matrices in a spinsymmetrical and translationally invariant system. It inherits from the class Vector
 * (because it is completely diagonal in k-space) and expands it with specific memberfunction and a knowledge of the nr of sp orbitals and particles.
 */

class sSPM {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << spm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << spm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de sSPM you want to print
    */
   friend ostream &operator<<(ostream &output,const sSPM &spm_p);

   public:
      
      //constructor
      sSPM();

      //copy constructor
      sSPM(const sSPM &);

      //destructor
      virtual ~sSPM();

      double &operator[](int);

      double operator[](int) const;

      void transform(const SPM &);

   private:

      double *sspm;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
