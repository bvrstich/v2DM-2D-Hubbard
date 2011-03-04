#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 11-05-2010\n\n
 * This class SPM was written for single particle matrices in a spinsymmetrical and translationally invariant system. It inherits from the class Vector
 * (because it is completely diagonal in k-space) and expands it with specific memberfunction and a knowledge of the nr of sp orbitals and particles.
 */

class SPM : public Vector {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << spm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << spm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM you want to print
    */
   friend ostream &operator<<(ostream &output,const SPM &spm_p);

   public:
      
      //constructor
      SPM();

      //copy constructor
      SPM(const SPM &);

      //TPM constructor
      SPM(double ,const TPM &);

      //PHM constructor
      SPM(double ,const PHM &);

      //destructor
      virtual ~SPM();

      using Vector::operator=;

      int gN() const;

      int gM() const;

      int gL() const;

      void bar(double,const TPM &);

      void bar(double,const PHM &);

      static void init(int,int);

   private:

      //!dimension of single particle space
      static int M;

      //!nr of particles
      static int N;

      //!dimension of the lattice
      static int L;

};

#endif
