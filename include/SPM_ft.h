#ifndef SPM_ft_H
#define SPM_ft_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

class PHM;
class PPHM;

/**
 * @author Brecht Verstichel
 * @date 17-09-2012\n\n
 * This class SPM_ft was written for single particle matrices in a spinsymmetrical and translationally invariant system. It is the fourier transform
 * of the momentum SPM class.
 */

class SPM_ft : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM_ft you want to print
    */
   friend ostream &operator<<(ostream &output,const SPM_ft &spm_p);

   public:
      
      //constructor
      SPM_ft();

      //copy constructor
      SPM_ft(const SPM_ft &);

      SPM_ft(const SPM &);

      //destructor
      virtual ~SPM_ft();

      using Matrix::operator=;

      using Matrix::operator();

   private:

};

#endif
