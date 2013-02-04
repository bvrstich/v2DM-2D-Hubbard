#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <cstdlib>

using std::ostream;

#include "Vector.h"

class TPM;

/**
 * @author Brecht Verstichel
 * @date 18-12-2012\n\n
 * This class contains a basis for the TPM matrix space. For testing purposes constructe here to compare with hessian program.
 */
class Basis{

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param matrix_p de Basis you want to print
    */
   friend ostream &operator<<(ostream &output,const Basis &basis_p);

   public:

      //constructor
      Basis();

      //copy constructor
      Basis(const Basis &);

      //destructor
      virtual ~Basis();

      TPM &operator[](int);

      const TPM &operator[](int) const;

   private:

      //!array of TPM's, which is the basis
      TPM **basis;

};

#endif
