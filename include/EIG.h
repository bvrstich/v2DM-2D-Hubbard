#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockVector.h"
#include "SUP.h"

#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __Q_CON
#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

/**
 * @author Brecht Verstichel
 * @date 06-05-2010\n\n
 * This class, EIG is a "block"-vector over the carrierspace's of the active condtions. It contains room
 * to store the eigenvalues and special member function that work with these eigenvalues.
 * This class should only be used when a SUP matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 */
class EIG{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << eig_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << eig_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,const EIG &eig_p);

   public:

      //constructor met initialisatie op 
      EIG(SUP &);

      //copy constructor
      EIG(const EIG &);

      //destructor
      ~EIG();

      //overload equality operator
      EIG &operator=(const EIG &);

      BlockVector<TPM> &gvI();

      const BlockVector<TPM> &gvI() const;

#ifdef __Q_CON
      BlockVector<TPM> &gvQ();

      const BlockVector<TPM> &gvQ() const;
#endif

#ifdef __G_CON
      BlockVector<PHM> &gvG();

      const BlockVector<PHM> &gvG() const;
#endif

#ifdef __T1_CON
      BlockVector<DPM> &gvT1();

      const BlockVector<DPM> &gvT1() const;
#endif

#ifdef __T2_CON
      BlockVector<PPHM> &gvT2();

      const BlockVector<PPHM> &gvT2() const;
#endif

      double min() const;

      double max() const;

      double lsfunc(double) const;

   private:

      //!pointer to a BlockVector<TPM> object, the eigenvalues of the P part of a SUP matrix are stored here
      BlockVector<TPM> *vI;

#ifdef __Q_CON
      BlockVector<TPM> *vQ;
#endif

#ifdef __G_CON
      //!single pointer to a BlockVector<PHM> object, the eigenvalues of G part of a SUP matrix will be stored here.
      BlockVector<PHM> *vG;
#endif

#ifdef __T1_CON
      //!single pointer to a BlockVector<DPM> object, the eigenvalues of T1 part of a SUP matrix will be stored here.
      BlockVector<DPM> *vT1;
#endif

#ifdef __T2_CON
      //!single pointer to a BlockVector<PPHM> object, the eigenvalues of T2 part of a SUP matrix will be stored here.
      BlockVector<PPHM> *vT2;
#endif

};

#endif
