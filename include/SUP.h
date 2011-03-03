#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;
using std::ofstream;
using std::ifstream;

#include "TPM.h"

class EIG;

/**
 * @author Brecht Verstichel
 * @date 09-03-2010\n\n
 * This class, SUP is a blockmatrix over the carrierspace's of active N-representability conditions. 
 * This class contains two TPM objects, and if compiled with the right option a PHM or DPM object, 
 * You have to remember that these matrices are independent of each other (by which I mean that TPM::Q(SUP_PQ::tpm (0))
 * is not neccesarily equal to SUP_PQ::tpm (1)) etc. .
 */
class SUP{
  
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param SZ_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,const SUP &SZ_p);

   public:

      //constructor
      SUP();

      //copy constructor
      SUP(const SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(const SUP &);

      //overload -= operator
      SUP &operator-=(const SUP &);

      //overload equality operator
      SUP &operator=(const SUP &);

      //overload equality operator
      SUP &operator=(double &);

      TPM &tpm(int i);

      const TPM &tpm(int i) const;

      //initialiseer S
      void init_S();

      //initialiseer Z
      void init_Z(double alpha,const TPM &ham,const SUP &u_0);

      int gN() const;

      int gM() const;

      int gL() const;

      int gn_tp() const;

      int gdim() const;

      double ddot(const SUP &) const;

      void invert();

      void dscal(double alpha);

      void proj_U();

      void proj_C(const TPM &);

      //maak de matrix D, nodig voor de hessiaan van het stelsel
      void D(const SUP &S,const SUP &Z);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(const SUP &,const SUP &);

      void daxpy(double alpha,const SUP &);

      void proj_C();

      SUP &mprod(const SUP &,const SUP &);

      void fill(const TPM &);

      void fill();

      int solve(SUP &B,const SUP &D);

      void H(const SUP &B,const SUP &D);

      double center_dev(const SUP &Z) const;

      double line_search(const SUP &DZ,const SUP &S,const SUP &Z,double max_dev) const;

      void fill_Random();

      void out(ofstream &) const;

      void in(ifstream &);

      static void init(int,int);

   private:

      //!double pointer of TPM's, will contain the P and Q block of the SUP in the first and second block.
      TPM **SZ_tp;

      //!number of sp orbitals
      static int M;

      //!nr of particles
      static int N;

      //!dimension of the lattice
      static int L;

      //!dimension of tp space
      static int n_tp;

      //!total dimension of the SUP matrix
      static int dim;

};

#endif
