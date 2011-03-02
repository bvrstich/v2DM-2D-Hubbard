#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::ifstream;
using std::vector;

#include "BlockMatrix.h"

class SUP;

/**
 * @author Brecht Verstichel
 * @date 10-05-2010\n\n
 * This class TPM is a class written for two particle matrices with spinsymmetry and translational symemtry included, it inherits alle the function from its mother 
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the tp basis.
 */
class TPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPM &tpm_p);

   public:
      
      //constructor
      TPM();

      //copy constructor
      TPM(const TPM &);

      //destructor
      virtual ~TPM();

      void constr_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode and blockindex
      double operator()(int B,int a,int b,int c,int d) const;

      //easy to access the numbers, in sp mode and with tp spin and momentum quantumnumber
      double operator()(int S,int K,int a,int b,int c,int d) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //geef L terug
      int gL() const;

      void hubbard(double U);

      //Q afbeelding en zijn inverse
      void Q(int option,const TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,const TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,const TPM &);

      void unit();

      void proj_Tr();

      //de hessiaan afbeelding:
      void H(const TPM &b,const SUP &D);

      //los het stelsel op
      int solve(TPM &b,const SUP &D);

      void min_unit(double scale);

      void collaps(int option,const SUP &);

      //return the spin
      double spin() const;

      //output to file
      void out_sp(const char *) const;

      //input from file
      void in(ifstream &);

      double trace_pair() const;

      void set_S_2();

      static void init(int,int);

      static void clear();

   private:

      static vector< vector<int> > *t2s;

      static int ***s2t;

      //!static list that takes a blockindex B and returns the tp spin S and the tp momenta K_x and K_y: S = block_char[B][0] , K_x = block_char[B][1], K_y = block_char[B][2]
      static int **block_char;

      //!list of 6j symbols needed.
      static double **_6j;

      //!static counter that counts the number of TPM objects running in the program
      static int counter;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!dimension of the lattice
      static int L;

};

#endif
