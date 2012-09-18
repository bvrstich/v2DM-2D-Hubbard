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
class PHM;
class DPM;
class PPHM;

/**
 * @author Brecht Verstichel
 * @date 10-05-2010\n\n
 * This class TPM is a class written for two particle matrices with spinsymmetry and translational symemtry included,
 * It inherits alle the function from its mother BlockMatrix.
 * Some special member functions and two lists that give the relationship between the sp and the tp basis.
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

      //easy to access the numbers, in sp mode and with tp spin and momentum quantumnumbers
      double operator()(int S,int a,int b,int c,int d) const;

      void hubbard(double U);

      //Q afbeelding en zijn inverse
      void Q(int option,const TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double,double,double,const TPM &);

      void unit();

      void proj_Tr();

      void min_unit(double scale);

      void collaps(int option,const SUP &);

      //return the spin
      double spin() const;

      //output to file
      void out_sp(const char *) const;

      void set_S_2();

      void G(const PHM &);

      void bar(const DPM &);

      void T(const DPM &);

      void T(const PPHM &);

      void bar(const PPHM &);

      int solve(TPM &,const SUP &);

      void H(const TPM &,const SUP &);

      void S(int,const TPM &);

      static void init_overlap();

      static void init();

      static void clear();

   private:

      //!static list that takes in a tp index i and a blockindex B, and returns two sp indices: a = t2s[B][i][0] and b = t2s[B][i][1]
      static vector< vector<int> > *t2s;

      //!static list that takes two sp indices a,b and a blockindex B, and returns a tp index i: i = s2t[B][a][b]
      static int ***s2t;

      //!static list that takes a blockindex B and returns the tp spin S and the tp momenta K_x and K_y
      static int **block_char;

      //!static list that returns the blockindex when given the S, K_x and K_y.
      static int ***char_block;

      //!overlapmatrix parameters
      static double Sa,Sc;

};

#endif
