#ifndef PHM_H
#define PHM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"

class PPHM;

/**
 * @author Brecht Verstichel
 * @date 11-05-2010\n\n
 * This class, PHM, is a class written for spinsymmetrical, translationally invariant particle-hole matrices
 * it inherits all the functions from its mother class
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the ph basis.
 */
class PHM : public BlockMatrix {

    /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << phm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << phm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param phm_p the PHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PHM &phm_p);

   public:
      
      //constructor
      PHM();

      //copy constructor
      PHM(const PHM &);

      //destructor
      virtual ~PHM();

      void constr_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int S,int a,int b,int c,int d) const;

      void G(const TPM &);

      void G(int,const TPM &);

      void out_sp(const char *) const;

      void bar(const PPHM &);

      static void init();

      static void clear();

   private:

      //!static list that takes in a blockindex B, a ph index i and returns two sp indices: a = ph2s[B][i][0] and b = ph2s[B][i][1]
      static vector< vector<int> > *ph2s;

      //!static list that takes in a blockindex B and two sp indices a,b and returns a ph index i: i = s2ph[B][a][b]
      static int ***s2ph;

      //!static list that takes a blockindex B and returns the ph spin S and the ph momenta K_x and K_y.
      static int **block_char;

      //!static list that returns the blockindex when given the S, K_x and K_y.
      static int ***char_block;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
