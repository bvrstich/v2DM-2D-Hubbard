#ifndef sTPM_H
#define sTPM_H

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
 * This class sTPM is a class written for two particle matrices with spinsymmetry and translational symemtry included,
 * it inherits alle the function from its mother BlockMatrix, 
 * some special member functions and two lists that give the relationship between the sp and the tp basis.
 */
class sTPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the sTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const sTPM &tpm_p);

   public:
      
      //constructor
      sTPM();

      //copy constructor
      sTPM(const sTPM &);

      //destructor
      virtual ~sTPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      void transform(const TPM &);

      void hubbard1D(double);

      static void init();

      static void clear();

   private:

      //!static list that takes in a tp index i and a blockindex B, and returns two sp indices: a = t2s[B][i][0] and b = t2s[B][i][1]
      static vector< vector<int> > *t2s;

      //!static list that takes two sp indices a,b and a blockindex B, and returns a tp index i: i = s2t[B][a][b]
      static int ***s2t;

      //!static list that takes a blockindex B and returns the tp spin S and the tp momenta K_x and K_y
      static int **block_char;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
