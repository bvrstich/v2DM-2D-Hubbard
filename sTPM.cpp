#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

vector< vector<int> > *sTPM::t2s;
int ***sTPM::s2t;

int **sTPM::block_char;

/**
 * static function that initializes the static variables
 */
void sTPM::init(){

   //allocate stuff
   t2s = new vector< vector<int> > [2 * Tools::gL()];

   s2t = new int ** [2 * Tools::gL()];

   for(int B = 0;B < 2 * Tools::gL();++B){

      s2t[B] = new int * [Tools::gL()];

      for(int a = 0;a < Tools::gL();++a)
         s2t[B][a] = new int [Tools::gL()];

   }

   block_char = new int * [2*Tools::gL()];

   for(int B = 0;B < 2 * Tools::gL();++B)
      block_char[B] = new int [2];

   vector<int> v(2);

   int block = 0;

   //tp index
   int t;

   for(int K = 0;K < Tools::gL();++K){

      t = 0;

      //S = 0
      block_char[block][0] = 0;
      block_char[block][1] = K;

      for(int a = 0;a < Tools::gL();++a)
         for(int b = a;b < Tools::gL();++b){

            if( (a + b)%Tools::gL() == K ){

               v[0] = a;
               v[1] = b;

               t2s[block].push_back(v);

               s2t[block][a][b] = t;
               s2t[block][b][a] = t;

               ++t;

            }

         }

      t = 0;

      //S = 1
      block_char[Tools::gL() + block][0] = 1;
      block_char[Tools::gL() + block][1] = K;

      for(int a = 0;a < Tools::gL();++a)
         for(int b = a + 1;b < Tools::gL();++b){

            if( (a + b)%Tools::gL() == K ) {

               v[0] = a;
               v[1] = b;

               t2s[Tools::gL() + block].push_back(v);

               s2t[Tools::gL() + block][a][b] = t;
               s2t[Tools::gL() + block][b][a] = t;

               ++t;

            }

         }

      ++block;

   }

}

/**
 * static function that deallocates the static variables
 */
void sTPM::clear(){

   delete [] t2s;

   for(int B = 0;B < 2*Tools::gL();++B){

      for(int a = 0;a < Tools::gL();++a)
         delete [] s2t[B][a];

      delete [] s2t[B];

      delete [] block_char[B];

   }

   delete [] s2t;

   delete [] block_char;

}

/**
 * standard constructor for a spinsymmetrical, translationally invariant tp matrix on a 1D lattice: 
 * constructs BlockMatrix object with 2 * L blocks
 */
sTPM::sTPM() : BlockMatrix(2*Tools::gL()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL();++B)//S = 0
      setMatrixDim(B,t2s[B].size(),1);

   for(int B = Tools::gL();B < 2*Tools::gL();++B)//S = 1
      setMatrixDim(B,t2s[B].size(),3);

}

/**
 * copy constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice:
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 * @param tpm_c The sTPM object to be copied into (*this)
 */
sTPM::sTPM(const sTPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor
 */
sTPM::~sTPM(){ }
