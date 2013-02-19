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

void sTPM::transform(const TPM &tpm){

   int a,b,c,d;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            for(int ay = 0;ay < Tools::gL();++ay)
               for(int by = 0;by < Tools::gL();++by)
                  for(int cy = 0;cy < Tools::gL();++cy)
                     for(int dy = 0;dy < Tools::gL();++dy)
                        (*this)(B,i,j) += tpm(S,Hamiltonian::gxy_a(a,ay),Hamiltonian::gxy_a(b,by),Hamiltonian::gxy_a(c,cy),Hamiltonian::gxy_a(d,dy));

            (*this)(B,i,j) /= (double) ( Tools::gL() * Tools::gL() );

         }
      }

   }

   this->symmetrize();

}

ostream &operator<<(ostream &output,const sTPM &stpm_p){

   int S,K;

   for(int B = 0;B < stpm_p.gnr();++B){

      S = stpm_p.block_char[B][0];
      K = stpm_p.block_char[B][1];

      output << "S =\t" << S << "\tK =\t" << K << "\tdimension =\t" << stpm_p.gdim(B) << "\tdegeneracy =\t" << stpm_p.gdeg(B) << std::endl;

      output << std::endl;

      for(int i = 0;i < stpm_p.gdim(B);++i)
         for(int j = 0;j < stpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << stpm_p.t2s[B][i][0] << "\t" << stpm_p.t2s[B][i][1]

               << "\t" << stpm_p.t2s[B][j][0] << "\t" << stpm_p.t2s[B][j][1] << "\t" << stpm_p(B,i,j) << endl;

         }

      output << std::endl;

   }

   return output;

}


/**
 * construct the spinsymmetrical hubbard hamiltonian in momentum space with on site repulsion U
 * @param U onsite repulsion term
 */
void sTPM::hubbard1D(double U){

   int k_a,k_b,k_c,k_d;//sp momentum 

   double ward = 1.0/(Tools::gN() - 1.0);

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         k_a = t2s[B][i][0];
         k_b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            k_c = t2s[B][j][0];
            k_d = t2s[B][j][1];

            //init
            (*this)(B,i,j) = 0;

            //hopping (kinetic energy):
            if(i == j)
               (*this)(B,i,i) = -2.0 * ward * ( cos(2.0*k_a*3.141592653589793238462 / (double) Tools::gL())  + cos(2.0*k_b*3.141592653589793238462 / (double) Tools::gL()) );

            //on-site repulsion
            if(S == 0){

               double ward = 2.0*U / (double) Tools::gL();

               if(k_a == k_b)
                  ward /= std::sqrt(2.0);

               if(k_c == k_d)
                  ward /= std::sqrt(2.0);

               (*this)(B,i,j) += ward;

            }

         }
      }

   }

   this->symmetrize();

}
