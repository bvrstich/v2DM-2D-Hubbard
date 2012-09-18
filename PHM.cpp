#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::endl;

#include "include.h"

vector< vector<int> > *PHM::ph2s;
int ***PHM::s2ph;

int **PHM::block_char;
int ***PHM::char_block;

/**
 * initializes the statics
 */
void PHM::init(){

   //allocate stuff
   ph2s = new vector< vector<int> > [Tools::gM()];

   s2ph = new int ** [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B){

      s2ph[B] = new int * [Tools::gL()*Tools::gL()];

      for(int a = 0;a < Tools::gL()*Tools::gL();++a)
         s2ph[B][a] = new int [Tools::gL()*Tools::gL()];

   }

   block_char = new int * [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [Tools::gL()];

      for(int x = 0;x < Tools::gL();++x)
         char_block[S][x] = new int [Tools::gL()];

   }

   vector<int> v(2);

   int block = 0;

   //ph index
   int ph;

   //loop over the K_x K_y blocks
   for(int K_x = 0;K_x < Tools::gL();++K_x)
      for(int K_y = 0;K_y < Tools::gL();++K_y){

         ph = 0;

         //S = 0
         block_char[block][0] = 0;
         block_char[block][1] = K_x;
         block_char[block][2] = K_y;

         char_block[0][K_x][K_y] = block;

         //S = 1
         block_char[Tools::gL()*Tools::gL() + block][0] = 1;
         block_char[Tools::gL()*Tools::gL() + block][1] = K_x;
         block_char[Tools::gL()*Tools::gL() + block][2] = K_y;

         char_block[1][K_x][K_y] = Tools::gL()*Tools::gL() + block;

         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = 0;b < Tools::gL()*Tools::gL();++b){

               if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % Tools::gL() == K_x )
                  
                     && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % Tools::gL() == K_y ) ){

                  v[0] = a;
                  v[1] = b;

                  ph2s[block].push_back(v);//S = 0
                  ph2s[Tools::gL()*Tools::gL() + block].push_back(v);//S = 1

                  s2ph[block][a][b] = ph;
                  s2ph[Tools::gL()*Tools::gL() + block][a][b] = ph;

                  ++ph;

               }

            }

         ++block;

      }

}

/**
 * static function that deallocates the static variables
 */
void PHM::clear(){

   delete [] ph2s;

   for(int B = 0;B < Tools::gM();++B){

      for(int a = 0;a < Tools::gL()*Tools::gL();++a)
         delete [] s2ph[B][a];

      delete [] s2ph[B];

      delete [] block_char[B];

   }

   delete [] s2ph;

   delete [] block_char;

   for(int S = 0;S < 2;++S){

      for(int x = 0;x < Tools::gL();++x)
         delete [] char_block[S][x];

      delete [] char_block[S];

   }

   delete [] char_block;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 0 and 1.
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 */
PHM::PHM() : BlockMatrix(Tools::gM()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL()*Tools::gL();++B)//S = 0
      setMatrixDim(B,ph2s[B].size(),1);

   for(int B = Tools::gL()*Tools::gL();B < Tools::gM();++B)//S = 1
      setMatrixDim(B,ph2s[B].size(),3); 

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks of dimension M*M/4 and copies the content of phm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : BlockMatrix(phm_c){ }

/**
 * destructor: if counter == 1 the memory for the static lists ph2s en s2ph will be deleted.
 */
PHM::~PHM(){ }

ostream &operator<<(ostream &output,const PHM &phm_p){

   int S,K_x,K_y;

   for(int B = 0;B < phm_p.gnr();++B){

      S = phm_p.block_char[B][0];
      K_x = phm_p.block_char[B][1];
      K_y = phm_p.block_char[B][2];

      output << "S =\t" << S << "\tK_x =\t" << K_x << "\tK_y =\t" << K_y << "\tdimension =\t" << phm_p.gdim(B) << "\tdegeneracy =\t" << phm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < phm_p.gdim(B);++i)
         for(int j = 0;j < phm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << phm_p.ph2s[B][i][0] << "\t" << phm_p.ph2s[B][i][1]

               << "\t" << phm_p.ph2s[B][j][0] << "\t" << phm_p.ph2s[B][j][1] << "\t" << phm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * access the elements of the matrix in sp mode, 
 * @param S The tp spin quantumnumber
 * @param a first sp index that forms the ph row index i in block B together with b
 * @param b second sp index that forms the ph row index i in block B together with a
 * @param c first sp index that forms the ph column index j in block B together with d
 * @param d second sp index that forms the ph column index j in block B together with c
 * @return the number on place PHM(i,j)
 */
double PHM::operator()(int S,int a,int b,int c,int d) const{

   int K_x =  ( Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) )%Tools::gL();
   int K_y =  ( Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) )%Tools::gL();

   //check if momentum is conserved
   if( K_x != ( Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(d,0) )%Tools::gL() )
      return 0;

   if( K_y != ( Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(d,1) )%Tools::gL() ) 
      return 0;

   int B = char_block[S][K_x][K_y];

   int i = s2ph[B][a][b];
   int j = s2ph[B][c][d];

   return (*this)(B,i,j);

}

/**
 * The G map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G(const TPM &tpm){

   //construct the SPM corresponding to the TPM
   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   int a,b,c,d;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         a = ph2s[B][i][0];

         //transform k_b to tpm sp-momentum:
         b = Hamiltonian::bar(ph2s[B][i][1]);

         for(int j = i;j < gdim(B);++j){

            c = ph2s[B][j][0];

            //transform k_d to tpm sp-momentum:
            d = Hamiltonian::bar(ph2s[B][j][1]);

            (*this)(B,i,j) = - Tools::g6j(0,0,0,S) * tpm(0,a,d,c,b) - 3.0 * Tools::g6j(0,0,1,S) * tpm(1,a,d,c,b);

            if(a == d)
               (*this)(B,i,j) *= std::sqrt(2.0);

            if(b == c)
               (*this)(B,i,j) *= std::sqrt(2.0);

         }

         (*this)(B,i,i) += spm[a];

      }

   }

   this->symmetrize();

}

/**
 * output the PHM object for test purposes in sp-mode to be used in spin_pd program.
 * @param filename the filename
 */
void PHM::out_sp(const char *filename) const{

   ofstream output(filename);
   output.precision(15);

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i)
         for(int j = 0;j < gdim(B);++j)
            output << block_char[B][0] << "\t" << ph2s[B][i][0] << "\t" << ph2s[B][i][1] << "\t" << ph2s[B][j][0] << "\t" << ph2s[B][j][1] 

               << "\t" << (*this)(B,i,j) << endl;

   }

}


/**
 * The bar function that maps a PPHM object onto a PHM object by tracing away the first pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void PHM::bar(const PPHM &pphm){

   int a,b,c,d;

   double ward,hard;

   int S;

   for(int B = 0;B < gnr();++B){//loop over the blocks PHM

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         a = ph2s[B][i][0];
         b = ph2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = ph2s[B][j][0];
            d = ph2s[B][j][1];

            //init
            (*this)(B,i,j) = 0.0;

            //first the S = 1/2 block of the PPHM matrix
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_de = 0;S_de < 2;++S_de){

                  ward = 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) ) * Tools::g6j(0,0,S,S_ab) * Tools::g6j(0,0,S,S_de);

                  for(int e = 0;e < Tools::gL() * Tools::gL();++e){

                     hard = ward * pphm.pph(0,S_ab,e,a,b,S_de,e,c,d);

                     //norms
                     if(e == a)
                        hard *= std::sqrt(2.0);

                     if(e == c)
                        hard *= std::sqrt(2.0);

                     (*this)(B,i,j) += hard;

                  }

               }

            //then the S = 3/2 block
            if(S == 1){

               for(int e = 0;e < Tools::gL() * Tools::gL();++e)
                  (*this)(B,i,j) += 4.0/3.0 * pphm.pph(1,1,e,a,b,1,e,c,d);

            }

         }
      }

   }

   this->symmetrize();

}

/* vim: set ts=3 sw=3 expandtab :*/
