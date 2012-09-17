#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::cout;
using std::endl;

#include "include.h"

vector< vector<int> > *PPHM::pph2s;
int *****PPHM::s2pph;

int **PPHM::block_char;
int ***PPHM::char_block;

/**
 * allocate and initialize the statics
 */
void PPHM::init(){

   //allocate
   pph2s = new vector< vector<int> > [Tools::gM()];

   s2pph = new int **** [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B){

      s2pph[B] = new int *** [2];

      for(int Z = 0;Z < 2;++Z){

         s2pph[B][Z] = new int ** [Tools::gL()*Tools::gL()];

         for(int a = 0;a < Tools::gL()*Tools::gL();++a){

            s2pph[B][Z][a] = new int * [Tools::gL()*Tools::gL()];

            for(int b = 0;b < Tools::gL()*Tools::gL();++b)
               s2pph[B][Z][a][b] = new int [Tools::gL()*Tools::gL()];

         }
      }
   }

   block_char = new int * [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [Tools::gL()];

      for(int K_x = 0;K_x < Tools::gL();++K_x)
         char_block[S][K_x] = new int [Tools::gL()];

   }

   int block = 0;

   int pph;

   vector<int> v(4);

   for(int K_x = 0;K_x < Tools::gL();++K_x)
      for(int K_y = 0;K_y < Tools::gL();++K_y){

         //S = 1/2
         block_char[block][0] = 0;//means 1/2
         block_char[block][1] = K_x;
         block_char[block][2] = K_y;

         char_block[0][K_x][K_y] = block;

         pph = 0;

         //S_ab = 0
         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a;b < Tools::gL()*Tools::gL();++b)
               for(int c = 0;c < Tools::gL()*Tools::gL();++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == K_y ){

                        v[0] = 0;
                        v[1] = a;
                        v[2] = b;
                        v[3] = c;

                        pph2s[block].push_back(v);

                        s2pph[block][0][a][b][c] = pph;

                        ++pph;

                     }

               }

         //S_ab = 1
         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b)
               for(int c = 0;c < Tools::gL()*Tools::gL();++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == K_y ){

                        v[0] = 1;
                        v[1] = a;
                        v[2] = b;
                        v[3] = c;

                        pph2s[block].push_back(v);

                        s2pph[block][1][a][b][c] = pph;

                        ++pph;

                     }

               }

         //S = 3/2
         block_char[Tools::gL()*Tools::gL() + block][0] = 1;//means 3/2
         block_char[Tools::gL()*Tools::gL() + block][1] = K_x;
         block_char[Tools::gL()*Tools::gL() + block][2] = K_y;

         char_block[1][K_x][K_y] = Tools::gL()*Tools::gL() + block;

         pph = 0;

         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b)
               for(int c = 0;c < Tools::gL()*Tools::gL();++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == K_y ){

                        v[0] = 1;
                        v[1] = a;
                        v[2] = b;
                        v[3] = c;

                        pph2s[Tools::gL()*Tools::gL() + block].push_back(v);

                        s2pph[Tools::gL()*Tools::gL() + block][1][a][b][c] = pph;

                        ++pph;

                     }

               }

         ++block;

      }

}

/**
 * static function that deallocates the static lists.
 */
void PPHM::clear(){

   delete [] pph2s;

   for(int B = 0;B < Tools::gM();++B){

      for(int Z = 0;Z < 2;++Z){

         for(int a = 0;a < Tools::gL()*Tools::gL();++a){

            for(int b = 0;b < Tools::gL()*Tools::gL();++b)
               delete [] s2pph[B][Z][a][b];

            delete [] s2pph[B][Z][a];

         }

         delete [] s2pph[B][Z];

      }

      delete [] s2pph[B];

   }

   delete [] s2pph;

   for(int B = 0;B < Tools::gM();++B)
      delete [] block_char[B];

   delete [] block_char;

   for(int S = 0;S < 2;++S){

      for(int K_x = 0;K_x < Tools::gL();++K_x)
         delete [] char_block[S][K_x];

      delete [] char_block[S];

   }

   delete [] char_block; 

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 1/2 and 3/2.
 */
PPHM::PPHM() : BlockMatrix(Tools::gM()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL()*Tools::gL();++B)//S = 1/2
      setMatrixDim(B,pph2s[B].size(),2);

   for(int B = Tools::gL()*Tools::gL();B < 2*Tools::gL()*Tools::gL();++B)//S = 3/2
      setMatrixDim(B,pph2s[B].size(),4);

}

/**
 * copy constructor: constructs BlockMatrix object with M blocks, M/2 for S=1/2 and M/2 for S=3/2, and copies the content of the pphm_c blocks into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and pph basis.
 * @param pphm_c PPHM object to be copied into (*this)
 */
PPHM::PPHM(const PPHM &pphm_c) : BlockMatrix(pphm_c) { }

/**
 * Destructor, if counter = 1 the lists will be deallocated.
 */
PPHM::~PPHM(){ }

/**
 * @param block block index
 * @return the spin of the block.
 */
int PPHM::gS(int block) const{

   return block_char[block][0];

}

/**
 * @param block block index
 * @return the x-momentum of the block.
 */
int PPHM::gK_x(int block) const{

   return block_char[block][1];

}

/**
 * @param block block index
 * @return the y-momentum of the block.
 */
int PPHM::gK_y(int block) const{

   return block_char[block][2];

}

ostream &operator<<(ostream &output,const PPHM &pphm_p){

   for(int B = 0;B < pphm_p.gnr();++B){

      output << pphm_p.gS(B) << "\t" << pphm_p.gK_x(B) << "\t" << pphm_p.gK_y(B) << "\t" << pphm_p.gdim(B) << "\t" << pphm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < pphm_p.gdim(B);++i)
         for(int j = 0;j < pphm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << 

               pphm_p.pph2s[B][i][0] << "\t" << pphm_p.pph2s[B][i][1] << "\t" << pphm_p.pph2s[B][i][2] << "\t" << pphm_p.pph2s[B][i][3] << 

               "\t" << pphm_p.pph2s[B][j][0] << "\t" << pphm_p.pph2s[B][j][1] << "\t" << pphm_p.pph2s[B][j][2] << "\t" << pphm_p.pph2s[B][j][3] 

               << "\t" << pphm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * @param S The pphm-spin index, when == 0 then access the block S = 1/2, for spinindex == 1 we access the S = 3/2.
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the pph row index i together with b, c and S_ab in block B
 * @param b second sp index that forms the pph row index i together with a, c and S_ab in block B
 * @param c third sp index that forms the pph row index i together with a, b and S_ab in block B
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the pph column index j together with e, z and S_de in block B
 * @param e second sp index that forms the pph column index j together with d, z and S_de in block B
 * @param z third sp index that forms the pph column index j together with d, e and S_de in block B
 * @return the number on place PPHM(B,i,j) with the right phase.
 */
double PPHM::operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   int K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL();
   int K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL();

   //check the momentum
   if( K_x != (Hamiltonian::ga_xy(d,0) + Hamiltonian::ga_xy(e,0) + Hamiltonian::ga_xy(z,0))%Tools::gL() )
      return 0;

   if( K_y != (Hamiltonian::ga_xy(d,1) + Hamiltonian::ga_xy(e,1) + Hamiltonian::ga_xy(z,1))%Tools::gL() )
      return 0;

   int i,j;

   int phase_i = get_inco(char_block[S][K_x][K_y],S_ab,a,b,c,i);

   if(phase_i == 0)
      return 0;

   int phase_j = get_inco(char_block[S][K_x][K_y],S_de,d,e,z,j);

   if(phase_j == 0)
      return 0;

   return phase_i*phase_j* (*this)(char_block[S][K_x][K_y],i,j);

}

/** 
 * Member function that gets the pph-index and phase corresponding to the sp indices S, K, S_ab, k_a, k_b, k_c.
 * @param B the block index
 * @param S_ab intermediate spincoupling of k_a and k_b. = 0 or 1
 * @param a first sp orbital
 * @param b second sp orbital
 * @param c third sp orbital
 * @param i the corresponding pph index will be stored in this int after calling the function
 * @return the phase needed to get to a normal ordering of indices that corresponds to a pph index i
 */
int PPHM::get_inco(int B,int S_ab,int a,int b,int c,int &i) const{

   if(B < Tools::gL()*Tools::gL()){//S = 1/2

      if(S_ab == 0){//symmetric in spatial sp's

         if(a <= b)
            i = s2pph[B][0][a][b][c];
         else
            i = s2pph[B][0][b][a][c];

         return 1;

      }
      else{//antisymmetric in spatial sp's

         if(a == b)
            return 0;

         if(a < b){

            i = s2pph[B][1][a][b][c];

            return 1;

         }
         else{

            i = s2pph[B][1][b][a][c];

            return -1;

         }

      }

   }
   else{//S = 3/2

      if(S_ab == 0)//no possibile for S = 3/2
         return 0;

      if(a == b)//no possibile for S = 3/2
         return 0;

      if(a < b){

         i = s2pph[B][1][a][b][c];

         return 1;

      }
      else{

         i = s2pph[B][1][b][a][c];

         return -1;

      }

   }

}

/**
 * The spincoupled, translationally invariant T2 map, maps a TPM onto a PPHM object. See notes for more info
 * be aware that the c and z in the T2 notation are holes and transform in TPM space (remember the G-map)
 * @param tpm input TPM matrix
 */
void PPHM::T(const TPM &tpm){

   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   int a,b,c,d,e,z;
   int S_ab,S_de;

   double norm_ab,norm_de;
   int sign_ab,sign_de;

   //first the S = 1/2 blocks, these should be the most difficult ones.
   for(int B = 0;B < Tools::gL()*Tools::gL();++B){

      for(int i = 0;i < gdim(B);++i){

         S_ab = pph2s[B][i][0];

         a = pph2s[B][i][1];
         b = pph2s[B][i][2];

         //change to tp-notation:
         c = Hamiltonian::bar(pph2s[B][i][3]);

         sign_ab = 1 - 2*S_ab;

         norm_ab = 1.0;

         if(a == b)
            norm_ab /= std::sqrt(2.0);

         for(int j = i;j < gdim(B);++j){

            S_de = pph2s[B][j][0];

            d = pph2s[B][j][1];
            e = pph2s[B][j][2];

            //change to tp-notation:
            z = Hamiltonian::bar(pph2s[B][j][3]);

            sign_de = 1 - 2*S_de;

            norm_de = 1.0;

            if(d == e)
               norm_de /= std::sqrt(2.0);

            //start the map: init
            (*this)(B,i,j) = 0.0;

            //sp term becomes diagonal here:
            if(i == j)
               (*this)(B,i,j) += spm[c];

            //tp(1)
            if(c == z)
               if(S_ab == S_de)
                  (*this)(B,i,j) += tpm(S_ab,a,b,d,e);

            //tp(2)
            if(a == d){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,e,z,b);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == e)
                  ward *= std::sqrt(2.0);

               if(z == b)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= ward;

            }

            //tp(3)
            if(b == d){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,e,z,a);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == e)
                  ward *= std::sqrt(2.0);

               if(z == a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_ab * ward;

            }

            //tp(4)
            if(a == e){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,d,z,b);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == d)
                  ward *= std::sqrt(2.0);

               if(z == b)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_de * ward;

            }

            //tp(5)
            if(b == e){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,d,z,a);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == d)
                  ward *= std::sqrt(2.0);

               if(z == a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_ab * sign_de * ward;

            }

         }

      }

   }

   //the easier S = 3/2 part:
   for(int B = Tools::gL()*Tools::gL();B < Tools::gM();++B){

      for(int i = 0;i < gdim(B);++i){

         a = pph2s[B][i][1];
         b = pph2s[B][i][2];

         //change to correct sp-momentum
         c = Hamiltonian::bar(pph2s[B][i][3]);

         for(int j = i;j < gdim(B);++j){

            d = pph2s[B][j][1];
            e = pph2s[B][j][2];

            //change to correct sp-momentum
            z = Hamiltonian::bar(pph2s[B][j][3]);

            //init
            (*this)(B,i,j) = 0.0;

            //sp part is diagonal
            if(i == j)
               (*this)(B,i,j) += spm[c];

            //tp(1)
            if(c == z)
               (*this)(B,i,j) += tpm(1,a,b,d,e);

            //tp(2)
            if(a == d){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,e,z,b);

               if(c == e)
                  ward *= std::sqrt(2.0);

               if(z == b)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= ward;

            }

            //tp(3)
            if(b == d){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,e,z,a);

               if(c == e)
                  ward *= std::sqrt(2.0);

               if(z == a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) += ward;

            }

            //tp(5)
            if(b == e){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,d,z,a);

               if(c == d)
                  ward *= std::sqrt(2.0);

               if(z == a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= ward;

            }

         }

      }

   }

   this->symmetrize();

}

/**
 * Output to file, to be read by the spin_pd program.
 * @param filename output file
 */
void PPHM::out_sp(const char *filename) const{

   ofstream output(filename);
   output.precision(15);

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i)
         for(int j = i;j < gdim(B);++j)
            output << block_char[B][0] << "\t" << pph2s[B][i][0] << "\t" << pph2s[B][i][1] << "\t" << pph2s[B][i][2] << "\t" << pph2s[B][i][3] << "\t"

               << pph2s[B][j][0] << "\t" << pph2s[B][j][1] << "\t" << pph2s[B][j][2] << "\t" << pph2s[B][j][3] << "\t" << (*this)(B,i,j) << endl;

   }

}

/* vim: set ts=3 sw=3 expandtab :*/
