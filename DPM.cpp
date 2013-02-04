#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using std::ostream;
using std::ofstream;
using std::cout;
using std::endl;
using std::vector;

#include "include.h"

vector< vector<int> >*DPM::dp2s;
int *****DPM::s2dp;

int **DPM::block_char;
int ***DPM::char_block;

/**
 * initialize the statics and allocate and construct all the lists
 */
void DPM::init(){

   dp2s = new vector< vector<int> > [Tools::gM()];

   s2dp = new int **** [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B){

      s2dp[B] = new int *** [2];

      for(int S = 0;S < 2;++S){

         s2dp[B][S] = new int ** [Tools::gL()*Tools::gL()];

         for(int a = 0;a < Tools::gL()*Tools::gL();++a){

            s2dp[B][S][a] = new int * [Tools::gL()*Tools::gL()];

            for(int b = 0;b < Tools::gL()*Tools::gL();++b)
               s2dp[B][S][a][b] = new int [Tools::gL()*Tools::gL()];

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

   vector<int> v(4);

   int dp;

   //loop over the blocks
   for(int K_x = 0;K_x < Tools::gL();++K_x)
      for(int K_y = 0;K_y < Tools::gL();++K_y){

         //S = 1/2

         block_char[block][0] = 0;//0 means spin 1/2
         block_char[block][1] = K_x;
         block_char[block][2] = K_y;

         char_block[0][K_x][K_y] = block;

         dp = 0;

         //first S = 1/2, S_ab = 0, a = b != c
         for(int a = 0;a < Tools::gL()*Tools::gL();++a){

            for(int b = 0;b < a;++b){

               if( (2*Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%Tools::gL() == K_x )
                  if( (2*Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%Tools::gL() == K_y ){

                     v[0] = 0;//S_ab
                     v[1] = a;
                     v[2] = a;
                     v[3] = b;

                     dp2s[block].push_back(v);

                     s2dp[block][0][a][a][b] = dp;

                     ++dp;

                  }

            }

            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b){

               if( (2*Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%Tools::gL() == K_x )
                  if( (2*Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%Tools::gL() == K_y ){

                     v[0] = 0;//S_ab
                     v[1] = a;
                     v[2] = a;
                     v[3] = b;

                     dp2s[block].push_back(v);

                     s2dp[block][0][a][a][b] = dp;

                     ++dp;

                  }

            }

         }

         //then S = 1/2, S_ab = 0, a < b < c
         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b)
               for(int c = b + 1;c < Tools::gL()*Tools::gL();++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == K_y ){

                        v[0] = 0;//S_ab
                        v[1] = a;
                        v[2] = b;
                        v[3] = c;

                        dp2s[block].push_back(v);

                        s2dp[block][0][a][b][c] = dp;

                        ++dp;

                     }

               }

         //then S = 1/2, S_ab = 1, a < b < c
         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b)
               for(int c = b + 1;c < Tools::gL()*Tools::gL();++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == K_y ){

                        v[0] = 1;//S_ab
                        v[1] = a;
                        v[2] = b;
                        v[3] = c;

                        dp2s[block].push_back(v);

                        s2dp[block][1][a][b][c] = dp;

                        ++dp;

                     }

               }

         //S = 3/2

         block_char[Tools::gL()*Tools::gL() + block][0] = 1;//0 means spin 1/2
         block_char[Tools::gL()*Tools::gL() + block][1] = K_x;
         block_char[Tools::gL()*Tools::gL() + block][2] = K_y;

         char_block[1][K_x][K_y] = Tools::gL()*Tools::gL() + block;

         dp = 0;

         //then S = 3/2, S_ab = 1, a < b < c
         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b)
               for(int c = b + 1;c < Tools::gL()*Tools::gL();++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == K_y ){

                        v[0] = 1;//S_ab
                        v[1] = a;
                        v[2] = b;
                        v[3] = c;

                        dp2s[Tools::gL()*Tools::gL() + block].push_back(v);

                        s2dp[Tools::gL()*Tools::gL() + block][1][a][b][c] = dp;

                        ++dp;

                     }

               }




         ++block;

      }

}

/**
 * deallocate the static lists
 */
void DPM::clear(){

   delete [] dp2s;

   for(int B = 0;B < Tools::gM();++B){

      for(int S = 0;S < 2;++S){

         for(int a = 0;a < Tools::gL()*Tools::gL();++a){

            for(int b = 0;b < Tools::gL()*Tools::gL();++b)
               delete [] s2dp[B][S][a][b];

            delete [] s2dp[B][S][a];

         }

         delete [] s2dp[B][S];

      }

      delete [] s2dp[B];

   }

   delete [] s2dp;

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
 * standard constructor: constructs BlockMatrix object with 2 * (Tools::gL()*Tools::gL()) blocks, (Tools::gL()*Tools::gL()) for S = 1/2 and (M/2) 3/2.
 */
DPM::DPM() : BlockMatrix(Tools::gM()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL()*Tools::gL();++B)//S = 1/2
      setMatrixDim(B,dp2s[B].size(),2);

   for(int B = Tools::gL()*Tools::gL();B < 2*Tools::gL()*Tools::gL();++B)//S = 3/2
      setMatrixDim(B,dp2s[B].size(),4);

}

/**
 * copy constructor: constructs BlockMatrix object with 2 * (M/2) blocks, on (M/2) for S=1/2 and (M/2) for S=3/2,
 * and copies the content of the dpm_c blocks into it.
 * @param dpm_c DPM to be copied into (*this)
 */
DPM::DPM(const DPM &dpm_c) : BlockMatrix(dpm_c) { }

/**
 * destructor: 
 */
DPM::~DPM(){ }

ostream &operator<<(ostream &output,const DPM &dpm_p){

   for(int B = 0;B < dpm_p.gnr();++B){

      output << B << "\t(" << dpm_p.gS(B) << "," << dpm_p.gK_x(B) << "," << dpm_p.gK_y(B) << ")\t" << dpm_p.gdim(B) << "\t" << dpm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < dpm_p.gdim(B);++i)
         for(int j = 0;j < dpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << 

               dpm_p.dp2s[B][i][0] << "\t" << dpm_p.dp2s[B][i][1] << "\t" << dpm_p.dp2s[B][i][2] << "\t" << dpm_p.dp2s[B][i][3] << 

               "\t" << dpm_p.dp2s[B][j][0] << "\t" << dpm_p.dp2s[B][j][1] << "\t" << 
               
               dpm_p.dp2s[B][j][2] << "\t" << dpm_p.dp2s[B][j][3] << "\t" << dpm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * @return the dp spin-index corresponding to the blockindex (0 for 1/2 and 1 for 3/2)
 */
int DPM::gS(int block) const{

   return block_char[block][0];

}

/**
 * @return the dp x-momentum index corresponding to the blockindex
 */
int DPM::gK_x(int block) const{

   return block_char[block][1];

}

/**
 * @return the dp y-momentum index corresponding to the blockindex
 */
int DPM::gK_y(int block) const{

   return block_char[block][2];

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * DPM(S,S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z) = sum_S_ac (some terms dependent on spin) DPM(S,K,S_ac,k_a,k_c,k_b,S_de,k_d,k_e,k_z) etc...
 * @param S dp-spin block quantumnumber: S = 0 means S == 1/2 and S = 1 means S == 3/2 (for simplicity)
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the dp row index i together with b, c and S_ab, with dp spin and momentum SK
 * @param b second sp index that forms the dp row index i together with a, c and S_ab, with dp spin and momentum SK
 * @param c third sp index that forms the dp row index i together with a, b and S_ab, with dp spin and momentum SK
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the dp column index j together with e, z and S_de, with dp spin and momentum SK
 * @param e second sp index that forms the dp column index j together with d, z and S_de, with dp spin and momentum SK
 * @param z third sp index that forms the dp column index j together with d, e and S_de, with dp spin and momentum SK
 * @return the number on place DPM(B,i,j) with the right phase and forefactor.
 */
double DPM::operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   int K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL();
   int K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL();

   //check if the momentum is correct:
   if( K_x != (Hamiltonian::ga_xy(d,0) + Hamiltonian::ga_xy(e,0) + Hamiltonian::ga_xy(z,0))%Tools::gL() )
      return 0.0;

   if( K_y != (Hamiltonian::ga_xy(d,1) + Hamiltonian::ga_xy(e,1) + Hamiltonian::ga_xy(z,1))%Tools::gL() )
      return 0.0;

   //blockindex:
   int B = char_block[S][K_x][K_y];

   int *i = new int [2];
   double *coef_i = new double [2];

   int dim_i = get_inco(B,S_ab,a,b,c,i,coef_i);

   if(dim_i == 0){

      delete [] i;
      delete [] coef_i;

      return 0.0;

   }

   int *j = new int [2];
   double *coef_j = new double [2];

   int dim_j = get_inco(B,S_de,d,e,z,j,coef_j);

   if(dim_j == 0){

      delete [] i;
      delete [] j;

      delete [] coef_i;
      delete [] coef_j;

      return 0.0;

   }

   double ward = 0.0;

   for(int I = 0;I < dim_i;++I)
      for(int J = 0;J < dim_j;++J)
         ward += coef_i[I] * coef_j[J] * (*this)(B,i[I],j[J]);

   delete [] i;
   delete [] j;

   delete [] coef_i;
   delete [] coef_j;

   return ward;

   return 0;

}

/** 
 * Static member function that gets the dp-indices and their coefficients of the (s and t )-p indices S_ab,a,b,c.
 * @param B block index of the state
 * @param S_ab intermediate spincoupling of a and b.
 * @param a first sp orbital
 * @param b second sp orbital
 * @param c third sp orbital
 * @param i pointer of dim 1 or 2 containing the indices occuring in the expansion of this particular dp state in the normal basis (a==b,c a < b < c).
 * @param coef pointer of dim 1 or 2 containing the coefficients occuring in the expansion.
 * @return the number of terms in the expansion (1 or 2), also the dim of pointers i and coef. When zero is returned this is not a valid element.
 */
int DPM::get_inco(int B,int S_ab,int a,int b,int c,int *i,double *coef) const{

   //they cannot all be equal
   if(a == b && b == c)
      return 0;

   if(B < Tools::gL()*Tools::gL()){//spin 1/2 block:

      //if normal basis:
      if(a == b){

         if(S_ab == 1)//spin has to be zero for a == b
            return 0;

         i[0] = s2dp[B][0][a][b][c];
         coef[0] = 1;

         return 1;

      }
      else if (a < b && b < c){

         i[0] = s2dp[B][S_ab][a][b][c];
         coef[0] = 1;

         return 1;

      }
      else{//anomal basis:

         int min,max,phase;

         //first order a and b for code saving reasons
         if(a < b){

            min = a;
            max = b;

            phase = 1;

         }
         else{

            min = b;
            max = a;

            phase = 1 - 2*S_ab;

            if(c > max){//we still have one simple dim = 1 term left: b < a < c

               i[0] = s2dp[B][S_ab][b][a][c];
               coef[0] = phase;

               return 1;

            }

         }

         //now we have four possibilities left:
         //don't forget to multiply every result by phase to get the right a and b for min and max!
         // 1) c < min < max
         // 2) c == min < max
         // 3) min < c < max
         // 4) min < max == c
         if(c < min){//c < min < max

            //the S_ca == 0 part:
            i[0] = s2dp[B][0][c][min][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            //the S_ca == 1 part:
            i[1] = s2dp[B][1][c][min][max];
            coef[1] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * std::sqrt(3.0) * Tools::g6j(0,0,1,S_ab);

            return 2;

         }
         else if(c == min){//c == min < max: this will also be a 1 dim list, because S_ac can only be 0 if a == c.

            i[0] = s2dp[B][0][c][min][max];
            coef[0] = std::sqrt(2.0) * phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            return 1;

         }
         else if(c < max){//min < c < max

            //S_ac == 0 part:
            i[0] = s2dp[B][0][min][c][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            //S_ac == 1 part:
            i[1] = s2dp[B][1][min][c][max];
            coef[1] = - phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * std::sqrt(3.0) * Tools::g6j(0,0,1,S_ab);

            return 2;

         }
         else{// min < c == max: also a 1 dim list, S_bc can only be 0 if b == c

            i[0] = s2dp[B][0][max][c][min];
            coef[0] = phase * std::sqrt(2.0) * std::sqrt(2.0*S_ab + 1.0) *Tools::g6j(0,0,0,S_ab);

            return 1;

         }

      }

   }
   else{//spin 3/2 block, totally antisymmetrical in the spatial sp orbs.

      //only S_ab == 1 can couple to 3/2's.
      if(S_ab == 0)
         return 0;

      //if any of the sp orbs are equal, antisymmetry leads to zero:
      if(a == b || b == c || c == a)
         return 0;

      if(a < b){

         if(b < c){//a < b < c

            i[0] = s2dp[B][1][a][b][c];
            coef[0] = 1;

         }
         else if(c < a){//c < a < b

            i[0] = s2dp[B][1][c][a][b];
            coef[0] = 1;

         }
         else{//a < c < b

            i[0] = s2dp[B][1][a][c][b];
            coef[0] = -1;

         }

      }
      else{//b < a

         if(a < c){//b < a < c

            i[0] = s2dp[B][1][b][a][c];
            coef[0] = -1;

         }
         else if(c < b){//c < b < a

            i[0] = s2dp[B][1][c][b][a];
            coef[0] = -1;

         }
         else{//b < c < a

            i[0] = s2dp[B][1][b][c][a];
            coef[0] = 1;

         }

      }

      return 1;

   }

}

/**
 * The spincoupled, translationally invariant T1-like (generalized T1) map: maps a TPM object (tpm) on a DPM object (*this)
 * @param A term before the tp part of the map
 * @param B term before the np part of the map
 * @param C term before the sp part of the map
 * @param tpm input TPM
 */
void DPM::T(double A,double B,double C,const TPM &tpm) {

   //make sp matrix out of tpm
   SPM spm;
   spm.bar(C,tpm);

   double ward = 2.0*B*tpm.trace();

   int a,b,c,d,e,z;
   int S_ab,S_de;

   int sign_ab,sign_de;

   double norm_ab,norm_de;

   double hard;

   //start with the S = 1/2 blocks, these are the most difficult:
   for(int B = 0;B < Tools::gL()*Tools::gL();++B){

      for(int i = 0;i < gdim(B);++i){

         S_ab = dp2s[B][i][0];

         a = dp2s[B][i][1];
         b = dp2s[B][i][2];
         c = dp2s[B][i][3];

         sign_ab = 1 - 2*S_ab;

         norm_ab = 1.0;

         if(a == b)
            norm_ab /= std::sqrt(2.0);

         for(int j = i;j < gdim(B);++j){

            S_de = dp2s[B][j][0];

            d = dp2s[B][j][1];
            e = dp2s[B][j][2];
            z = dp2s[B][j][3];

            sign_de = 1 - 2*S_de;

            norm_de = 1.0;

            if(d == e)
               norm_de /= std::sqrt(2.0);

            hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * Tools::g6j(0,0,S_ab,S_de);

            //init
            (*this)(B,i,j) = 0.0;

            //the np + sp part
            if(i == j)
               (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

            //other parts are a bit more difficult.

            //tp(1)
            if(c == z)
               if(S_ab == S_de)
                  (*this)(B,i,j) += A * tpm(S_ab,a,b,d,e);

            //tp(2)
            if(b == z){

               if(a == c)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);
               else
                  (*this)(B,i,j) += A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);

            }

            //tp(3)
            if(a == z){

               if(b == c)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);
               else
                  (*this)(B,i,j) += A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);

            }

            //tp(4)
            if(c == e){

               if(d == z)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);
               else
                  (*this)(B,i,j) += A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);

            }

            //tp(5)
            if(b == e){

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,a,c,d,z);

               //correct for norms of the tpm
               if(a == c)
                  hulp *= std::sqrt(2.0);

               if(d == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * hulp;

            }

            //tp(6)
            if(a == e){

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,b,c,d,z);

               if(b == c)
                  hulp *= std::sqrt(2.0);

               if(d == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            }

            //tp(7)
            if(c == d){

               if(e == z)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);
               else
                  (*this)(B,i,j) += A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);

            }

            //tp(8)
            if(b == d){

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,a,c,e,z);

               if(a == c)
                  hulp *= std::sqrt(2.0);

               if(e == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            }

            //tp(9)
            if(a == d){

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,b,c,e,z);

               if(b == c)
                  hulp *= std::sqrt(2.0);

               if(e == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            }

         }
      }

   }

   //then the S = 3/2 blocks, this should be easy, totally antisymmetrical 
   for(int B = Tools::gL()*Tools::gL();B < Tools::gM();++B){

      for(int i = 0;i < gdim(B);++i){

         a = dp2s[B][i][1];
         b = dp2s[B][i][2];
         c = dp2s[B][i][3];

         for(int j = i;j < gdim(B);++j){

            d = dp2s[B][j][1];
            e = dp2s[B][j][2];
            z = dp2s[B][j][3];

            (*this)(B,i,j) = 0.0;

            //np + sp part:
            if(i == j)
               (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

            //tp(1)
            if(c == z)
               (*this)(B,i,j) += A * tpm(1,a,b,d,e);

            //tp(2)
            if(b == z)
               (*this)(B,i,j) -= A * tpm(1,a,c,d,e);

            //tp(4)
            if(c == e)
               (*this)(B,i,j) -= A * tpm(1,a,b,d,z);

            //tp(5)
            if(b == e)
               (*this)(B,i,j) += A * tpm(1,a,c,d,z);

            //tp(7)
            if(c == d)
               (*this)(B,i,j) += A * tpm(1,a,b,e,z);

            //tp(8)
            if(b == d)
               (*this)(B,i,j) -= A * tpm(1,a,c,e,z);

            //tp(9)
            if(a == d)
               (*this)(B,i,j) += A * tpm(1,b,c,e,z);

         }
      }

   }

   this->symmetrize();

}

void DPM::T(int option,const TPM &tpm){

   double A = 1.0;
   double B = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double C = 1.0/(Tools::gN() - 1.0);

   if(option == 1){

      //make sp matrix out of tpm
      SPM spm;
      spm.bar(C,tpm);

      double ward = 2.0*B*tpm.trace();

      int a,b,c,d,e,z;
      int S_ab,S_de;

      int sign_ab,sign_de;

      double norm_ab,norm_de;

      double hard;

      //start with the S = 1/2 blocks, these are the most difficult:
      for(int B = 0;B < Tools::gL()*Tools::gL();++B){

         for(int i = 0;i < gdim(B);++i){

            S_ab = dp2s[B][i][0];

            a = dp2s[B][i][1];
            b = dp2s[B][i][2];
            c = dp2s[B][i][3];

            sign_ab = 1 - 2*S_ab;

            norm_ab = 1.0;

            if(a == b)
               norm_ab /= std::sqrt(2.0);

            for(int j = i;j < gdim(B);++j){

               S_de = dp2s[B][j][0];

               d = dp2s[B][j][1];
               e = dp2s[B][j][2];
               z = dp2s[B][j][3];

               sign_de = 1 - 2*S_de;

               norm_de = 1.0;

               if(d == e)
                  norm_de /= std::sqrt(2.0);

               hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * Tools::g6j(0,0,S_ab,S_de);

               //init
               (*this)(B,i,j) = 0.0;

               //the np + sp part
               if(i == j)
                  (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

               //other parts are a bit more difficult.

               //tp(1)
               if(c == z)
                  if(S_ab == S_de)
                     (*this)(B,i,j) += A * tpm(S_ab,a,b,d,e);

               //tp(2)
               if(b == z){

                  if(a == c)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);
                  else
                     (*this)(B,i,j) += A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);

               }

               //tp(3)
               if(a == z){

                  if(b == c)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);
                  else
                     (*this)(B,i,j) += A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);

               }

               //tp(4)
               if(c == e){

                  if(d == z)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);
                  else
                     (*this)(B,i,j) += A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);

               }

               //tp(5)
               if(b == e){

                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,a,c,d,z);

                  //correct for norms of the tpm
                  if(a == c)
                     hulp *= std::sqrt(2.0);

                  if(d == z)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * hulp;

               }

               //tp(6)
               if(a == e){

                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,b,c,d,z);

                  if(b == c)
                     hulp *= std::sqrt(2.0);

                  if(d == z)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

               }

               //tp(7)
               if(c == d){

                  if(e == z)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);
                  else
                     (*this)(B,i,j) += A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);

               }

               //tp(8)
               if(b == d){

                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,a,c,e,z);

                  if(a == c)
                     hulp *= std::sqrt(2.0);

                  if(e == z)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

               }

               //tp(9)
               if(a == d){

                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,b,c,e,z);

                  if(b == c)
                     hulp *= std::sqrt(2.0);

                  if(e == z)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

               }

            }
         }

      }

      //then the S = 3/2 blocks, this should be easy, totally antisymmetrical 
      for(int B = Tools::gL()*Tools::gL();B < Tools::gM();++B){

         for(int i = 0;i < gdim(B);++i){

            a = dp2s[B][i][1];
            b = dp2s[B][i][2];
            c = dp2s[B][i][3];

            for(int j = i;j < gdim(B);++j){

               d = dp2s[B][j][1];
               e = dp2s[B][j][2];
               z = dp2s[B][j][3];

               (*this)(B,i,j) = 0.0;

               //np + sp part:
               if(i == j)
                  (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

               //tp(1)
               if(c == z)
                  (*this)(B,i,j) += A * tpm(1,a,b,d,e);

               //tp(2)
               if(b == z)
                  (*this)(B,i,j) -= A * tpm(1,a,c,d,e);

               //tp(4)
               if(c == e)
                  (*this)(B,i,j) -= A * tpm(1,a,b,d,z);

               //tp(5)
               if(b == e)
                  (*this)(B,i,j) += A * tpm(1,a,c,d,z);

               //tp(7)
               if(c == d)
                  (*this)(B,i,j) += A * tpm(1,a,b,e,z);

               //tp(8)
               if(b == d)
                  (*this)(B,i,j) -= A * tpm(1,a,c,e,z);

               //tp(9)
               if(a == d)
                  (*this)(B,i,j) += A * tpm(1,b,c,e,z);

            }
         }

      }

      this->symmetrize();

   }
   else{

      //make sp matrix out of tpm
      SPM spm;
      spm.bar(C,tpm);

      double ward = 2.0*B*tpm.trace();

      int a,b,c,d,e,z;
      int S_ab,S_de;

      int sign_ab,sign_de;

      double norm_ab,norm_de;

      double hard;

      //start with the S = 1/2 blocks, these are the most difficult:
      for(int B = 0;B < Tools::gL()*Tools::gL();++B){

         for(int i = 0;i < gdim(B);++i){

            S_ab = dp2s[B][i][0];

            a = dp2s[B][i][1];
            b = dp2s[B][i][2];
            c = dp2s[B][i][3];

            sign_ab = 1 - 2*S_ab;

            norm_ab = 1.0;

            if(a == b)
               norm_ab /= std::sqrt(2.0);

            for(int j = i;j < gdim(B);++j){

               S_de = dp2s[B][j][0];

               d = dp2s[B][j][1];
               e = dp2s[B][j][2];
               z = dp2s[B][j][3];

               sign_de = 1 - 2*S_de;

               norm_de = 1.0;

               if(d == e)
                  norm_de /= std::sqrt(2.0);

               hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * Tools::g6j(0,0,S_ab,S_de);

               //init
               (*this)(B,i,j) = 0.0;

               //the np + sp part
               if(i == j)
                  (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

               //other parts are a bit more difficult.

               //tp(1)
               if(c == z)
                  if(S_ab == S_de)
                     (*this)(B,i,j) += A * tpm(S_ab,a,b,d,e);

               //tp(2)
               if(b == z){

                  if(a == c)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);
                  else
                     (*this)(B,i,j) += A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);

               }

               //tp(3)
               if(a == z){

                  if(b == c)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);
                  else
                     (*this)(B,i,j) += A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);

               }

               //tp(4)
               if(c == e){

                  if(d == z)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);
                  else
                     (*this)(B,i,j) += A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);

               }

               //tp(5)
               if(b == e){

                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,a,c,d,z);

                  //correct for norms of the tpm
                  if(a == c)
                     hulp *= std::sqrt(2.0);

                  if(d == z)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * hulp;

               }

               //tp(6)
               if(a == e){

                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,b,c,d,z);

                  if(b == c)
                     hulp *= std::sqrt(2.0);

                  if(d == z)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

               }

               //tp(7)
               if(c == d){

                  if(e == z)
                     (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);
                  else
                     (*this)(B,i,j) += A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);

               }

               //tp(8)
               if(b == d){

                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,a,c,e,z);

                  if(a == c)
                     hulp *= std::sqrt(2.0);

                  if(e == z)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

               }

               //tp(9)
               if(a == d){

                  double hulp = 0.0;

                  //sum over intermediate spin
                  for(int Z = 0;Z < 2;++Z)
                     hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,b,c,e,z);

                  if(b == c)
                     hulp *= std::sqrt(2.0);

                  if(e == z)
                     hulp *= std::sqrt(2.0);

                  (*this)(B,i,j) += A * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

               }

            }
         }

      }

      //then the S = 3/2 blocks, this should be easy, totally antisymmetrical 
      for(int B = Tools::gL()*Tools::gL();B < Tools::gM();++B){

         for(int i = 0;i < gdim(B);++i){

            a = dp2s[B][i][1];
            b = dp2s[B][i][2];
            c = dp2s[B][i][3];

            for(int j = i;j < gdim(B);++j){

               d = dp2s[B][j][1];
               e = dp2s[B][j][2];
               z = dp2s[B][j][3];

               (*this)(B,i,j) = 0.0;

               //np + sp part:
               if(i == j)
                  (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

               //tp(1)
               if(c == z)
                  (*this)(B,i,j) += A * tpm(1,a,b,d,e);

               //tp(2)
               if(b == z)
                  (*this)(B,i,j) -= A * tpm(1,a,c,d,e);

               //tp(4)
               if(c == e)
                  (*this)(B,i,j) -= A * tpm(1,a,b,d,z);

               //tp(5)
               if(b == e)
                  (*this)(B,i,j) += A * tpm(1,a,c,d,z);

               //tp(7)
               if(c == d)
                  (*this)(B,i,j) += A * tpm(1,a,b,e,z);

               //tp(8)
               if(b == d)
                  (*this)(B,i,j) -= A * tpm(1,a,c,e,z);

               //tp(9)
               if(a == d)
                  (*this)(B,i,j) += A * tpm(1,b,c,e,z);

            }
         }

      }

      this->symmetrize();

   }


}

/**
 * The T1-map: maps a TPM object (tpm) on a DPM object (*this). 
 * @param tpm input TPM
 */
void DPM::T(const TPM &tpm){

   double a = 1.0;
   double b = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double c = 1.0/(Tools::gN() - 1.0);

   this->T(a,b,c,tpm);

}

/** 
 * The hat function maps a TPM object tpm to a DPM object (*this) so that bar(this) = tpm,
 * The inverse of the TPM::bar function. It is a T1-like map.
 * @param tpm input TPM
 */
void DPM::hat(const TPM &tpm){

   double a = 1.0/(Tools::gM() - 4.0);
   double b = 1.0/((Tools::gM() - 4.0)*(Tools::gM() - 3.0)*(Tools::gM() - 2.0));
   double c = 1.0/((Tools::gM() - 4.0)*(Tools::gM() - 3.0));

   this->T(a,b,c,tpm);

}

/**
 * Output to file, to be read by the spin_pd program.
 * @param filename output file
 */
void DPM::out_sp(const char *filename) const{

   ofstream output(filename);
   output.precision(15);

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i)
         for(int j = i;j < gdim(B);++j)
            output << block_char[B][0] << "\t" << dp2s[B][i][0] << "\t" << dp2s[B][i][1] << "\t" << dp2s[B][i][2] << "\t" << dp2s[B][i][3] << "\t"

               << dp2s[B][j][0] << "\t" << dp2s[B][j][1] << "\t" << dp2s[B][j][2] << "\t" << dp2s[B][j][3] << "\t" << (*this)(B,i,j) << endl;

   }

}

/* vim: set ts=3 sw=3 expandtab :*/
