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

double **DPM::_6j;

int **DPM::block_char;
int ***DPM::char_block;

int DPM::M;
int DPM::N;
int DPM::L;

/**
 * initialize the statics and allocate and construct all the lists
 * @param L_in dimension of the lattice
 * @param N_in nr of particles
 */
void DPM::init(int L_in,int N_in){

   L = L_in;
   N = N_in;

   M = L*L*2;

   dp2s = new vector< vector<int> > [M];

   s2dp = new int **** [M];

   for(int B = 0;B < M;++B){

      s2dp[B] = new int *** [2];

      for(int S = 0;S < 2;++S){

         s2dp[B][S] = new int ** [L*L];

         for(int a = 0;a < L*L;++a){

            s2dp[B][S][a] = new int * [L*L];

            for(int b = 0;b < L*L;++b)
               s2dp[B][S][a][b] = new int [L*L];

         }
      }
   }

   block_char = new int * [M];

   for(int B = 0;B < M;++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [L];

      for(int K_x = 0;K_x < L;++K_x)
         char_block[S][K_x] = new int [L];

   }

   int block = 0;

   vector<int> v(4);

   int dp;

   //loop over the blocks
   for(int K_x = 0;K_x < L;++K_x)
      for(int K_y = 0;K_y < L;++K_y){

         //S = 1/2

         block_char[block][0] = 0;//0 means spin 1/2
         block_char[block][1] = K_x;
         block_char[block][2] = K_y;

         char_block[0][K_x][K_y] = block;

         dp = 0;

         //first S = 1/2, S_ab = 0, a = b != c
         for(int a = 0;a < L*L;++a){

            for(int b = 0;b < a;++b){

               if( (2*Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L == K_x )
                  if( (2*Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L == K_y ){

                     v[0] = 0;//S_ab
                     v[1] = a;
                     v[2] = a;
                     v[3] = b;

                     dp2s[block].push_back(v);

                     s2dp[block][0][a][a][b] = dp;

                     ++dp;

                  }

            }

            for(int b = a + 1;b < L*L;++b){

               if( (2*Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L == K_x )
                  if( (2*Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L == K_y ){

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
         for(int a = 0;a < L*L;++a)
            for(int b = a + 1;b < L*L;++b)
               for(int c = b + 1;c < L*L;++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%L == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%L == K_y ){

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
         for(int a = 0;a < L*L;++a)
            for(int b = a + 1;b < L*L;++b)
               for(int c = b + 1;c < L*L;++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%L == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%L == K_y ){

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

         block_char[L*L + block][0] = 1;//0 means spin 1/2
         block_char[L*L + block][1] = K_x;
         block_char[L*L + block][2] = K_y;

         char_block[1][K_x][K_y] = L*L + block;

         dp = 0;

         //then S = 3/2, S_ab = 1, a < b < c
         for(int a = 0;a < L*L;++a)
            for(int b = a + 1;b < L*L;++b)
               for(int c = b + 1;c < L*L;++c){

                  if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%L == K_x )
                     if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%L == K_y ){

                        v[0] = 1;//S_ab
                        v[1] = a;
                        v[2] = b;
                        v[3] = c;

                        dp2s[L*L + block].push_back(v);

                        s2dp[L*L + block][1][a][b][c] = dp;

                        ++dp;

                     }

               }




         ++block;

      }

   //allocate
   _6j = new double * [2];

   for(int S = 0;S < 2;++S)
      _6j[S] = new double [2];

   //initialize
   _6j[0][0] = -0.5;
   _6j[0][1] = 0.5;
   _6j[1][0] = 0.5;
   _6j[1][1] = 1.0/6.0;

}

/**
 * deallocate the static lists
 */
void DPM::clear(){

   delete [] dp2s;

   for(int B = 0;B < M;++B){

      for(int S = 0;S < 2;++S){

         for(int a = 0;a < L*L;++a){

            for(int b = 0;b < L*L;++b)
               delete [] s2dp[B][S][a][b];

            delete [] s2dp[B][S][a];

         }

         delete [] s2dp[B][S];

      }

      delete [] s2dp[B];

   }

   delete [] s2dp;

   for(int B = 0;B < M;++B)
      delete [] block_char[B];

   delete [] block_char;

   for(int S = 0;S < 2;++S){

      for(int K_x = 0;K_x < L;++K_x)
         delete [] char_block[S][K_x];

      delete [] char_block[S];

   }

   delete [] char_block;

   for(int S = 0;S < 2;++S)
      delete [] _6j[S];

   delete [] _6j;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 * (L*L) blocks, (L*L) for S = 1/2 and (M/2) 3/2.
 */
DPM::DPM() : BlockMatrix(M) {

   //set the dimension of the blocks

   for(int B = 0;B < L*L;++B)//S = 1/2
      setMatrixDim(B,dp2s[B].size(),2);

   for(int B = L*L;B < 2*L*L;++B)//S = 3/2
      setMatrixDim(B,dp2s[B].size(),4);

}

/**
 * copy constructor: constructs BlockMatrix object with 2 * (M/2) blocks, on (M/2) for S=1/2 and (M/2) for S=3/2, and copies the content of the dpm_c blocks into it,
 * @param dpm_c DPM to be copied into (*this)
 */
DPM::DPM(const DPM &dpm_c) : BlockMatrix(dpm_c) {

}

/**
 * destructor: 
 */
DPM::~DPM(){

}

/**
 * @return number of particles
 */
int DPM::gN() const {

   return N;

}

/**
 * @return number of single particle oribals
 */
int DPM::gM() const{

   return M;

}

/**
 * @return dimension of the lattice
 */
int DPM::gL() const{

   return L;

}

ostream &operator<<(ostream &output,const DPM &dpm_p){

   for(int B = 0;B < dpm_p.gnr();++B){

      output << B << "\t(" << dpm_p.gS(B) << "," << dpm_p.gK_x(B) << "," << dpm_p.gK_y(B) << ")\t" << dpm_p.gdim(B) << "\t" << dpm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < dpm_p.gdim(B);++i)
         for(int j = 0;j < dpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << 

               dpm_p.dp2s[B][i][0] << "\t" << dpm_p.dp2s[B][i][1] << "\t" << dpm_p.dp2s[B][i][2] << "\t" << dpm_p.dp2s[B][i][3] << 

               "\t" << dpm_p.dp2s[B][j][0] << "\t" << dpm_p.dp2s[B][j][1] << "\t" << dpm_p.dp2s[B][j][2] << "\t" << dpm_p.dp2s[B][j][3] << "\t" << dpm_p(B,i,j) << endl;

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
 * DPM(B,S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z) = sum_S_ac (some terms dependent on spin) DPM(B,S_ac,k_a,k_c,k_b,S_de,k_d,k_e,k_z) etc...
 * @param B The block index
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the dp row index i of block B together with b, c and S_ab
 * @param b second sp index that forms the dp row index i of block B together with a, c and S_ab
 * @param c third sp index that forms the dp row index i of block B together with a, b and S_ab
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the dp column index j of block B together with e, z and S_de
 * @param e second sp index that forms the dp column index j of block B together with d, z and S_de
 * @param z third sp index that forms the dp column index j of block B together with d, e and S_de
 * @return the number on place DPM(B,i,j) with the right phase and forefactor.
 */
double DPM::operator()(int B,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   return (*this)(block_char[B][0],block_char[B][1],block_char[B][2],S_ab,a,b,c,S_de,d,e,z);

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * DPM(S,K,S_ab,k_a,k_b,k_c,S_de,k_d,k_e,k_z) = sum_S_ac (some terms dependent on spin) DPM(S,K,S_ac,k_a,k_c,k_b,S_de,k_d,k_e,k_z) etc...
 * @param S dp-spin block quantumnumber: S = 0 means S == 1/2 and S = 1 means S == 3/2 (for simplicity)
 * @param K_x dp x-momentum block quantumnumber: K_x = 0 -> L-1 
 * @param K_y dp y-momentum block quantumnumber: K_y = 0 -> L-1 
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
double DPM::operator()(int S,int K_x,int K_y,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   //check if the momentum is correct:
   if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%L != K_x)
      return 0.0;

   if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%L != K_y)
      return 0.0;

   if( (Hamiltonian::ga_xy(d,0) + Hamiltonian::ga_xy(e,0) + Hamiltonian::ga_xy(z,0))%L != K_x)
      return 0.0;

   if( (Hamiltonian::ga_xy(d,1) + Hamiltonian::ga_xy(e,1) + Hamiltonian::ga_xy(z,1))%L != K_y)
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

   if(B < L*L){//spin 1/2 block:

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
            coef[0] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * _6j[0][S_ab];

            //the S_ca == 1 part:
            i[1] = s2dp[B][1][c][min][max];
            coef[1] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * std::sqrt(3.0) * _6j[1][S_ab];

            return 2;

         }
         else if(c == min){//c == min < max: this will also be a 1 dim list, because S_ac can only be 0 if a == c.

            i[0] = s2dp[B][0][c][min][max];
            coef[0] = std::sqrt(2.0) * phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * _6j[0][S_ab];

            return 1;

         }
         else if(c < max){//min < c < max

            //S_ac == 0 part:
            i[0] = s2dp[B][0][min][c][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * _6j[0][S_ab];

            //S_ac == 1 part:
            i[1] = s2dp[B][1][min][c][max];
            coef[1] = - phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * std::sqrt(3.0) * _6j[1][S_ab];

            return 2;

         }
         else{// min < c == max: also a 1 dim list, S_bc can only be 0 if b == c

            i[0] = s2dp[B][0][max][c][min];
            coef[0] = phase * std::sqrt(2.0) * std::sqrt(2.0*S_ab + 1.0) *_6j[0][S_ab];

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
   SPM spm(C,tpm);

   double ward = 2.0*B*tpm.trace();

   int a,b,c,d,e,z;
   int S_ab,S_de;

   int K_x,K_y;

   int sign_ab,sign_de;

   double norm_ab,norm_de;

   double hard;

   //start with the S = 1/2 blocks, these are the most difficult:
   for(int B = 0;B < L*L;++B){

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

            hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * _6j[S_ab][S_de];

            //init
            (*this)(B,i,j) = 0.0;

            //the np + sp part
            if(i == j)
               (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

            //other parts are a bit more difficult.

            //tp(1)
            if(c == z)
               if(S_ab == S_de){

                  K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L;
                  K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L;

                  (*this)(B,i,j) += A * tpm(S_ab,K_x,K_y,a,b,d,e);

               }

            //tp(2)
            if(b == z){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(c,1))%L;

               if(a == c)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,K_x,K_y,a,c,d,e);
               else
                  (*this)(B,i,j) += A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,K_x,K_y,a,c,d,e);

            }

            //tp(3)
            if(a == z){

               K_x = (Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%L;

               if(b == c)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_de * hard * tpm(S_de,K_x,K_y,b,c,d,e);
               else
                  (*this)(B,i,j) += A * norm_ab * sign_de * hard * tpm(S_de,K_x,K_y,b,c,d,e);

            }

            //tp(4)
            if(c == e){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L;

               if(d == z)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,K_x,K_y,a,b,d,z);
               else
                  (*this)(B,i,j) += A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,K_x,K_y,a,b,d,z);

            }

            //tp(5)
            if(b == e){

               double hulp = 0.0;

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(c,1))%L;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,K_x,K_y,a,c,d,z);

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

               K_x = (Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%L;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,K_x,K_y,b,c,d,z);

               if(b == c)
                  hulp *= std::sqrt(2.0);

               if(d == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            }

            //tp(7)
            if(c == d){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L;

               if(e == z)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * hard * tpm(S_ab,K_x,K_y,a,b,e,z);
               else
                  (*this)(B,i,j) += A * norm_de * sign_ab * hard * tpm(S_ab,K_x,K_y,a,b,e,z);

            }

            //tp(8)
            if(b == d){

               double hulp = 0.0;

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(c,1))%L;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,K_x,K_y,a,c,e,z);

               if(a == c)
                  hulp *= std::sqrt(2.0);

               if(e == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            }

            //tp(9)
            if(a == d){

               K_x = (Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%L;

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,K_x,K_y,b,c,e,z);

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
   for(int B = L*L;B < M;++B){

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
            if(c == z){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L;

               (*this)(B,i,j) += A * tpm(1,K_x,K_y,a,b,d,e);

            }

            //tp(2)
            if(b == z){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(c,1))%L;

               (*this)(B,i,j) -= A * tpm(1,K_x,K_y,a,c,d,e);

            }

            //tp(4)
            if(c == e){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L;

               (*this)(B,i,j) -= A * tpm(1,K_x,K_y,a,b,d,z);

            }

            //tp(5)
            if(b == e){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(c,1))%L;

               (*this)(B,i,j) += A * tpm(1,K_x,K_y,a,c,d,z);

            }

            //tp(7)
            if(c == d){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1))%L;

               (*this)(B,i,j) += A * tpm(1,K_x,K_y,a,b,e,z);

            }

            //tp(8)
            if(b == d){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(c,1))%L;

               (*this)(B,i,j) -= A * tpm(1,K_x,K_y,a,c,e,z);

            }

            //tp(9)
            if(a == d){

               K_x = (Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%L;
               K_y = (Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%L;

               (*this)(B,i,j) += A * tpm(1,K_x,K_y,b,c,e,z);

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * The T1-map: maps a TPM object (tpm) on a DPM object (*this). 
 * @param tpm input TPM
 */

void DPM::T(const TPM &tpm){

   double a = 1.0;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->T(a,b,c,tpm);

}

/** 
 * The hat function maps a TPM object tpm to a DPM object (*this) so that bar(this) = tpm,
 * The inverse of the TPM::bar function. It is a T1-like map.
 * @param tpm input TPM
 */
void DPM::hat(const TPM &tpm){

   double a = 1.0/(M - 4.0);
   double b = 1.0/((M - 4.0)*(M - 3.0)*(M - 2.0));
   double c = 1.0/((M - 4.0)*(M - 3.0));

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
