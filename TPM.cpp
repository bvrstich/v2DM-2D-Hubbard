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

vector< vector<int> > *TPM::t2s;

int ***TPM::s2t;

int **TPM::block_char;
int ***TPM::char_block;

double **TPM::_6j;

int TPM::N;
int TPM::L;
int TPM::M;

double TPM::Sa = 1.0;
double TPM::Sc = 0.0;

/**
 * static function that initializes the static variables
 * @param L_in dimension of the lattice
 * @param N_in nr of particles
 */
void TPM::init(int L_in,int N_in){

   L = L_in;
   N = N_in;
   M = L*L*2;

   //allocate stuff
   t2s = new vector< vector<int> > [M];

   s2t = new int ** [M];

   for(int B = 0;B < M;++B){

      s2t[B] = new int * [L*L];

      for(int a = 0;a < L*L;++a)
         s2t[B][a] = new int [L*L];

   }

   block_char = new int * [M];

   for(int B = 0;B < M;++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [L];

      for(int x = 0;x < L;++x)
         char_block[S][x] = new int [L];

   }

   vector<int> v(2);

   int block = 0;

   //tp index
   int t;

   //loop over the K_x K_y blocks
   for(int K_x = 0;K_x < L;++K_x)
      for(int K_y = 0;K_y < L;++K_y){

         t = 0;

         //S = 0
         block_char[block][0] = 0;
         block_char[block][1] = K_x;
         block_char[block][2] = K_y;

         char_block[0][K_x][K_y] = block;

         for(int a = 0;a < L*L;++a)
            for(int b = a;b < L*L;++b){

               if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % L == K_x ) && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % L == K_y ) ){

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
         block_char[L*L + block][0] = 1;
         block_char[L*L + block][1] = K_x;
         block_char[L*L + block][2] = K_y;

         char_block[1][K_x][K_y] = L*L + block;

         for(int a = 0;a < L*L;++a)
            for(int b = a + 1;b < L*L;++b){

               if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % L == K_x ) && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % L == K_y ) ){

                  v[0] = a;
                  v[1] = b;

                  t2s[L*L + block].push_back(v);

                  s2t[L*L + block][a][b] = t;
                  s2t[L*L + block][b][a] = t;

                  ++t;

               }

            }

         ++block;

      }

   //allocate 6j
   _6j = new double * [2];

   for(int S = 0;S < 2;++S)
      _6j[S] = new double [2]; 

   //initialize
   _6j[0][0] = -0.5;
   _6j[0][1] = 0.5;
   _6j[1][0] = 0.5;
   _6j[1][1] = 1.0/6.0;

   //construct the parametrs needed for the overlapmatrix
   constr_overlap();

}

/**
 * Fucntion that construct the variables needed for the overlapmatrix
 */
void TPM::constr_overlap(){

   Sa += 1.0;
   Sc += (2.0*N - M)/((N - 1.0)*(N - 1.0));

}

/**
 * static function that deallocates the static variables
 */
void TPM::clear(){

   delete [] t2s;

   for(int B = 0;B < M;++B){

      for(int a = 0;a < L*L;++a)
         delete [] s2t[B][a];

      delete [] s2t[B];

      delete [] block_char[B];

   }

   delete [] s2t;

   delete [] block_char;

   for(int S = 0;S < 2;++S)
      delete [] _6j[S];

   delete [] _6j;

   for(int S = 0;S < 2;++S){

      for(int x = 0;x < L;++x)
         delete [] char_block[S][x];

      delete [] char_block[S];

   }

   delete [] char_block;

}

/**
 * standard constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice: 
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 */
TPM::TPM() : BlockMatrix(L*L*2) {

   //set the dimension of the blocks

   for(int B = 0;B < L*L;++B)//S = 0
      setMatrixDim(B,t2s[B].size(),1);

   for(int B = L*L;B < 2*L*L;++B)//S = 1
      setMatrixDim(B,t2s[B].size(),3);

}

/**
 * copy constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice:
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 * @param tpm_c The TPM object to be copied into (*this)
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){

}

/**
 * destructor
 */
TPM::~TPM(){

}

ostream &operator<<(ostream &output,const TPM &tpm_p){

   int S,K_x,K_y;

   for(int B = 0;B < tpm_p.gnr();++B){

      S = tpm_p.block_char[B][0];
      K_x = tpm_p.block_char[B][1];
      K_y = tpm_p.block_char[B][2];

      output << "S =\t" << S << "\tK_x =\t" << K_x << "\tK_y =\t" << K_y << "\tdimension =\t" << tpm_p.gdim(B) << "\tdegeneracy =\t" << tpm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < tpm_p.gdim(B);++i)
         for(int j = 0;j < tpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << tpm_p.t2s[B][i][0] << "\t" << tpm_p.t2s[B][i][1]

               << "\t" << tpm_p.t2s[B][j][0] << "\t" << tpm_p.t2s[B][j][1] << "\t" << tpm_p(B,i,j) << endl;

         }

      output << std::endl;

   }

   return output;

}

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param B The block-index
 * @param a first sp index that forms the tp row index i of block B, together with b
 * @param b second sp index that forms the tp row index i of block B, together with a
 * @param c first sp index that forms the tp column index j of block B, together with d
 * @param d second sp index that forms the tp column index j of block B, together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int B,int a,int b,int c,int d) const{

   int S = block_char[B][0];
   int K_x = block_char[B][1];
   int K_y = block_char[B][2];

   //check if momentum is ok
   if( ( Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) )%L != K_x)
      return 0;

   if( ( Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) )%L != K_y)
      return 0;

   if( ( Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(d,0) )%L != K_x)
      return 0;

   if( ( Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(d,1) )%L != K_y)
      return 0;

   if(S == 0){

      int i = s2t[B][a][b];
      int j = s2t[B][c][d];

      return (*this)(B,i,j);

   }
   else{

      if( (a == b) || (c == d) )
         return 0;
      else{

         int i = s2t[B][a][b];
         int j = s2t[B][c][d];

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase * (*this)(B,i,j);

      }

   }

}
/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S The tp spin quantumnumber
 * @param K_x The tp x-momentum
 * @param K_y The tp y-momentum
 * @param a first sp index that forms the tp row index i of block B(S,K_x,K_y), together with b
 * @param b second sp index that forms the tp row index i of block B(S,K_x,K_y), together with a
 * @param c first sp index that forms the tp column index j of block B(S,K_x,K_y), together with d
 * @param d second sp index that forms the tp column index j of block B(S,K_x,K_y), together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int K_x,int K_y,int a,int b,int c,int d) const{

   //check if momentum is ok
   if( ( Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) )%L != K_x)
      return 0;

   if( ( Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) )%L != K_y)
      return 0;

   if( ( Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(d,0) )%L != K_x)
      return 0;

   if( ( Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(d,1) )%L != K_y)
      return 0;

   int B = char_block[S][K_x][K_y];

   if(S == 0){

      int i = s2t[B][a][b];
      int j = s2t[B][c][d];

      return (*this)(B,i,j);

   }
   else{//S = 1

      if( (a == b) || (c == d) )
         return 0;
      else{

         int i = s2t[B][a][b];
         int j = s2t[B][c][d];

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase * (*this)(B,i,j);

      }

   }

}

/**
 * @return number of particles
 */
int TPM::gN() const{

   return N;

}

/**
 * @return number of sp orbitals
 */
int TPM::gM() const{

   return M;

}

/**
 * @return dimension of the lattice
 */
int TPM::gL() const{

   return L;

}

/**
 * construct the spinsymmetrical hubbard hamiltonian in momentum space with on site repulsion U
 * @param U onsite repulsion term
 */
void TPM::hubbard(double U){

   int a,b,c,d;//sp indices

   double ward = 1.0/(N - 1.0);

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            //init
            (*this)(B,i,j) = 0;

            //hopping (kinetic energy):
            if(i == j){

               (*this)(B,i,i) = -2.0 * ward * ( cos( 2.0 * Hamiltonian::ga_xy(a,0) * 3.141592653589793238462 / (double) L)  

                     + cos( 2.0 * Hamiltonian::ga_xy(a,1) * 3.141592653589793238462 / (double) L) 

                     + cos( 2.0 * Hamiltonian::ga_xy(b,0) * 3.141592653589793238462 / (double) L) 

                     + cos( 2.0 * Hamiltonian::ga_xy(b,1) * 3.141592653589793238462 / (double) L) );

            }

            //on-site repulsion
            if(B < L*L){//only spin zero

               double hard = 2.0*U / (double) (L*L);

               if(a == b)
                  hard /= std::sqrt(2.0);

               if(c == d)
                  hard /= std::sqrt(2.0);

               (*this)(B,i,j) += hard;

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * Output to a "non-translationally invariant" file, for input in spin_pd.
 * @param filename name and location of the file you want to print to.
 */
void TPM::out_sp(const char *filename) const{

   ofstream output(filename);
   output.precision(15);

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i)
         for(int j = 0;j < gdim(B);++j)
            output << block_char[B][0] << "\t" << t2s[B][i][0] << "\t" << t2s[B][i][1] << "\t" << t2s[B][j][0] << "\t" << t2s[B][j][1] 

               << "\t" << (*this)(B,i,j) << endl;

   }

}


/**
 * The spincoupled Q map
 * @param option = 1, regular Q map , = -1 inverse Q map
 * @param tpm_d the TPM of which the Q map is taken and saved in this.
 */
void TPM::Q(int option,const TPM &tpm_d){

   double a = 1;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * The spincoupled Q-like map: see primal-dual.pdf for more info (form: Q^S(A,B,C)(TPM) )
 * @param option = 1, regular Q-like map , = -1 inverse Q-like map
 * @param A factor in front of the two particle piece of the map
 * @param B factor in front of the no particle piece of the map
 * @param C factor in front of the single particle piece of the map
 * @param tpm_d the TPM of which the Q-like map is taken and saved in this.
 */
void TPM::Q(int option,double A,double B,double C,const TPM &tpm_d){

   //for inverse
   if(option == -1){

      B = (B*A + B*C*M - 2.0*C*C)/( A * (C*(M - 2.0) -  A) * ( A + B*M*(M - 1.0) - 2.0*C*(M - 1.0) ) );
      C = C/(A*(C*(M - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm(C,tpm_d);

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   int a,b;

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         //tp part is only nondiagonal part
         for(int j = i;j < gdim(B);++j)
            (*this)(B,i,j) = A * tpm_d(B,i,j);

         (*this)(B,i,i) += ward - spm[a] - spm[b];

      }

   }

   this->symmetrize();

}

/**
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = N*(N - 1.0)/(M*(M - 1.0));

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         (*this)(B,i,i) = ward;

         for(int j = i + 1;j < gdim(B);++j)
            (*this)(B,i,j) = (*this)(B,j,i) = 0.0;

      }
   }

}

/**
 * orthogonal projection onto the space of traceless matrices
 */
void TPM::proj_Tr(){

   double ward = (2.0 * this->trace())/(M*(M - 1));

   this->min_unit(ward);

}

/**
 * Primal hessian map:\n\n
 * Hb = D_1 b D_1 + D_2 Q(b) D_2 + D_3 G(b) D_3 + D_4 T1(b) D_4 + D_5 T2(b) D5 \n\n
 * with D_1, D_2, D_3, D_4 and D_5 the P, Q, G, T1 and T2 blocks of the SUP D. 
 * @param b TPM domain matrix, hessian will act on it and the image will be put in this
 * @param D SUP matrix that defines the structure of the hessian map. (see primal-dual.pdf for more info)
 */
void TPM::H(const TPM &b,const SUP &D){

   this->L_map(D.tpm(0),b);

   //maak Q(b)
   TPM Qb;
   Qb.Q(1,b);

   TPM hulp;

   hulp.L_map(D.tpm(1),Qb);

   Qb.Q(1,hulp);

   *this += Qb;

   this->proj_Tr();

}

/**
 * Implementation of a linear conjugate gradient algoritm for the solution of the primal Newton equations\n\n
 * H(*this) =  b\n\n 
 * in which H represents the hessian map.
 * @param b righthandside of the equation
 * @param D SUP matrix that defines the structure of the hessian
 * @return return number of iterations needed to converge to the desired accuracy
 */
int TPM::solve(TPM &b,const SUP &D){

   *this = 0;

   //de r initialiseren op b
   TPM r(b);

   double rr = r.ddot(r);
   double rr_old,ward;

   //nog de Hb aanmaken ook, maar niet initialiseren:
   TPM Hb;

   int cg_iter = 0;

   while(rr > 1.0e-15){

      ++cg_iter;

      Hb.H(b,D);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      b.dscal(rr/rr_old);

      //dan r er bijtellen
      b += r;

   }

   return cg_iter;

}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param option = 1 direct overlapmatrix-map is used , = -1 inverse overlapmatrix map is used
 * @param tpm_d the input TPM
 */

void TPM::S(int option,const TPM &tpm_d){

   this->Q(option,Sa,0,Sc,tpm_d);

}

/**
 * Deduct the unitmatrix times a constant (scale) from this.\n\n
 * this -= scale* 1
 * @param scale the constant
 */
void TPM::min_unit(double scale){

   for(int B = 0;B < gnr();++B)
      for(int i = 0;i < gdim(B);++i)
         (*this)(B,i,i) -= scale;

}

/**
 * Collaps a SUP matrix S onto a TPM matrix like this:\n\n
 * sum_i Tr (S u^i)f^i = this
 * @param option = 0, project onto full symmetric matrix space, = 1 project onto traceless symmetric matrix space
 * @param S input SUP
 */

void TPM::collaps(int option,const SUP &S){

   *this = S.tpm(0);

   TPM hulp;

   hulp.Q(1,S.tpm(1));

   *this += hulp;

   if(option == 1)
      this->proj_Tr();

}

/**
 * @return The expectation value of the total spin for the TPM.
 */

double TPM::spin() const{

   double ward = 0.0;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      if(S == 0){

         for(int i = 0;i < gdim(B);++i)
            ward += -1.5 * (N - 2.0)/(N - 1.0) * (*this)(B,i,i);

      }
      else{

         for(int i = 0;i < this->gdim(B);++i)
            ward += 3.0 * ( -1.5 * (N - 2.0)/(N - 1.0) + 2.0 ) * (*this)(B,i,i);

      }

   }

   return ward;

}


/**
 * Fill a TPM object from a file.
 * @param input The ifstream object, corresponding to the file containing the TPM
 */
void TPM::in(ifstream &input){

   double block,dim,deg;
   int I,J;

   for(int B = 0;B < gnr();++B){

      input >> block >> dim >> deg;

      for(int i = 0;i < gdim(B);++i)
         for(int j = 0;j < gdim(B);++j)
            input >> I >> J >> (*this)(B,i,j);

   }

}

/**
 * Fill (*this) with the S^2 operator
 */

void TPM::set_S_2(){

   *this = 0.0;

   //S = 0 blocks
   for(int B = 0;B < L*L;++B)
      for(int i = 0;i < gdim(B);++i)
         (*this)(B,i,i) = -1.5 * (N - 2.0)/(N - 1.0);

   //S = 1 blocks
   for(int B = L*L;B < M;++B)
      for(int i = 0;i < gdim(B);++i)
         (*this)(B,i,i) = -1.5 * (N - 2.0)/(N - 1.0) + 2.0;

}
