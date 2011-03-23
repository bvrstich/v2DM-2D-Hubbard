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

#ifdef __G_CON

   Sa += 4.0;
   Sc += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __T1_CON

   Sa += M - 4.0;
   Sc -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );

#endif

#ifdef __T2_CON

   Sa += 5.0*M - 8.0;
   Sc += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));

#endif

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

#pragma omp parallel for private(a,b)
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

#ifdef __G_CON

   //maak G(b)
   PHM Gb;
   Gb.G(b);

   PHM hulpje;

   hulpje.L_map(D.phm(),Gb);

   hulp.G(hulpje);

   *this += hulp;

#endif

#ifdef __T1_CON

   DPM T1b;
   T1b.T(b);

   DPM hulp_T1;

   hulp_T1.L_map(D.dpm(),T1b);

   hulp.T(hulp_T1);

   *this += hulp;

#endif

#ifdef __T2_CON

   PPHM T2b;
   T2b.T(b);

   PPHM hulp_T2;

   hulp_T2.L_map(D.pphm(),T2b);

   hulp.T(hulp_T2);

   *this += hulp;

#endif

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

   while(rr > 1.0e-10){

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

#pragma omp parallel for
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

#ifdef __G_CON

   hulp.G(S.phm());

   *this += hulp;

#endif

#ifdef __T1_CON

   hulp.T(S.dpm());

   *this += hulp;

#endif

#ifdef __T2_CON

   hulp.T(S.pphm());

   *this += hulp;

#endif

   if(option == 1)
      this->proj_Tr();

}

/**
 * @return The expectation value of the total spin for the TPM.
 */

double TPM::spin() const{

   double ward = 0.0;

   int S;

#pragma omp parallel for private(S) reduction(+:ward)
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

/**
 * The G down map, maps a PHM object onto a TPM object using the G map.
 * @param phm input PHM
 */
void TPM::G(const PHM &phm){

   SPM spm(1.0/(N - 1.0),phm);

   int a,b,c,d;

   //the conjugated indices
   int a_,b_,c_,d_;

   int S,K_x,K_y;

   int sign;

#pragma omp parallel for private(a,b,c,d,a_,b_,c_,d_,S,K_x,K_y,sign)
   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      sign = 1 - 2*S;

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         a_ = Hamiltonian::bar(a);
         b_ = Hamiltonian::bar(b);

         //tp part is only nondiagonal part
         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            c_ = Hamiltonian::bar(c);
            d_ = Hamiltonian::bar(d);

            (*this)(B,i,j) = 0.0;

            //four ph terms:
            //1)
            K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(d_,0))%L;
            K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(d_,1))%L;

            for(int Z = 0;Z < 2;++Z)
               (*this)(B,i,j) -= (2.0*Z + 1.0) * _6j[S][Z] * phm(Z,K_x,K_y,a,d_,c,b_);

            //2)
            K_x = (Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c_,0))%L;
            K_y = (Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c_,1))%L;

            for(int Z = 0;Z < 2;++Z)
               (*this)(B,i,j) -= (2.0*Z + 1.0) * _6j[S][Z] * phm(Z,K_x,K_y,b,c_ ,d, a_);

            //3)
            K_x = (Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(d_,0))%L;
            K_y = (Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(d_,1))%L;

            for(int Z = 0;Z < 2;++Z)
               (*this)(B,i,j) -= sign * (2.0*Z + 1.0) * _6j[S][Z] * phm(Z,K_x,K_y,b,d_ ,c, a_ );

            //4)
            K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(c_,0))%L;
            K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(c_,1))%L;

            for(int Z = 0;Z < 2;++Z)
               (*this)(B,i,j) -= sign * (2.0*Z + 1.0) * _6j[S][Z] * phm(Z,K_x,K_y,a,c_,d,b_);

            //norm:
            if(a == b)
               (*this)(B,i,j) /= std::sqrt(2.0);

            if(c == d)
               (*this)(B,i,j) /= std::sqrt(2.0);

         }

         (*this)(B,i,i) += spm[a] + spm[b];

      }

   }

   this->symmetrize();

}

/**
 * Construct a spincoupled, translationally invariant TPM matrix out of a spincoupled, translationally invariant DPM matrix.
 * For the definition and derivation see symmetry.pdf
 * @param dpm input DPM
 */
void TPM::bar(const DPM &dpm){

   int a,b,c,d;
   int K_x,K_y;

   double ward;

#pragma omp parallel
{
   //first the S = 0 part, easiest:
#pragma omp for private(a,b,c,d,K_x,K_y,ward) nowait
   for(int B = 0;B < L*L;++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            //only total S = 1/2 can remain because cannot couple to S = 3/2 with intermediate S = 0
            for(int e = 0;e < L*L;++e){

               K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(e,0))%L;
               K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(e,1))%L;

               (*this)(B,i,j) += 2.0 * dpm(0,K_x,K_y,0,a,b,e,0,c,d,e);

            }

         }
      }

   }

   //then the S = 1 part:
#pragma omp for private(a,b,c,d,K_x,K_y,ward) nowait
   for(int B = M/2;B < M;++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            for(int Z = 0;Z < 2;++Z){//loop over the dpm blocks: S = 1/2 and 3/2 = Z + 1/2

               ward = 0.0;

               for(int e = 0;e < L*L;++e){

                  K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(e,0))%L;
                  K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(e,1))%L;

                  ward += dpm(Z,K_x,K_y,1,a,b,e,1,c,d,e);

               }

               ward *= (2 * (Z + 0.5) + 1.0)/3.0;

               (*this)(B,i,j) += ward;

            }

         }
      }

   }
}

   this->symmetrize();

}

/** 
 * The T1-down map that maps a DPM on TPM. This is just a Q-like map using the TPM::bar (dpm) as input.
 * @param dpm the input DPM matrix
 */
void TPM::T(const DPM &dpm){

   TPM tpm;
   tpm.bar(dpm);

   double a = 1;
   double b = 1.0/(3.0*N*(N - 1.0));
   double c = 0.5/(N - 1.0);

   this->Q(1,a,b,c,tpm);

}

/**
 * The spincoupled T2-down map that maps a PPHM on a TPM object.
 * @param pphm Input PPHM object
 */
void TPM::T(const PPHM &pphm){

   //first make the bar tpm
   TPM tpm;
   tpm.bar(pphm);

   //then make the bar phm
   PHM phm;
   phm.bar(pphm);

   //also make the bar spm with the correct scale factor
   SPM spm;
   spm.bar(0.5/(N - 1.0),pphm);

   int a,b,c,d;
   int a_,b_,c_,d_;
   int sign;

   double norm;

   int S;
   int K_x,K_y;

#pragma omp parallel for private(a,b,c,d,a_,b_,c_,d_,sign,norm,S,K_x,K_y)
   for(int B = 0;B < gnr();++B){//loop over the blocks

      S = block_char[B][0];

      sign = 1 - 2*S;

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         //and for access to the hole elements:
         a_ = Hamiltonian::bar(a);
         b_ = Hamiltonian::bar(b);

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            c_ = Hamiltonian::bar(c);
            d_ = Hamiltonian::bar(d);

            //determine the norm for the basisset
            norm = 1.0;

            if(S == 0){

               if(a == b)
                  norm /= std::sqrt(2.0);

               if(c == d)
                  norm /= std::sqrt(2.0);

            }

            //first the tp part
            (*this)(B,i,j) = tpm(B,i,j);

            //sp part is diagonal for translationaly invariance
            if(i == j)
               (*this)(B,i,j) += spm[a_] + spm[b_];

            for(int Z = 0;Z < 2;++Z){

               K_x = (Hamiltonian::ga_xy(d,0) + Hamiltonian::ga_xy(a_,0))%L;
               K_y = (Hamiltonian::ga_xy(d,1) + Hamiltonian::ga_xy(a_,1))%L;

               (*this)(B,i,j) -= norm * (2.0 * Z + 1.0) * _6j[S][Z] * phm(Z,K_x,K_y,d,a_,b,c_);

               K_x = (Hamiltonian::ga_xy(d,0) + Hamiltonian::ga_xy(b_,0))%L;
               K_y = (Hamiltonian::ga_xy(d,1) + Hamiltonian::ga_xy(b_,1))%L;

               (*this)(B,i,j) -=  norm * (2.0 * Z + 1.0) * _6j[S][Z] * sign * phm(Z,K_x,K_y,d,b_,a,c_);

               K_x = (Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(a_,0))%L;
               K_y = (Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(a_,1))%L;

               (*this)(B,i,j) -=  norm * (2.0 * Z + 1.0) * _6j[S][Z] * sign * phm(Z,K_x,K_y,c,a_,b,d_);

               K_x = (Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(b_,0))%L;
               K_y = (Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(b_,1))%L;

               (*this)(B,i,j) -=  norm * (2.0 * Z + 1.0) * _6j[S][Z] * phm(Z,K_x,K_y,c,b_,a,d_);

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * The bar function that maps a PPHM object onto a TPM object by tracing away the last pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void TPM::bar(const PPHM &pphm){

   int a,b,c,d;
   int Z;
   int K_x,K_y;

   double ward;

#pragma omp parallel for private(a,b,c,d,Z,K_x,K_y,ward)
   for(int B = 0;B < gnr();++B){//loop over the tp blocks

      Z = block_char[B][0];//spin of the TPM - block

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            for(int S = 0;S < 2;++S){//loop over three particle spin: 1/2 and 3/2

               ward = (2.0*(S + 0.5) + 1.0)/(2.0*Z + 1.0);

               for(int e = 0;e < M/2;++e){

                  K_x = (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(e,0))%L;
                  K_y = (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(e,1))%L;

                  (*this)(B,i,j) += ward * pphm(S,K_x,K_y,Z,a,b,e,Z,c,d,e);

               }

            }

         }
      }

   }

   this->symmetrize();

}

/* vim: set ts=3 sw=3 expandtab :*/
