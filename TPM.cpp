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

double TPM::Sa = 1;
double TPM::Sc = 0;

/**
 * static function that initializes the static variables
 */
void TPM::init(){

   //allocate stuff
   t2s = new vector< vector<int> > [Tools::gM()];

   s2t = new int ** [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B){

      s2t[B] = new int * [Tools::gL()*Tools::gL()];

      for(int a = 0;a < Tools::gL()*Tools::gL();++a)
         s2t[B][a] = new int [Tools::gL()*Tools::gL()];

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

   //tp index
   int t;

   //loop over the K_x K_y blocks
   for(int K_x = 0;K_x < Tools::gL();++K_x)
      for(int K_y = 0;K_y < Tools::gL();++K_y){

         t = 0;

         //S = 0
         block_char[block][0] = 0;
         block_char[block][1] = K_x;
         block_char[block][2] = K_y;

         char_block[0][K_x][K_y] = block;

         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a;b < Tools::gL()*Tools::gL();++b){

               if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % Tools::gL() == K_x ) 

                     && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % Tools::gL() == K_y ) ){

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
         block_char[Tools::gL()*Tools::gL() + block][0] = 1;
         block_char[Tools::gL()*Tools::gL() + block][1] = K_x;
         block_char[Tools::gL()*Tools::gL() + block][2] = K_y;

         char_block[1][K_x][K_y] = Tools::gL()*Tools::gL() + block;

         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b){

               if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % Tools::gL() == K_x ) 

                     && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % Tools::gL() == K_y ) ){

                  v[0] = a;
                  v[1] = b;

                  t2s[Tools::gL()*Tools::gL() + block].push_back(v);

                  s2t[Tools::gL()*Tools::gL() + block][a][b] = t;
                  s2t[Tools::gL()*Tools::gL() + block][b][a] = t;

                  ++t;

               }

            }

         ++block;

      }

   init_overlap();

}

/**
 * initialize the overlapmatrix parameters
 */
void TPM::init_overlap(){

#ifdef __Q_CON
   Sa += 1.0;
   Sc += (2.0*Tools::gN() - Tools::gM())/((Tools::gN() - 1.0)*(Tools::gN() - 1.0));
#endif

#ifdef __G_CON
   Sa += 4.0;
   Sc += (2.0*Tools::gN() - Tools::gM() - 2.0)/((Tools::gN() - 1.0)*(Tools::gN() - 1.0));
#endif

#ifdef __T1_CON
   Sa += Tools::gM() - 4.0;
   Sc -= (Tools::gM()*Tools::gM() + 2.0*Tools::gN()*Tools::gN() - 4.0*Tools::gM()*Tools::gN() - Tools::gM() + 8.0*Tools::gN() - 4.0)
   
      /( 2.0*(Tools::gN() - 1.0)*(Tools::gN() - 1.0) );
#endif

#ifdef __T2_CON
   Sa += 5.0*Tools::gM() - 4.0;
   Sc += (2.0*Tools::gN()*Tools::gN() + (Tools::gM() - 2.0)*(4.0*Tools::gN() - 3.0) - Tools::gM()*Tools::gM() - 2.0)

      /(2.0*(Tools::gN() - 1.0)*(Tools::gN() - 1.0));
#endif

}

/**
 * static function that deallocates the static variables
 */
void TPM::clear(){

   delete [] t2s;

   for(int B = 0;B < Tools::gM();++B){

      for(int a = 0;a < Tools::gL()*Tools::gL();++a)
         delete [] s2t[B][a];

      delete [] s2t[B];

      delete [] block_char[B];

   }

   delete [] s2t;

   delete [] block_char;

   for(int S = 0;S < 2;++S){

      for(int x = 0;x < Tools::gL();++x)
         delete [] char_block[S][x];

      delete [] char_block[S];

   }

   delete [] char_block;

}

/**
 * standard constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice: 
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 */
TPM::TPM() : BlockMatrix(Tools::gM()) {

   //set the dimension of the blocks

   for(int B = 0;B < Tools::gL()*Tools::gL();++B)//S = 0
      setMatrixDim(B,t2s[B].size(),1);

   for(int B = Tools::gL()*Tools::gL();B < 2*Tools::gL()*Tools::gL();++B)//S = 1
      setMatrixDim(B,t2s[B].size(),3);

}

/**
 * copy constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice:
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 * @param tpm_c The TPM object to be copied into (*this)
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor
 */
TPM::~TPM(){ }

ostream &operator<<(ostream &output,const TPM &tpm_p){

   int S,K_x,K_y;

   for(int B = 0;B < tpm_p.gnr();++B){

      S = tpm_p.block_char[B][0];
      K_x = tpm_p.block_char[B][1];
      K_y = tpm_p.block_char[B][2];

      output << "S =\t" << S << "\tK_x =\t" << K_x << "\tK_y =\t" << K_y <<

         "\tdimension =\t" << tpm_p.gdim(B) << "\tdegeneracy =\t" << tpm_p.gdeg(B) << std::endl;

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
 * @param S The tp spin quantumnumber
 * @param a first sp index that forms the tp row index i of block B(S,K_x,K_y), together with b
 * @param b second sp index that forms the tp row index i of block B(S,K_x,K_y), together with a
 * @param c first sp index that forms the tp column index j of block B(S,K_x,K_y), together with d
 * @param d second sp index that forms the tp column index j of block B(S,K_x,K_y), together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int a,int b,int c,int d) const{

   int K_x = ( Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) )%Tools::gL();
   int K_y = ( Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) )%Tools::gL();

   //check if momentum is ok
   if( K_x != ( Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(d,0) )%Tools::gL() )
      return 0;

   if( K_y !=  ( Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(d,1) )%Tools::gL() )
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
 * construct the spinsymmetrical hubbard hamiltonian in momentum space with on site repulsion U
 * @param U onsite repulsion term
 */
void TPM::hubbard(double U){

   int a,b,c,d;//sp indices

   double ward = 1.0/(Tools::gN() - 1.0);

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

               (*this)(B,i,i) = -2.0 * ward * ( cos( 2.0 * Hamiltonian::ga_xy(a,0) * 3.141592653589793238462 / (double) Tools::gL())  

                     + cos( 2.0 * Hamiltonian::ga_xy(a,1) * 3.141592653589793238462 / (double) Tools::gL()) 

                     + cos( 2.0 * Hamiltonian::ga_xy(b,0) * 3.141592653589793238462 / (double) Tools::gL()) 

                     + cos( 2.0 * Hamiltonian::ga_xy(b,1) * 3.141592653589793238462 / (double) Tools::gL()) );

            }

            //on-site repulsion
            if(B < Tools::gL()*Tools::gL()){//only spin zero

               double hard = 2.0*U / (double) (Tools::gL()*Tools::gL());

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
   double b = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double c = 1.0/(Tools::gN() - 1.0);

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

      B = (B*A + B*C*Tools::gM() - 2.0*C*C)/( A * (C*(Tools::gM() - 2.0) -  A) * ( A + B*Tools::gM()*(Tools::gM() - 1.0) - 2.0*C*(Tools::gM() - 1.0) ) );
      C = C/(A*(C*(Tools::gM() - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm;
   spm.bar(C,tpm_d);

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

   double ward = Tools::gN()*(Tools::gN() - 1.0)/(Tools::gM()*(Tools::gM() - 1.0));

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

   double ward = (2.0 * this->trace())/(Tools::gM()*(Tools::gM() - 1));

   this->min_unit(ward);

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
 * @param sup input SUP
 */
void TPM::collaps(int option,const SUP &sup){

   *this = sup.gI();

#ifdef __Q_CON
   TPM hulp;

   hulp.Q(1,sup.gQ());

   *this += hulp;
#endif

#ifdef __G_CON
   hulp.G(sup.gG());

   *this += hulp;
#endif

#ifdef __T1_CON
   hulp.T(sup.gT1());

   *this += hulp;
#endif

#ifdef __T2_CON
   hulp.T(sup.gT2());

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

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      if(S == 0){

         for(int i = 0;i < gdim(B);++i)
            ward += -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) * (*this)(B,i,i);

      }
      else{

         for(int i = 0;i < this->gdim(B);++i)
            ward += 3.0 * ( -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0 ) * (*this)(B,i,i);

      }

   }

   return ward;

}

/**
 * Fill (*this) with the S^2 operator
 */
void TPM::set_S_2(){

   *this = 0.0;

   //S = 0 blocks
   for(int B = 0;B < Tools::gL()*Tools::gL();++B)
      for(int i = 0;i < gdim(B);++i)
         (*this)(B,i,i) = -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0);

   //S = 1 blocks
   for(int B = Tools::gL()*Tools::gL();B < Tools::gM();++B)
      for(int i = 0;i < gdim(B);++i)
         (*this)(B,i,i) = -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0;

}

/**
 * The G down map, maps a PHM object onto a TPM object using the G map.
 * @param phm input PHM
 */
void TPM::G(const PHM &phm){

   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),phm);

   int a,b,c,d;

   //the conjugated indices
   int a_,b_,c_,d_;

   int S;

   int sign;

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
            for(int Z = 0;Z < 2;++Z)
               (*this)(B,i,j) -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * ( phm(Z,a,d_,c,b_) + phm(Z,b,c_ ,d, a_) 

                     + sign * ( phm(Z,b,d_ ,c, a_ ) + phm(Z,a,c_,d,b_) ) );

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
 * The T1-down map that maps a DPM on TPM. This is just a Q-like map using the TPM::bar (dpm) as input.
 * @param dpm the input DPM matrix
 */
void TPM::T(const DPM &dpm){

   TPM tpm;
   tpm.bar(dpm);

   double a = 1;
   double b = 1.0/(3.0*Tools::gN()*(Tools::gN() - 1.0));
   double c = 0.5/(Tools::gN() - 1.0);

   this->Q(1,a,b,c,tpm);

}

/**
 * Construct a spincoupled, translationally invariant TPM matrix out of a spincoupled, translationally invariant DPM matrix.
 * For the definition and derivation see symmetry.pdf
 * @param dpm input DPM
 */
void TPM::bar(const DPM &dpm){

   int a,b,c,d;

   double ward;

   //first the S = 0 part, easiest:
   for(int B = 0;B < Tools::gL()*Tools::gL();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            //only total S = 1/2 can remain because cannot couple to S = 3/2 with intermediate S = 0
            for(int e = 0;e < Tools::gL()*Tools::gL();++e)
               (*this)(B,i,j) += 2.0 * dpm(0,0,a,b,e,0,c,d,e);

         }
      }

   }

   //then the S = 1 part:
   for(int B = Tools::gL()*Tools::gL();B < Tools::gM();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            for(int Z = 0;Z < 2;++Z){//loop over the dpm blocks: S = 1/2 and 3/2 = Z + 1/2

               ward = 0.0;

               for(int e = 0;e < Tools::gL()*Tools::gL();++e)
                  ward += dpm(Z,1,a,b,e,1,c,d,e);

               ward *= (2 * (Z + 0.5) + 1.0)/3.0;

               (*this)(B,i,j) += ward;

            }

         }
      }

   }

   this->symmetrize();

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
   spm.bar(0.5/(Tools::gN() - 1.0),pphm);

   int a,b,c,d;
   int a_,b_,c_,d_;
   int sign;

   double norm_ab,norm_cd;

   int S;

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
            norm_ab = 1.0;
            norm_cd = 1.0;

            if(S == 0){

               if(a == b)
                  norm_ab /= std::sqrt(2.0);

               if(c == d)
                  norm_cd /= std::sqrt(2.0);

            }

            //first the tp part
            (*this)(B,i,j) = tpm(B,i,j) - std::sqrt(1.0/(S + 0.5)) * ( norm_cd * (pphm.w(S,a,b,c_,d) + sign * pphm.w(S,a,b,d_,c))
            
                  + norm_ab * (pphm.w(S,c,d,a_,b) + sign * pphm.w(S,c,d,b_,a)) );

            //sp part is diagonal for translational invariance
            if(i == j)
               (*this)(B,i,j) += spm[a_] + spm[b_] + ( pphm.sp(a) + pphm.sp(b) )/(Tools::gN() - 1.0);

            for(int Z = 0;Z < 2;++Z){

               (*this)(B,i,j) -= norm_ab * norm_cd * (2.0 * Z + 1.0) * Tools::g6j(0,0,S,Z) * ( phm(Z,d,a_,b,c_) + sign * phm(Z,d,b_,a,c_) 
               
                     + sign * phm(Z,c,a_,b,d_) +  phm(Z,c,b_,a,d_) );

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

   double ward;

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

               for(int e = 0;e < Tools::gL()*Tools::gL();++e)
                  (*this)(B,i,j) += ward * pphm.pph(S,Z,a,b,e,Z,c,d,e);

            }

         }

      }

   }

   this->symmetrize();

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
 * Primal hessian map:\n\n
 * Hb = D_1 b D_1 + D_2 Q(b) D_2 + D_3 G(b) D_3 + D_4 T1(b) D_4 + D_5 T2(b) D5 \n\n
 * with D_1, D_2, D_3, D_4 and D_5 the P, Q, G, T1 and T2 blocks of the SUP D. 
 * @param b TPM domain matrix, hessian will act on it and the image will be put in this
 * @param D SUP matrix that defines the structure of the hessian map. (see primal-dual.pdf for more info)
 */
void TPM::H(const TPM &b,const SUP &D){

   this->L_map(D.gI(),b);

#ifdef __Q_CON
   //maak Q(b)
   TPM Qb;
   Qb.Q(1,b);

   TPM hulp;

   hulp.L_map(D.gQ(),Qb);

   Qb.Q(1,hulp);

   *this += Qb;
#endif

#ifdef __G_CON
   //maak G(b)
   PHM Gb;
   Gb.G(b);

   PHM hulpje;

   hulpje.L_map(D.gG(),Gb);

   hulp.G(hulpje);

   *this += hulp;
#endif

#ifdef __T1_CON
   DPM T1b;
   T1b.T(b);

   DPM hulp_T1;

   hulp_T1.L_map(D.gT1(),T1b);

   hulp.T(hulp_T1);

   *this += hulp;
#endif

#ifdef __T2_CON
   PPHM T2b;
   T2b.T(b);

   PPHM hulp_T2;

   hulp_T2.L_map(D.gT2(),T2b);

   hulp.T(hulp_T2);

   *this += hulp;
#endif

   this->proj_Tr();

}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param option = 1 direct overlapmatrix-map is used , = -1 inverse overlapmatrix map is used
 * @param tpm_d the input TPM
 */
void TPM::S(int option,const TPM &tpm_d){

   this->Q(option,Sa,0.0,Sc,tpm_d);

}
