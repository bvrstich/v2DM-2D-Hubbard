#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * standard constructor\n
 * Allocates two TPM matrices and optionally a PHM, DPM or PPHM matrix.
 * @param M number of sp orbitals
 * @param N number of particles
 */
SUP::SUP(int M,int N){

   this->M = M;
   this->N = N;
   this->n_tp = M*(M - 1)/2;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   this->dim = 2*n_tp;

}

/**
 * standard constructor\n
 * Allocates two TPM matrices and optionally a PHM, DPM or PPHM matrix, then copies the content of
 * input SUP SZ_c into it.
 * @param SZ_c input SUP
 */
SUP::SUP(const SUP &SZ_c){

   this->M = SZ_c.M;
   this->N = SZ_c.N;
   this->n_tp = SZ_c.n_tp;
   this->dim = 2*n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

}

/**
 * Destructor
 */
SUP::~SUP(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

}

/**
 * Overload += operator
 * @param SZ_pl The SUP matrix that has to be added to this
 */
SUP &SUP::operator+=(const SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

   return *this;

}

/**
 * Overload -= operator
 * @param SZ_pl The SUP that will be deducted from this
 */
SUP &SUP::operator-=(const SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

   return *this;

}

/**
 * Overload equality operator, copy SZ_c into this
 * @param SZ_c SUP_PQ to be copied into this
 */
SUP &SUP::operator=(const SUP &SZ_c){

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP &SUP::operator=(double &a){

   (*SZ_tp[0]) = a;
   (*SZ_tp[1]) = a;

   return *this;

}

/**
 * @param i which block you want to have the pointer to.
 * @return pointer to the individual TPM blocks: SZ_tp[i]
 */
TPM &SUP::tpm(int i){

   return *SZ_tp[i];

}

/**
 * const version
 * @param i which block you want to have the pointer to.
 * @return pointer to the individual TPM blocks: SZ_tp[i]
 */
const TPM &SUP::tpm(int i) const{

   return *SZ_tp[i];

}

/**
 * Initialization of the SUP matrix S, is just u^0: see primal_dual.pdf for more information
 */
void SUP::init_S(){

   (*SZ_tp[0]).unit();

   this->fill();

}

ostream &operator<<(ostream &output,const SUP &SZ_p){

   output << (*SZ_p.SZ_tp[0]) << std::endl;
   output << (*SZ_p.SZ_tp[1]);

   return output;

}

/**
 * Fill the SUP matrix with random elements, watch out, not necessarily positive definite
 */
void SUP::fill_Random(){

   SZ_tp[0]->fill_Random();
   SZ_tp[1]->fill_Random();

}

/**
 * Initialisation for dual SUP matrix Z, see primal_dual.pdf for info.
 */
void SUP::init_Z(double alpha,const TPM &ham,const SUP &u_0){

   this->fill_Random();

   this->proj_C(ham);

   //nog een eenheidsmatrix maal constante bijtellen zodat Z positief definiet is:
   this->daxpy(alpha,u_0); 

}

/**
 * @return number of particles
 */
int SUP::gN() const {

   return N;

}

/**
 * @return dimension of sp space
 */
int SUP::gM() const{

   return M;

}

/**
 * @return dimension of tp space
 */
int SUP::gn_tp() const{

   return n_tp;

}

/**
 * @return total dimension of SUP (carrier) space
 */
int SUP::gdim() const{

   return dim;

}

/**
 * @param SZ_i input SUP_PQ SZ_i
 * @return inproduct between this and input matrix SZ_i, defined as Tr(this SZ_i)
 */
double SUP::ddot(const SUP &SZ_i) const{

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 * Makes use of cholesky decomposition, so only positive definite matrices can be used as input!
 */
void SUP::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

}

/**
 * Orthogonal projection of a general SUP matrix diag[ M M_Q ( M_G M_T1 M_T2 ) ] onto U space: diag[ M_u Q(M_u) ( G(M_u) T1(M_u) T2(M_u) ) ]
 * for more information see primal_dual.pdf
 */
void SUP::proj_U(){
  
   //eerst M_Gamma + Q(M_Q) + ( G(M_G) + T1(M_T1) + T2(M_T2) ) in O stoppen
   TPM O(M,N);

   O.collaps(1,*this);

   //dan de inverse overlapmatrix hierop laten inwerken en in this[0] stoppen
   SZ_tp[0]->S(-1,O);

   //fill up the rest with the right maps
   this->fill();

}

/**
 * Project the general SUP matrix (*this) orthogonally onto the linear space for which\n\n
 * Tr(Z u^i) = h^i      with h^i = Tr(tpm f^i)\n\n
 * is valid.
 * @param tpm input TPM (mostly the hamiltonian of the problem)
 */
void SUP::proj_C(const TPM &tpm){

   TPM hulp(M,N);

   hulp.collaps(1,*this);

   hulp -= tpm;

   //Z_res is the orthogonal piece of this that will be deducted,
   //so the piece of this in the U-space - ham
   SUP Z_res(M,N);

   //apply iverse S to it and put it in Z_res.tpm(0)
   (Z_res.tpm(0)).S(-1,hulp);

   //and fill it up Johnny
   Z_res.fill();

   *this -= Z_res;

}

/**
 * Construct the D matrix and put it in this, D is the matrix matrix of the hessian, see primal_dual.pdf for more information
 * @param S The primal SUP matrix S
 * @param Z The dual SUP matrix Z
 */
void SUP::D(const SUP &S,const SUP &Z){

   //positieve vierkantswortel uit Z
   SUP Z_copy(Z);

   Z_copy.sqrt(1);

   //links en rechts vermenigvuldigen met wortel Z
   SUP hulp(M,N);

   hulp.L_map(Z_copy,S);

   hulp.sqrt(1);

   //negatieve vierkantswortel uit Z
   Z_copy = Z;

   Z_copy.sqrt(-1);

   //en links en rechts hulp vermenigvuldigen met die wortel, en in this steken:
   this->L_map(Z_copy,hulp);

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

}

/**
 * Multiply symmetric SUP blockmatrix object left en right with symmetric SUP blockmatrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map SUP that will be multiplied to the left en to the right of matrix object
 * @param object central SUP
 */
void SUP::L_map(const SUP &map,const SUP &object){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->L_map(map.tpm(i),object.tpm(i));

}

/**
 * add the SUP SZ_p times the constant alpha to this
 * @param alpha the constant to multiply the SZ_p with
 * @param SZ_p the SUP to be multiplied by alpha and added to (*this)
 */
void SUP::daxpy(double alpha,const SUP &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,SZ_p.tpm(i));

}

/**
 * Orthogonal projection of a general SUP matrix [ M M_Q ( M_G M_T1 M_T2 ) ] onto the orthogonal complement of the U space (C space)
 * See primal_dual.pdf for more information
 */
void SUP::proj_C(){

   SUP Z_copy(*this);

   //projecteer op de U ruimte
   Z_copy.proj_U();

   //en het orthogonaal complement nemen:
   *this -= Z_copy;

}

/**
 * General matrixproduct between two SUP matrices, act with Matrix::mprod on every block
 * 
 * @param A left hand matrix
 * @param B right hand matrix
 * @return The product AB
 */
SUP &SUP::mprod(const SUP &A,const SUP &B){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->mprod(A.tpm(i),B.tpm(i));

   return *this;

}

/**
 * Fill the SUP matrix (*this) with a TPM matrix like: this = diag[tpm  Q(tpm)  ( G(tpm) T1(tpm) T2(tpm) ) ]
 * @param tpm input TPM
 */
void SUP::fill(const TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(1,tpm);

}

/**
 * fill the SUP matrix with the TPM matrix stored in the first block:\n\n
 * this = diag[this->tpm(0) Q(this->tpm(0)) ( G(this->tpm(0)) T1(this->tpm(0)) T2(this-tpm(0)) ) ]
 */
void SUP::fill(){

   SZ_tp[1]->Q(1,*SZ_tp[0]);

}

/**
 * Implementation of the linear conjugate gradient algorithm for the solution of the dual Newton system\n\n
 * H(*this) = B in which H is the dual hessian map
 * @param B right hand side of the equation
 * @param D SUP matrix that defines the structure of the hessian map (the metric) (inverse of the primal Newton equation hessian)
 * @return return the number of iteration required to converge
 */
int SUP::solve(SUP &B,const SUP &D){

   SUP HB(M,N);
   HB.H(*this,D);

   B -= HB;

   //de r initialiseren op B - H DZ
   SUP r(B);

   double rr = r.ddot(r);
   double rr_old,ward;

   int cg_iter = 0;

   while(rr > 1.0e-10){

      ++cg_iter;

      HB.H(B,D);

      ward = rr/B.ddot(HB);

      //delta += ward*b
      this->daxpy(ward,B);

      //r -= ward*HB
      r.daxpy(-ward,HB);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      B.dscal(rr/rr_old);

      //dan r er bijtellen
      B += r;

   }
   
   return cg_iter;

}

/**
 * The dual hessian map:\n\n
 * HB = DBD (dus SUP::L_map), projected onto C-space (SUP::proj_C)
 * @param B SUP matrix onto which the hessian works.
 * @param D SUP matrix that defines the structure of the map (metric)
 */
void SUP::H(const SUP &B,const SUP &D){

   this->L_map(D,B);

   this->proj_C();

}

/**
 * @return Deviation from the central path measured through the logarithmic potential, it's a measure for
 * the deviation of the product of the primal with the dual matrix (SZ) from the unit matrix.\n
 * Usage of the function: S.center_dev(Z) gives returns the deviation.\n\n
 * (*this) = S = primal matrix of the problem
 * @param Z = dual matrix of the problem
 */
double SUP::center_dev(const SUP &Z) const{

   SUP sqrt_S(*this);

   sqrt_S.sqrt(1);

   SUP SZ(M,N);
   SZ.L_map(sqrt_S,Z);

   EIG eig(SZ);

   return eig.center_dev();

}

/**
 * Line search function that checks how large a step you can take in a given primal dual predictor direction (DS,DZ), starting from 
 * the current primal dual point (S,Z), before deviating beyond max_dev from the central path.\n\n
 * (*this) = DS --> primal search direction
 * @param DZ dual search direction
 * @param S Current primal point
 * @param Z Current dual point
 * @param max_dev number (double) input by which you can tell the function how far you want to deviate from the central path after the step.
 */
double SUP::line_search(const SUP &DZ,const SUP &S,const SUP &Z,double max_dev) const{

   //eerst de huidige deviatie van het centraal pad nemen:
   double center_dev = S.center_dev(Z);

   //eigenwaarden zoeken van S^{-1/2} DS S^{-1/2} en Z^{-1/2} DZ Z^{-1/2}

   //kopieer S in de zogeheten wortel:
   SUP wortel(S);

   //maak negatieve vierkantswortel uit S
   wortel.sqrt(-1);

   //de L_map
   SUP hulp(M,N);
   hulp.L_map(wortel,*this);

   //eigenwaarden in eigen_S stoppen
   EIG eigen_S(hulp);

   //nu idem voor Z
   wortel = Z;

   wortel.sqrt(-1);

   hulp.L_map(wortel,DZ);

   EIG eigen_Z(hulp);

   //nog c_S en c_Z uitrekenen:
   double pd_gap = S.ddot(Z);

   //c_S = Tr (DS Z)/Tr (SZ)
   double c_S = this->ddot(Z)/pd_gap;

   //c_Z = Tr (S DZ)/Tr (SZ)
   double c_Z = S.ddot(DZ)/pd_gap;

   //waar zitten de singulariteiten: tot waar mag ik zoeken?
   double a_max = -1.0/eigen_S.min();
   double b_max = -1.0/eigen_Z.min();

   //a_max is de waarde tot waar ik zal zoeken:
   if(b_max < a_max)
      a_max = b_max;

   double a = 0.0;
   double b = a_max;

   double c = (a + b)/2.0;

   //bissectiemethode om stapgrootte te bepalen:
   while(b - a > 1.0e-5){

      c = (a + b)/2.0;

      if( (center_dev + eigen_S.centerpot(c,eigen_Z,c_S,c_Z) - max_dev) < 0.0 )
         a = c;
      else
         b = c;

   }

   return c;

}

/*
 * Output the SUP matrix to a file, associated with the ofstream object output
 * @param output The ofstream object you will dump the SUP in.
 */
void SUP::out(ofstream &output) const{

   output.precision(15);

   //P
   for(int B = 0;B < this->tpm(0).gnr();++B){

      for(int i = 0;i < this->tpm(0).gdim(B);++i)
         for(int j = 0;j < this->tpm(0).gdim(B);++j)
            output << i << "\t" << j << "\t" << this->tpm(0)(B,i,j) << std::endl;

   }

   //Q
   for(int B = 0;B < this->tpm(1).gnr();++B){

      for(int i = 0;i < this->tpm(1).gdim(B);++i)
         for(int j = 0;j < this->tpm(1).gdim(B);++j)
            output << i << "\t" << j << "\t" << this->tpm(1)(B,i,j) << std::endl;

   }

}

/*
 * Input a SUP matrix from a file
 * @param input ifstream object association with input file.
 */
void SUP::in(ifstream &input){

   int I,J;

   //first P
   for(int B = 0;B < (this->tpm(0)).gnr();++B){

      for(int i = 0;i < (this->tpm(0)).gdim(B);++i)
         for(int j = 0;j < (this->tpm(0)).gdim(B);++j)
            input >> I >> J >> (this->tpm(0))(B,i,j);

   }

   //then Q
   for(int B = 0;B < (this->tpm(1)).gnr();++B){

      for(int i = 0;i < (this->tpm(1)).gdim(B);++i)
         for(int j = 0;j < (this->tpm(1)).gdim(B);++j)
            input >> I >> J >> (this->tpm(1))(B,i,j);

   }

}
