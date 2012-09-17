/**
 * @mainpage 
 * This is an implementation of a boundary point method to solve a semidefinite program:
 * we optimizing the second order density matrix using the P Q G T1 and T2 N-representability conditions. 
 * This is the same program but with the boundary point method.
 * At compile time you can decide which condtions will be active compile with make PQ, PQG, PQGT1, PQGT2 or PQGT=(for all conditions).
 * @author Brecht Verstichel, Ward Poelmans
 * @date 21-01-2011
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>

using std::cout;
using std::endl;
using std::ofstream;

#include "include.h"

/**
 * 
 * In the main the actual program is run.\n 
 */

int main(int argc,char **argv)
{
   cout.precision(10);

   // these are the default values
   int L = 3;//dim sp hilbert space
   int N = 9;//nr of particles

   double U = 1;//onsite interaction strength

   Tools::init(L,N);

   Hamiltonian::init();
   TPM::init();

#ifdef __G_CON
   PHM::init();
#endif

#ifdef __T1_CON
   DPM::init(L,N);
#endif

#ifdef __T2_CON
   PPHM::init(L,N);
#endif

   SPM::init(L,N);
   SUP::init(L,N);
   EIG::init(L,N);

   //hamiltoniaan
   TPM ham;
   ham.hubbard(U);

   TPM ham_copy(ham);

   //only traceless hamiltonian needed in program.
   ham.proj_Tr();

   //primal
   SUP X;

   //dual
   SUP Z;

   //Lagrange multiplier
   SUP V;

   //just dubya
   SUP W;

   SUP u_0;

   //little help
   TPM hulp;

   u_0.tpm(0).unit();

   u_0.fill();

   X = 0.0;
   Z = 0.0;

   //what does this do?
   double sigma = 1.0;

   double tolerance = 1.0e-7;

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   // mazziotti uses 1.6 for this
   double mazzy = 1.0;

   int iter_dual,iter_primal(0);
   int max_iter = 1;

   while(P_conv > tolerance || D_conv > tolerance || fabs(convergence) > tolerance){

      ++iter_primal;

      D_conv = 1.0;

      iter_dual = 0;

      while(D_conv > tolerance  && iter_dual <= max_iter)
      {

         ++iter_dual;

         //solve system
         SUP B(Z);

         B -= u_0;

         B.daxpy(mazzy/sigma,X);

         TPM b;

         b.collaps(1,B);

         b.daxpy(-mazzy/sigma,ham);

         hulp.S(-1,b);

         //hulp is the matrix containing the gamma_i's
         hulp.proj_Tr();

         //construct W
         W.fill(hulp);

         W += u_0;

         W.daxpy(-1.0/sigma,X);

         //update Z and V with eigenvalue decomposition:
         W.sep_pm(Z,V);

         V.dscal(-sigma);

         //check infeasibility of the primal problem:
         TPM v;

         v.collaps(1,V);

         v -= ham;

         D_conv = sqrt(v.ddot(v));

         //cout << "D\t\t\t" << D_conv << endl;

     }

      //update primal:
      X = V;

      //check dual feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= Z;

      P_conv = sqrt(W.ddot(W));

      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;

      convergence = Z.tpm(0).ddot(ham) + u_0.ddot(X);

      cout << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << "\t" << ham_copy.ddot(Z.tpm(0)) << endl;

   }

   cout << endl;
   cout << "Energy: " << ham_copy.ddot(Z.tpm(0)) << endl;
   cout << "pd gap: " << Z.ddot(X) << endl;
   cout << "dual conv: " << D_conv << endl;
   cout << "primal conv: " << P_conv << endl;

#ifdef __T2_CON
   PPHM::clear();
#endif

#ifdef __T1_CON
   DPM::clear();
#endif

#ifdef __G_CON
   PHM::clear();
#endif

   TPM::clear();
   Hamiltonian::clear();
   Tools::clear();

   return 0;

}

/* vim: set ts=3 sw=3 expandtab :*/
