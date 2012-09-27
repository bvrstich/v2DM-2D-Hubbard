/**
 * @mainpage 
 * This is an implementation of a boundary point method
 * for optimizing the second order density matrix using the P Q G T1 and T2 N-representability conditions for the special
 * case of the 1 dimensional Hubbard model with periodic boundary condition, which means that all the symmetries for this case have been 
 * implemented.
 * The method used is a path following algorithm with predictor corrector steps.
 * At compile time you can decide which condtions will be active compile with make PQ, PQG, PQGT1, PQGT2 or PQGT=(for all conditions).
 * @author Brecht Verstichel
 * @date 10-05-2010
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

using std::cout;
using std::endl;
using std::ofstream;

#include "include.h"

/**
 * 
 * In the main the actual program is run.\n 
 * a boundary point method
 */
int main(int argc,char *argv[]){

   //initialize the random nr generator
   srand(time(NULL));

   cout.precision(10);

   int L = atoi(argv[1]);//dimension of the lattice, nr of sites
   int N = atoi(argv[2]);//nr of particles

   double U = atof(argv[3]);//onsite repulsion

   Tools::init(L,N);
   
   Hamiltonian::init();

   TPM::init();

#ifdef __G_CON
   PHM::init();
#endif

#ifdef __T1_CON
   DPM::init();
#endif

#ifdef __T2_CON
   PPHM::init();
#endif

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

   u_0.gI().unit();

   u_0.fill();

   X = 0.0;
   Z = 0.0;

   //what does this do?
   double sigma = 1.0;

   double tolerance = 1.0e-7;

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   double mazzy = 1.6;

   int iter;
   int max_iter = 1;

   int tot_iter;

   while(D_conv > tolerance || P_conv > tolerance || fabs(convergence) > tolerance){

      D_conv = 1.0;

      iter = 0;

      while(D_conv > tolerance && iter <= max_iter){

         tot_iter++;

         ++iter;

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

         //check infeasibility of the dual problem:
         TPM v;

         v.collaps(1,V);

         v -= ham;

         D_conv = sqrt(v.ddot(v));

     }

      //update primal:
      X = V;

      //check primal feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= Z;

      P_conv = sqrt(W.ddot(W));

      convergence = Z.gI().ddot(ham) + X.ddot(u_0);

      cout << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << "\t" << Z.gI().ddot(ham_copy) << endl;

      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;

   }

   cout << endl;
   cout << "Energy: " << ham_copy.ddot(Z.gI()) << endl;
   cout << "pd gap: " << Z.ddot(X) << endl;
   cout << "dual conv: " << D_conv << endl;
   cout << "primal conv: " << P_conv << endl;

   cout << endl;
   cout << tot_iter << endl;
   cout << endl;

   SPM spm;
   spm.bar(1.0/(N - 1.0),Z.gI());

   cout << spm;

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

   Tools::clear();

   return 0;

}
