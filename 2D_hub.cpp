/**
 * @mainpage 
 * This is an implementation of a primal dual interior point method
 * for optimizing the second order density matrix using the P Q G T1 and T2 N-representability conditions for the special
 * case of the 2 dimensional Hubbard model with periodic boundary condition, which means that all the symmetries for this case have been 
 * implemented.
 * The method used is a path following algorithm with predictor corrector steps.
 * At compile time you can decide which condtions will be active compile with make PQ, PQG, PQGT1, PQGT2 or PQGT=(for all conditions).
 * @author Brecht Verstichel
 * @date 02-03-2011
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
 * Part 1: An easy initial point is taken and then centered to the required precision (flag == 0)\n
 * Part 2: When the primal dual point is sufficiently centered steps are taken to reduce the
 * primal dual gap and take a large step in that direction (predictor) (flag == 1)\n
 * After each step a correcting step (flag == 2) is taken that brings the primal dual point closer to
 * the central path.\n
 * Part 3: When the primal dual gap is smaller that the required accuracy exit the while. (flag == 3)\n
 * For more information on the actual method, see primal_dual.pdf
 */
int main(int argc,char *argv[]){

   //initialize the random nr generator
   srand(time(NULL));

   cout.precision(10);

   int L = atof(argv[1]);
   int N = atof(argv[2]);

   double U = atof(argv[3]);

   Hamiltonian::init(L);
   TPM::init(L,N);
   SPM::init(L,N);

   TPM tpm;
   tpm.fill_Random();

   tpm.out_sp("../spin_pd-bright/tpm.in");

   TPM Q;
   Q.Q(1,tpm);

   BlockVector<TPM> v(Q);
   cout << v;

/*
   SUP S(M,N);
   S.init_S();

   SUP Z(M,N);
   Z.init_Z(10000.0,ham,S);

   int dim = Z.gdim();

   //eerste primal dual gap:
   double pd_gap = S.ddot(Z);
   double energy = (S.tpm(0)).ddot(ham);

   double center_dev = S.center_dev(Z);

   //eerst centering
   double gamma = 1.0;

   double tolerance = 1.0e-10;

   //flag == 0 : initiele centering run (tot op tolerance)
   //flag == 1 : doe een stap met gamma = 0
   //flag == 2 : doe een stap met gamma = 1
   //flag == 3 : game over man
   int flag = 0;

   double a;//stapgrootte

   while(flag != 3){

      cout << (S.tpm(0)).trace() << "\t" << pd_gap << "\t" << center_dev << "\t" << energy << "\t" << S.tpm(0).spin() << "\t";

      //matrix D aanmaken voor de hessiaan van het duale stelsel
      SUP D(M,N);
      D.D(S,Z);

      //D inverteren voor de hessiaan van het primale stelsel
      SUP D_inv(D);
      D_inv.invert();

      //rechterlid maken van stelsel dat moet worden opgelost:
      SUP B(S);

      //invert B
      B.invert();

      //schalen met 
      B.dscal(gamma*pd_gap/dim);

      B -= Z;

      //collaps B onto b to construct the right hand side of the primal Newton equation
      TPM b(M,N);

      b.collaps(1,B);

      //dit wordt de stap:
      TPM delta(M,N);

      //los het stelsel op, geeft aantal iteraties nodig terug:
      cout << delta.solve(b,D_inv) << "\t";

      //nog updaten van S en Z
      SUP DS(M,N);

      DS.fill(delta);

      //DZ is B - D^{-1}*DS*D^{-1}
      SUP DZ(B);

      //eerst D^{-1}*DS*D^{-1} in DZ stoppen
      B.L_map(D_inv,DS);

      DZ -= B;

      //voor de zekerheid nog projecteren op juiste subruimte:
      DZ.proj_C();

      //met deze 'ansatz' het Z stelsel proberen op te lossen
      //eerst rechterlid B maken
      B = Z;

      B.invert();

      B.dscal(gamma*pd_gap/dim);

      B -= S;

      B.proj_C();

      //los het stelsel op, geeft aantal duale iteraties nodig terug:
      cout << DZ.solve(B,D) << endl;

      //welke stapgrootte moet ik nemen?
      if(flag == 0 || flag == 2){//voor centering

         S += DS;
         Z += DZ;

      }
      else{

         //zoek de ideale afstand (geef ook een waarde mee voor de maximale afwijking van het centraal pad):
         a = DS.line_search(DZ,S,Z,2.0);

         S.daxpy(a,DS);
         Z.daxpy(a,DZ);

      }

      //update van enkele belangrijke variabelen
      pd_gap = S.ddot(Z);
      energy = (S.tpm(0)).ddot(ham);
      center_dev = S.center_dev(Z);

      //keuze voor volgende iteratie:
      if(flag == 0){

         //als hij voldoende gecenterd is, exit.
         if(center_dev < tolerance){

            flag = 1;
            gamma = 0.0;

         }

      }
      else if(flag == 1){

         if(pd_gap < tolerance)//exit when converged
            flag = 3;
         else{//center when not convergence

            flag = 2;
            gamma = 1.0;

         }

      }
      else{//flag == 2: dus na een centering stap

         if(pd_gap < tolerance)//exit when converged
            flag = 3;
         else{//take another step downwards when not converged

            flag = 1;
            gamma = 0;

         }

      }

   }

   cout << endl;
   cout << "FINAL RESULT " << endl;
   cout << endl;
   cout << "E_0 = " << energy << " with accuracy of " << pd_gap << " and a deviation from centrality of " << center_dev << endl;
   cout << endl;
   cout << "<S^2>\t=\t" << S.tpm(0).spin() << endl;

   //print density matrix to file
   //(S.tpm(0)).out("rdm.out");

   TPM kinetic(M,N);
   kinetic.hubbard(0.0);

   cout << endl;
   cout << "Square root of the trace pair:\t" << S.tpm(0).trace_pair() << endl;
   cout << endl;
   cout << "expectation value of the kinetic energy\t" << S.tpm(0).ddot(kinetic) << endl;
*/
   TPM::clear();
   Hamiltonian::clear();

   return 0;

}
