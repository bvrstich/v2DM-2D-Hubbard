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
   int L = 6;//dim sp hilbert space
   int N = 36;//nr of particles

   double U = 4;//onsite interaction strength

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

   sTPM::init();

   ifstream in("/home/bright/bestanden/results/2D_hub/PQG/DM_out/6x6/N34/U8.dm");

   TPM tpm;

   for(int B = 0;B < tpm.gnr();++B)
      for(int i = 0;i < tpm.gdim(B);++i)
         for(int j = i;j < tpm.gdim(B);++j)
            in >> B >> i >> j >> tpm(B,i,j);

   tpm.symmetrize();

   Stripe stripe(tpm);

   cout << stripe.gN() << "\t" << stripe.gnpairs() << endl;

   Stripe ham;
   ham.hubbard1D(8.0);

   cout << ham.ddot(stripe) << endl;

   sTPM::clear();

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
