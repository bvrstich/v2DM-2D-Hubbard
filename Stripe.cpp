#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * constructor 
 * @param n dimension of the stripe
 */
Stripe::Stripe(){

   sspm = new sSPM();

   ttpm = new sTPM();

}

/**
 * construct a stripe object from a 2DM on the full lattice
 * @param ttpm input TPM
 */
Stripe::Stripe(const TPM &tpm){

   //first construct the full system 1DM
   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   sspm = new sSPM();
   sspm->transform(spm);

   ttpm = new sTPM();
   ttpm->transform(tpm);

}

/**
 * copy constructor 
 * @param mat_copy The stripe you want to be copied into the object you are constructing
 */
Stripe::Stripe(const Stripe &stripe_c){

   sspm = new sSPM(g1DM());
   ttpm = new sTPM(g2DM());

}

/**
 * Destructor
 */
Stripe::~Stripe(){

   delete sspm;
   delete ttpm;

}

/**
 * @return the 1DM
 */
const sSPM &Stripe::g1DM() const {

   return *sspm;

}

/**
 * @return the 2DM
 */
const sTPM &Stripe::g2DM() const {

   return *ttpm;

}

ostream &operator<<(ostream &output,const Stripe &stripe_p){
   
   output << std::endl;
   output << "1DM" << std::endl;
   output << std::endl;
   output << stripe_p.g1DM() << std::endl;
   output << std::endl;
   output << "2DM" << std::endl;
   output << std::endl;
   output << stripe_p.g2DM() << std::endl;

   return output;

}

/**
 * @return the number of susbsystem particles
 */
double Stripe::gN() const{

   double ward = 0.0;

   for(int i = 0;i < Tools::gL();++i)
      ward += (*sspm)[i];

   return 2.0 * ward;

}

/**
 * @return the number of susbsystem pairs
 */
double Stripe::gnpairs() const{

   return ttpm->trace();

}

/**
 * construct the 1D hubbard subsystem hamiltonian
 * @param U on-site repulsion
 */
void Stripe::hubbard1D(double U){

   sspm->hubbard1D_kin();
   ttpm->hubbard1D_rep(U);

}

/**
 * @return the inproduct of the 1DM and 2DM parts of two Stripe objects
 */
double Stripe::ddot(const Stripe &stripe_i) const {

   return sspm->ddot(stripe_i.g1DM()) + ttpm->ddot(stripe_i.g2DM());

}
