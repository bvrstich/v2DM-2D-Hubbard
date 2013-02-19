#ifndef STRIPE_H
#define STRIPE_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class Vector;

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * This is a class written to contain the 1 and 2DM of a stripe of a 2D lattice. These two objects
 * suffice to describe the fractional N system completely.
 */

class Stripe{

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param stripe_p de Stripe you want to print
    */
   friend ostream &operator<<(ostream &output,const Stripe &stripe_p);

   public:

      //constructor
      Stripe();

      Stripe(const TPM &);

      //copy constructor
      Stripe(const Stripe &);

      //destructor
      virtual ~Stripe();

      void hubbard1D(double);

      double gN() const;
      
      double gnpairs() const;

      const sSPM &g1DM() const;

      const sTPM &g2DM() const;

      double ddot(const Stripe &) const;

   private:
      
      //!striped SPM object
      sSPM *sspm;

      //!striped TPM object
      sTPM *ttpm;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
