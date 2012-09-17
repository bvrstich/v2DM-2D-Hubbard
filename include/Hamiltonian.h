#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <cstdlib>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 28-02-2011\n\n
 * This is a class written for the 2 dimension hubbard model, it translates the physical degrees of freedom of the
 * 2D hubbard model: (x,y,sigma) to the one particle basis of the regular program.
 */

class Hamiltonian{

   public:

      //initializes the lists.
      static void init();

      //clears the lists;
      static void clear();

      //print the list
      static void print();

      //access the lists from outside the class
      static int gxy_a(int,int);

      //access the lists from outside the class
      static int ga_xy(int,int);

      static int bar(int);

   private:
      
      //!static list that translates the two indices of the 2D hubbard model to one sp index.
      static int **xy_a;

      //!static list that translates the sp index alpha to the two physical indices of the 2D hubbard model.
      static int **a_xy;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
