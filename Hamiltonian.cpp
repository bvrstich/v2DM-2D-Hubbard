#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int **Hamiltonian::xy_a;
int **Hamiltonian::a_xy;

int Hamiltonian::L;
int Hamiltonian::M;

/**
 * function that allocates and constructs the lists.
 * @param L_in the dimension of the square lattice.
 */
void Hamiltonian::init(int L_in){

   L = L_in;
   M = L*L;

   //allocate
   xy_a = new int * [L];

   for(int x = 0;x < L;++x)
      xy_a[x] = new int [L];

   a_xy = new int * [M];

   for(int i = 0;i < M;++i)
      a_xy[i] = new int [2];

   //construct:
   int a = 0;

   for(int x = 0;x < L;++x)
      for(int y = 0;y < L;++y){

         a_xy[a][0] = x;
         a_xy[a][1] = y;

         xy_a[x][y] = a;

         a++;

      }

}

/**
 * deallocates the lists.
 */
void Hamiltonian::clear(){

   //delete xy_a
   for(int x = 0;x < L;++x)
      delete [] xy_a[x];

   delete [] xy_a;

   //delete a_xy
   for(int i = 0;i < M;++i)
      delete [] a_xy[i];

   delete [] a_xy;

}

/**
 * print the list
 */
void Hamiltonian::print(){

   for(int i = 0;i < M;++i)
      std::cout << i << "\t" << a_xy[i][0]<< "\t" << a_xy[i][1] << std::endl;

}

/**
 * access to the a_xy list from outside this class
 * @param a the sp-index
 * @param option can be 0 or 1
 * @return k_x if option == 0 , k_y if option == 1
 */
int Hamiltonian::ga_xy(int a,int option){

   return a_xy[a][option];

}

/**
 * access to the list xy_a from outside this class
 * @param k_x the x momentum
 * @param k_y the y momentum
 * @return the sp index corresponding to k_x and k_y
 */
int Hamiltonian::gxy_a(int k_x,int k_y){
   
   return xy_a[k_x][k_y];

}
