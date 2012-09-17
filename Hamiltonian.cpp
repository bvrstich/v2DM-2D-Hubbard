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

/**
 * function that allocates and constructs the lists.
 */
void Hamiltonian::init(){

   //allocate
   xy_a = new int * [Tools::gL()];

   for(int x = 0;x < Tools::gL();++x)
      xy_a[x] = new int [Tools::gL()];

   a_xy = new int * [Tools::gM()];

   for(int i = 0;i < Tools::gM();++i)
      a_xy[i] = new int [2];

   //construct:
   int a = 0;

   for(int x = 0;x < Tools::gL();++x)
      for(int y = 0;y < Tools::gL();++y){

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
   for(int x = 0;x < Tools::gL();++x)
      delete [] xy_a[x];

   delete [] xy_a;

   //delete a_xy
   for(int i = 0;i < Tools::gM();++i)
      delete [] a_xy[i];

   delete [] a_xy;

}

/**
 * print the list
 */
void Hamiltonian::print(){

   for(int i = 0;i < Tools::gM();++i)
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

/**
 * transform the "particle momentum" sp-index to the "hole momentum" sp-index
 * @param a the input sp-index
 */
int Hamiltonian::bar(int a){

   return xy_a[(-a_xy[a][0] + Tools::gL())%Tools::gL()][(-a_xy[a][1] + Tools::gL())%Tools::gL()];

}
