#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

using std::ostream;
using std::cout;
using std::endl;
using std::complex;

#include "include.h"

/**
 * constructor, makes empty Matrix of dimension Tools::gL()*Tools::gL()
 */
SPM_ft::SPM_ft() : Matrix(Tools::gL()*Tools::gL()) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM_ft::SPM_ft(const SPM_ft &spm_copy) : Matrix(spm_copy) { }

/**
 * construct a SPM_ft by Fourier Transforming a regular spm object
 * @param SPM object in momentum space
 */
SPM_ft::SPM_ft(const SPM &spm_i) : Matrix(Tools::gL()*Tools::gL()){

   complex<double> I(0.0,1.0);

   double ba_x;
   double ba_y;

   for(int a = 0;a < gn();++a)
      for(int b = 0;b < gn();++b){

         ba_x = Hamiltonian::ga_xy(b,0) - Hamiltonian::ga_xy(a,0);
         ba_y = Hamiltonian::ga_xy(b,1) - Hamiltonian::ga_xy(a,1);

         (*this)(a,b) = 0.0;

         complex<double> ward(0.0,0.0);

         for(int kx = 0;kx < Tools::gL();++kx){

            double kx_ = kx/(double)Tools::gL();

            for(int ky = 0;ky < Tools::gL();++ky){

               double ky_ = ky/(double)Tools::gL();

               ward += exp( 2.0 * M_PI * I *  kx_ * ba_x) * exp( 2.0 * M_PI * I * ky_ * ba_y) * spm_i[Hamiltonian::gxy_a(kx,ky)];

            }

         }

         cout << a << "\t" << b << "\t" << ward << endl;

      }

}

/**
 * destructor
 */
SPM_ft::~SPM_ft(){ }

ostream &operator<<(ostream &output,const SPM_ft &spm_p){

   for(int a = 0;a < spm_p.gn();++a)
      for(int b = 0;b < spm_p.gn();++b){

         output << a << "\t" << b << "\t|\t(" << Hamiltonian::ga_xy(a,0) << "," << Hamiltonian::ga_xy(a,1)
         
            << ")\t(" << Hamiltonian::ga_xy(b,0) << "," << Hamiltonian::ga_xy(b,1) << ")\t" << spm_p(a,b) << endl;

      }

   return output;

}
