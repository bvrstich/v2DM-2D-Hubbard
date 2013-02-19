#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <cstring>

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
   cout.precision(15);

   char *filename = 0;

   struct option long_options[] =
   {
      {"file",  required_argument, 0, 'f'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;
   while( (j = getopt_long (argc, argv, "hf:", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -f, --file=file.h5           file to read (in HDF5 format)\n"
               "    -h, --help                   Display this help\n"
               "\n"
               "By default, the program rebuilds the hamiltonian matrix unless -s is given. The file\n"
               "and setupfile can be the same.\n"
               "\n";
            return 0;
            break;
         case 'f':
            filename = optarg;
            if( strlen(optarg) <= 0 )
            {
               std::cerr << "Couldn't find a filename. Please specify the file to read." << endl;
               return -1;
            }
            break;
      }

   if(!filename)
   {
      cout << "Usage: " << argv[0] << " [OPTIONS]\n"
         "\n"
         "    -f, --file=file.h5           file to read (in HDF5 format)\n"
         "    -h, --help                   Display this help\n"
         "\n";
      return 0;
   }

   cout << "Reading: " << filename << endl;

   int L,N;
   double U;

   TPM::ReadInitfromFile(filename,L,N,U);

   U = 4;

   cout << "Found L=" << L << " N=" << N << " U=" << U << endl;
   cout << endl;

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

   TPM ham;
   ham.hubbard(U);

   TPM rdm;
   rdm.ReadfromFile(filename);

   cout << "Energy: " << ham.ddot(rdm) << endl;

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
