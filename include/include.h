//nog enkele definities:
#ifdef PQG

#define __G_CON

#endif

#ifdef PQGT1

#define __G_CON
#define __T1_CON

#endif

#include "lapack.h"
#include "Matrix.h"
#include "BlockMatrix.h"
#include "Vector.h"
#include "BlockVector.h"
#include "TPM.h"
#include "SPM.h"
#include "PHM.h"
#include "DPM.h"

#include "SUP.h"
#include "EIG.h"

#include "Hamiltonian.h"
