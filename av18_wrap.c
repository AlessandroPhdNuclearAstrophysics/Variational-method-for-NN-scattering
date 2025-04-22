#include "/usr/local/Wolfram/Mathematica/13.1/SystemFiles/IncludeFiles/C/WolframLibrary.h"
#include <stdio.h>

/* required by every LibraryLink library */
DLLEXPORT mint WolframLibrary_getVersion() {
  return WolframLibraryVersion;
}
DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) {
  return 0;
}

/* the new wrapper */
DLLEXPORT int av18pw90_LL(
    WolframLibraryData libData,
    mint              Argc,
    MArgument         *Args,
    MArgument         Res
) {
  /* pull out the 7 ints and the real scalar */
  mint lpot  = MArgument_getInteger(Args[0]);
  mint l     = MArgument_getInteger(Args[1]);
  mint s     = MArgument_getInteger(Args[2]);
  mint j     = MArgument_getInteger(Args[3]);
  mint t     = MArgument_getInteger(Args[4]);
  mint t1z   = MArgument_getInteger(Args[5]);
  mint t2z   = MArgument_getInteger(Args[6]);
  mreal r    = MArgument_getReal(Args[7]);

  /* grab the MTensor for the 2×2 output */
  MTensor vpwTensor = MArgument_getMTensor(Args[8]);
  mreal *vpwData    = libData->MTensor_getRealData(vpwTensor);

  /* call the Fortran routine (pass pointers!) */
  av18pw90_(&lpot, &l, &s, &j, &t, &t1z, &t2z, &r, vpwData);

  /* since you passed the MTensor in with "Shared", it’s already visible
     on the Mathematica side – no need to MArgument_setMTensor(Res,…) for Void return. */
  return LIBRARY_NO_ERROR;
}
