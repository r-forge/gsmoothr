
#include "R_ext/Rdynload.h"
#include "FIRMAGene.h"

static const R_CMethodDef cMethods[]  = {
  {"muf",(DL_FUNC)&muf,3},
  {NULL, NULL, 0}
};

void R_init_flagme(DllInfo *info){
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
