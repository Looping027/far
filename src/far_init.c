#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void CVkerfon(void *, void *, void *, void *, void *, void *, void *);
extern void prevkerfon(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"CVkerfon",   (DL_FUNC) &CVkerfon,   7},
  {"prevkerfon", (DL_FUNC) &prevkerfon, 7},
  {NULL, NULL, 0}
};

void R_init_far(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
