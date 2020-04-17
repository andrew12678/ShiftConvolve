#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void computeMGF(void *, void *, void *, void *, void *, void *);
extern void fftconvPairs(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fullconvolvePaired(void *, void *, void *);
extern void fullconvolvePairedLog(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"computeMGF",            (DL_FUNC) &computeMGF,            6},
    {"fftconvPairs",          (DL_FUNC) &fftconvPairs,          9},
    {"fullconvolvePaired",    (DL_FUNC) &fullconvolvePaired,    3},
    {"fullconvolvePairedLog", (DL_FUNC) &fullconvolvePairedLog, 3},
    {NULL, NULL, 0}
};

void R_init_ShiftConvolvePoibin(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
