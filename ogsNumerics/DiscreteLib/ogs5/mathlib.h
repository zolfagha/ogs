/***************************************************************************
   ROCKFLOW - Modul: mathlib.h

   Aufgabe:
   ROCKFLOW-Schnittstelle fuer alle mathematischen Operationen, die nicht
   Standard (nicht in math.h enthalten) sind.
   mathlib.h benutzt externe Bibliotheken, die bei Bedarf ausgetauscht
   werden koennen, ohne die ROCKFLOW-Schnittstelle aendern zu muessen.

 **************************************************************************/

#ifndef mathlib_INC
#define mathlib_INC
/* Schutz gegen mehrfaches Einfuegen */

#include <cstddef>

#define noTESTMATH

/* Andere oeffentlich benutzte Module */
//#include "test.h"

/* Die Schnittstellen der Gleichungsloeser und der Speichertechnik werden
   sozusagen durchgeschleift: */
/* #include "matrix.h" */
/* Speichertechnik fuer Matrix des Gesamtgleichungssystems */
/* Iterative GLS-Loeser auf Basis der Speichertechnik aus 'matrix.h' (herkoemmliche Verfahren) */

// C++
//#include <vector>

/*##########################################################################
      Mathematische Funktionen
 ########################################################################*/


/*   MNulleVec             - Setze angegebenen Vektor = 0.0 */
extern void MNulleVec (double* vec, long g);
/*   MNulleMat             - Setze angegebene Matrix = 0.0 */
extern void MNulleMat (double* vec, long m, long n);

/*  Berechnet Gradient der 3D Testfunktionen */
extern double MXPGaussPkt(long grd, long pkt);
/* Punkte fuer die X Punkt Gauss-Integration */
extern double MXPGaussFkt(long grd, long pkt);

extern void  realCoordTriHQ(double* x, const double* XY, const double* u );

// Family of  element interpolation. WW
extern void ShapeFunctionLine(double* N1, const double* u);
extern void ShapeFunctionLineHQ(double* N1, const double* u);
extern void SamplePointTriHQ(const int nsample, double* SPoints);
extern void SamplePointTet5(const int nsample, double* SPoints);
extern void SamplePointTet15(const int nsample, double* SPoints);
extern void SamplePointPyramid5(const int nsample, double* SPoints);
extern void SamplePointPyramid8(const int i, double* SPoint);
extern void SamplePointPyramid13(const int nsample, double* SPoints);
extern void ShapeFunctionTri(double* N3, const double* u);
extern void ShapeFunctionTriHQ(double* N6, const double* u);
extern void ShapeFunctionQuad(double* N4, const double* u);
extern void ShapeFunctionQuadHQ8(double* N8, const double* u);
extern void ShapeFunctionQuadHQ(double* N9, const double* u);
extern void ShapeFunctionTet(double* Nt4, const double* u);
extern void ShapeFunctionTetHQ(double* N10, const double* u);
extern void ShapeFunctionHex(double* N8, const double* x);
extern void ShapeFunctionHexHQ(double* N9, const double* u);
extern void ShapeFunctionPri(double* N, const double* x);
extern void ShapeFunctionPriHQ(double* N, const double* u);
extern void ShapeFunctionPyra(double* N, const double* x);
extern void ShapeFunctionPyraHQ13(double* N, const double* u);
// Gradient of ...
extern void GradShapeFunctionLine(double* dN1, const double* u);
extern void GradShapeFunctionLineHQ(double* dN1, const double* u);
extern void GradShapeFunctionTri(double* dN3, const double* u);
extern void GradShapeFunctionTriHQ(double* dN3, const double* u);
extern void GradShapeFunctionQuad(double* dN4, const double* u);
extern void GradShapeFunctionQuadHQ(double* dN9, const double* u);
extern void GradShapeFunctionQuadHQ8(double* dN8, const double* u);
extern void GradShapeFunctionTet(double* dNt4, const double* u);
extern void GradShapeFunctionTetHQ(double* dN10, const double* u);
extern void GradShapeFunctionHex(double* N8, const double* x);
extern void GradShapeFunctionHexHQ(double* dN9, const double* u);
extern void GradShapeFunctionPri(double* dN, const double* x);
extern void GradShapeFunctionPriHQ(double* dN, const double* u);
extern void GradShapeFunctionPyra(double* dN, const double* x);
extern void GradShapeFunctionPyraHQ13(double* dN, const double* u);

extern double pai;

extern long binarySearch(long* arr, long target, long start, long end);
//WW Cubic spline
//WW
double ComputeDetTri(const double* x1, const double* x2,
                     const double* x3);
double ComputeDetTex(const double* x1, const double* x2,
                     const double* x3, const double* x4);
void CrossProduction(const double* x, const double* y, double* z);
double NormalizeVector(double* x, size_t n);

//extern double MVectorlength(double dx, double dy, double dz);
extern double PointProduction(double* x, double* y);
extern void VCopy(double* x, const double* y, const int n);

//NW
extern double MLangevin(double v);
extern double MinMod(double v1, double v2);
extern double SuperBee(double v1, double v2);
extern double GetFCTADiff(double v1, double v2);

#endif                                            /* gehoert zum Schutz gegen mehrfaches Einfuegen */

/*##########################################################################
    Ende von: ROCKFLOW - Modul mathlib.h
 ########################################################################*/
