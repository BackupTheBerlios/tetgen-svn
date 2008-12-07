// g++ -o tsample -I../src -g tsample.cxx -L../src -ltet 
/*
  samples    Given a set of inequalities: Ax < b (or Ax <= b) on a domain,
             generate a set of addmissible points in file out.smesh.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tetgen.h"

#define REAL double

REAL ccwtest(REAL *pa, REAL *pb, REAL *pc, REAL *pd)
{
  REAL adx, bdx, cdx;
  REAL ady, bdy, cdy;
  REAL adz, bdz, cdz;

  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];
  adz = pa[2] - pd[2];
  bdz = pb[2] - pd[2];
  cdz = pc[2] - pd[2];

  return adx * (bdy * cdz - bdz * cdy)
       + bdx * (cdy * adz - cdz * ady)
       + cdx * (ady * bdz - adz * bdy);
}

int main(int argc, char *argv[])
{
  double xmin, ymin, zmin; /* The domain */
  double xmax, ymax, zmax; 
  double A[3], B[3], C[3], O[3]; /* Coordinates */
  double P[3], Q[3], R[3];
  double stepx, stepy;
  double s4, s5, s6, s7, s8, s9;
  int i, j;

  tetgenmesh::arraypool *ptarray; /* Addimissible points */
  double *coord;

  xmin = 0;
  ymin = -11;
  xmax = 13;
  ymax = 0;

  /* A, B, C are fixed. */
  A[0] = 4.5;
  A[1] = -6;
  A[2] = 0;

  B[0] = 8.5;
  B[1] = -7;
  B[2] = 0;

  C[0] = 7.5;
  C[1] = -3;
  C[2] = 0;

  O[0] = 7;
  O[1] = -7;
  O[2] = 5;

  /* Spacing between two sampled points */
  stepx = 0.2;
  stepy = 0.2;

  /* The array of addimisible points. Block size 1024.*/
  ptarray = new tetgenmesh::arraypool(sizeof(double) * 3, 10);

  R[2] = 0;

  /* Pick a point P. */
  P[0] = 1.4;
  P[1] = -6.2;
  P[2] = 0;

  /* Let Q be fixed for case (1). */
  Q[0] = 1.4;
  Q[1] = -8;
  Q[2] = 0;

  printf("  Sample R-points for P (%g, %g) and Q (%g, %g).\n", P[0], P[1],
         Q[0], Q[1]);

  i = 0;
  do {
    R[1] = ymin + ((double) i) * stepy;
    j = 0;
    do {
      R[0] = xmin + ((double) j) * stepx;
      s4 = ccwtest(P, Q, R, O);
      if (s4 < 0) {
        s5 = ccwtest(A, B, O, Q);
        if (s5 < 0) {
          s6 = ccwtest(A, B, O, R);
          if (s6 >= 0) {
            s7 = ccwtest(Q, R, O, A);
            if (s7 < 0) {
              s8 = ccwtest(Q, R, O, C);
              if (s8 >= 0) {
                s9 = ccwtest(C, A, O, R);
                if (s9 >= 0) {
                  ptarray->newindex((void **) &coord);
                  coord[0] = R[0];
                  coord[1] = R[1];
                  coord[2] = R[2];
                }
              }
            }
          }
        } else {
          printf("  !! Q (%g, %g) is not addmissible for case (1).\n", 
                 Q[0], Q[1]);
          exit(1);
        }
      }
      j++;
    } while (R[0] < xmax);
    i++;
  } while (R[1] < ymax);

  printf("  Found %ld addimissible R-points.\n", ptarray->objects);

  printf("  Writing points to samples.smesh\n");
  FILE * fout;
  fout = fopen("tsamples.smesh", "w");
  fprintf(fout, "%d  3  0  0\n", (int) (9 + ptarray->objects));

  /* Output the bounding box points. */
  fprintf(fout, "1  %g %g %g  # A\n", xmin, ymin, 0);
  fprintf(fout, "2  %g %g %g  # A\n", xmax, ymin, 0);
  fprintf(fout, "3  %g %g %g  # A\n", xmax, ymax, 0);
  fprintf(fout, "4  %g %g %g  # A\n", xmin, ymax, 0);

  /* Output A, B, C, P, Q */
  fprintf(fout, "5  %g %g %g  # A\n", A[0], A[1], A[2]);
  fprintf(fout, "6  %g %g %g  # B\n", B[0], B[1], B[2]);
  fprintf(fout, "7  %g %g %g  # C\n", C[0], C[1], C[2]);
  fprintf(fout, "8  %g %g %g  # P\n", P[0], P[1], P[2]);
  fprintf(fout, "9  %g %g %g  # Q\n", Q[0], Q[1], Q[2]);

  /* Output sample points */
  for (i = 0; i < ptarray->objects; i++) {
    coord = (double *) fastlookup(ptarray, i);
    fprintf(fout, "%d  %g %g %g.\n", i + 10, coord[0], coord[1], coord[2]);
  }
  
  fprintf(fout, "\n");
  fprintf(fout, "1  0\n");
  fprintf(fout, "3  5 6 7\n");
  fprintf(fout, "0\n");
  fprintf(fout, "0\n");

  fclose(fout);

  delete ptarray; 
}
