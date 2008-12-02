#ifndef geomCXX
#define geomCXX

#include "tetgen.h"

REAL tetgenmesh::PI = 3.14159265358979323846264338327950288419716939937510582;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_edge_2d()    Triangle-edge coplanar intersection test.                //
//                                                                           //
// This routine takes a triangle T (with vertices A, B, C) and an edge E (P, //
// Q) in a plane in 3D, and tests if they intersect each other.  Return 1 if //
// they are intersected, i.e., T \cap E is not empty, otherwise, return 0.   //
//                                                                           //
// If the point 'R' is not NULL, it lies strictly above T [A, B, C].         //
//                                                                           //
// If T1 and T2 intersect each other (return 1), they may intersect in diff- //
// erent ways. If 'level' > 0, their intersection type will be reported in   //
// combinations of 'types' and 'pos'.                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::tri_edge_2d(point A, point B, point C, point P, point Q, 
  point R, int level, int *types, int *pos)
{
  point U[3], V[3];  // The permuted vectors of points.
  int pu[3], pv[3];  // The original positions of points.
  REAL sA, sB, sC;
  REAL s1, s2, s3, s4;
  int z1;

  if (R == NULL) {
    REAL n[3], len;
    // Calculate a lift point, saved in dummypoint.
    facenormal(A, B, C, n, 1);
    len = sqrt(DOT(n, n));
    n[0] /= len;
    n[1] /= len;
    n[2] /= len;
    len = DIST(A, B);
    len += DIST(B, C);
    len += DIST(C, A);
    len /= 3.0;
    R = dummypoint;
    R[0] = A[0] + len * n[0];
    R[1] = A[1] + len * n[1];
    R[2] = A[2] + len * n[2];
  }

  // Test A's, B's, and C's orientations wrt plane PQR. 
  sA = orient3d(P, Q, R, A);
  sB = orient3d(P, Q, R, B);
  sC = orient3d(P, Q, R, C);
  orient3dcount+=3;

  if (b->epsilon) {
    // Re-evaluate the sign with respect to the tolerance.
    if ((sA != 0) && iscoplanar(P, Q, R, A, sA)) sA = 0;
    if ((sB != 0) && iscoplanar(P, Q, R, B, sB)) sB = 0;
    if ((sC != 0) && iscoplanar(P, Q, R, C, sC)) sC = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-edge-2d (%d %d %d)-(%d %d)-(%d) (%c%c%c)", pointmark(A),
      pointmark(B), pointmark(C), pointmark(P), pointmark(Q), pointmark(R),
      sA > 0 ? '+' : (sA < 0 ? '-' : '0'), sB>0 ? '+' : (sB<0 ? '-' : '0'),
      sC>0 ? '+' : (sC<0 ? '-' : '0'));
  }
  triedgcopcount++;

  if (sA < 0) {
    if (sB < 0) {
      if (sC < 0) { // (---).
        return 0; 
      } else {
        if (sC > 0) { // (--+).
          // All points are in the right positions.
          SETVECTOR3(U, A, B, C);  // I3
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 0, 1, 2);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else { // (--0).
          SETVECTOR3(U, A, B, C);  // I3
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 0, 1, 2);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        }
      }
    } else { 
      if (sB > 0) {
        if (sC < 0) { // (-+-).
          SETVECTOR3(U, C, A, B);  // PT = ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 2, 0, 1);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else {
          if (sC > 0) { // (-++).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 0;
          } else { // (-+0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          }
        }
      } else {
        if (sC < 0) { // (-0-).
          SETVECTOR3(U, C, A, B);  // PT = ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 2, 0, 1);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        } else {
          if (sC > 0) { // (-0+).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 2;
          } else { // (-00).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 3; 
          }
        }
      }
    }
  } else {
    if (sA > 0) {
      if (sB < 0) {
        if (sC < 0) { // (+--).
          SETVECTOR3(U, B, C, A);  // PT = ST x ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 1, 2, 0);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else {
          if (sC > 0) { // (+-+).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 0;
          } else { // (+-0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 2;
          }
        }
      } else { 
        if (sB > 0) {
          if (sC < 0) { // (++-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 0;
          } else {
            if (sC > 0) { // (+++).
              return 0; 
            } else { // (++0).
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 0, 1, 2);
              SETVECTOR3(pv, 1, 0, 2);
              z1 = 1; 
            }
          }
        } else { // (+0#)
          if (sC < 0) { // (+0-).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          } else {
            if (sC > 0) { // (+0+).
              SETVECTOR3(U, C, A, B);  // PT = ST
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 2, 0, 1);
              SETVECTOR3(pv, 1, 0, 2);
              z1 = 1;
            } else { // (+00).
              SETVECTOR3(U, B, C, A);  // PT = ST x ST
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 1, 2, 0);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            }
          }
        }
      }
    } else { 
      if (sB < 0) {
        if (sC < 0) { // (0--).
          SETVECTOR3(U, B, C, A);  // PT = ST x ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 1, 2, 0);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        } else {
          if (sC > 0) { // (0-+).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          } else { // (0-0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 3; 
          }
        }
      } else { 
        if (sB > 0) {
          if (sC < 0) { // (0+-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 2;
          } else {
            if (sC > 0) { // (0++).
              SETVECTOR3(U, B, C, A);  // PT = ST x ST
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 1, 2, 0);
              SETVECTOR3(pv, 1, 0, 2);
              z1 = 1;
            } else { // (0+0).
              SETVECTOR3(U, C, A, B);  // PT = ST
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 2, 0, 1);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            }
          }
        } else { // (00#)
          if (sC < 0) { // (00-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2);
            z1 = 3; 
          } else {
            if (sC > 0) { // (00+).
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 0, 1, 2);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            } else { // (000)
              // Not possible unless ABC is degenerate.
              z1 = 4;
            }
          }
        }
      }
    }
  }

  s1 = orient3d(U[0], U[2], R, V[1]);  // A, C, R, Q
  s2 = orient3d(U[1], U[2], R, V[0]);  // B, C, R, P
  orient3dcount+=2;

  if (b->epsilon) {
    if ((s1 != 0) && iscoplanar(U[0], U[2], R, V[1], s1)) s1 = 0;
    if ((s2 != 0) && iscoplanar(U[1], U[2], R, V[0], s2)) s2 = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-edge-2d (%d %d %d)-(%d %d %d) (%d) (%c%c)\n",
      pointmark(U[0]), pointmark(U[1]), pointmark(U[2]), pointmark(V[0]),
      pointmark(V[1]), pointmark(V[2]), z1, s1>0 ? '+' : (s1<0 ? '-' : '0'),
      s2>0 ? '+' : (s2<0 ? '-' : '0'));
  }
  assert(z1 != 4); // SELF_CHECK

  if (s1 > 0) {
    return 0;
  }
  if (s2 < 0) {
    return 0;
  }

  if (level == 0) {
    return 1;  // They are intersected.
  }

  if (z1 == 1) {
    if (s1 == 0) {  // (0###)
      // C = Q.
      types[0] = (int) SHAREVERT;
      pos[0] = pu[2]; // C
      pos[1] = pv[1]; // Q
      types[1] = (int) DISJOINT; // No 2nd intersection
    } else {
      if (s2 == 0) { // (#0##)
        // C = P.
        types[0] = (int) SHAREVERT;
        pos[0] = pu[2]; // C
        pos[1] = pv[0]; // P
        types[1] = (int) DISJOINT; // No 2nd intersection
      } else { // (-+##)
        // C in [P, Q].
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[2]; // C
        pos[1] = pv[0]; // [P, Q]
        types[1] = (int) DISJOINT; // No 2nd intersection
      }
    }
    return 1;
  }

  s3 = orient3d(U[0], U[2], R, V[0]);  // A, C, R, P
  s4 = orient3d(U[1], U[2], R, V[1]);  // B, C, R, Q
  orient3dcount+=2;
      
  if (b->epsilon) {
    if ((s3 != 0) && iscoplanar(U[0], U[2], R, V[0], s3)) s3 = 0;
    if ((s4 != 0) && iscoplanar(U[1], U[2], R, V[1], s4)) s4 = 0;
  }

  if (z1 == 0) {  // (tritri-03)
    if (s1 < 0) {
      if (s3 > 0) {
        assert(s2 > 0); // SELF_CHECK
        if (s4 > 0) {
          // [P, Q] overlaps [k, l] (-+++).
          types[0] = (int) ACROSSEDGE;
          pos[0] = pu[2]; // [C, A]
          pos[1] = pv[0]; // [P, Q]
          types[1] = (int) TOUCHFACE;
          pos[2] = 3;     // [A, B, C]
          pos[3] = pv[1]; // Q
        } else {
          if (s4 == 0) {
            // Q = l, [P, Q] contains [k, l] (-++0).
            types[0] = (int) ACROSSEDGE;
            pos[0] = pu[2]; // [C, A]
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) TOUCHEDGE;
            pos[2] = pu[1]; // [B, C]
            pos[3] = pv[1]; // Q
          } else { // s4 < 0
            // [P, Q] contains [k, l] (-++-).
            types[0] = (int) ACROSSEDGE;
            pos[0] = pu[2]; // [C, A]
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) ACROSSEDGE;
            pos[2] = pu[1]; // [B, C]
            pos[3] = pv[0]; // [P, Q]
          }
        }
      } else {
        if (s3 == 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // P = k, [P, Q] in [k, l] (-+0+).
            types[0] = (int) TOUCHEDGE;
            pos[0] = pu[2]; // [C, A]
            pos[1] = pv[0]; // P
            types[1] = (int) TOUCHFACE;
            pos[2] = 3;     // [A, B, C]
            pos[3] = pv[1]; // Q
          } else {
            if (s4 == 0) {
              // [P, Q] = [k, l] (-+00).
              types[0] = (int) TOUCHEDGE;
              pos[0] = pu[2]; // [C, A]
              pos[1] = pv[0]; // P
              types[1] = (int) TOUCHEDGE;
              pos[2] = pu[1]; // [B, C]
              pos[3] = pv[1]; // Q
            } else {
              // P = k, [P, Q] contains [k, l] (-+0-).
              types[0] = (int) TOUCHEDGE;
              pos[0] = pu[2]; // [C, A]
              pos[1] = pv[0]; // P
              types[1] = (int) ACROSSEDGE;
              pos[2] = pu[1]; // [B, C]
              pos[3] = pv[0]; // [P, Q]
            }
          }
        } else { // s3 < 0
          if (s2 > 0) {
            if (s4 > 0) {
              // [P, Q] in [k, l] (-+-+).
              types[0] = (int) TOUCHFACE;
              pos[0] = 3;     // [A, B, C]
              pos[1] = pv[0]; // P
              types[1] = (int) TOUCHFACE;
              pos[2] = 3;     // [A, B, C]
              pos[3] = pv[1]; // Q
            } else {
              if (s4 == 0) {
                // Q = l, [P, Q] in [k, l] (-+-0).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // [A, B, C]
                pos[1] = pv[0]; // P
                types[1] = (int) TOUCHEDGE;
                pos[2] = pu[1]; // [B, C]
                pos[3] = pv[1]; // Q
              } else { // s4 < 0
                // [P, Q] overlaps [k, l] (-+--).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // [A, B, C]
                pos[1] = pv[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = pu[1]; // [B, C]
                pos[3] = pv[0]; // [P, Q]
              }
            }
          } else { // s2 == 0
            // P = l (#0##).
            types[0] = (int) TOUCHEDGE;
            pos[0] = pu[1]; // [B, C]
            pos[1] = pv[0]; // P
            types[1] = (int) DISJOINT; // No 2nd intersection
          }
        }
      }
    } else { // s1 == 0
      // Q = k (0####)
      types[0] = (int) TOUCHEDGE;
      pos[0] = pu[2]; // [C, A]
      pos[1] = pv[1]; // Q
      types[1] = (int) DISJOINT; // No 2nd intersection
    }
  } else if (z1 == 2) {  // (tritri-23)
    if (s1 < 0) {
      if (s3 > 0) {
        assert(s2 > 0); // SELF_CHECK
        if (s4 > 0) {
          // [P, Q] overlaps [A, l] (-+++).
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[0]; // A
          pos[1] = pv[0]; // [P, Q]
          types[1] = (int) TOUCHFACE;
          pos[2] = 3;     // [A, B, C]
          pos[3] = pv[1]; // Q
        } else {
          if (s4 == 0) {
            // Q = l, [P, Q] contains [A, l] (-++0).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) TOUCHEDGE;
            pos[2] = pu[1]; // [B, C]
            pos[3] = pv[1]; // Q
          } else { // s4 < 0
            // [P, Q] contains [A, l] (-++-).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) ACROSSEDGE;
            pos[2] = pu[1]; // [B, C]
            pos[3] = pv[0]; // [P, Q]
          }
        }
      } else {
        if (s3 == 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // P = A, [P, Q] in [A, l] (-+0+).
            types[0] = (int) SHAREVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // P
            types[1] = (int) TOUCHFACE;
            pos[2] = 3;     // [A, B, C]
            pos[3] = pv[1]; // Q
          } else {
            if (s4 == 0) {
              // [P, Q] = [A, l] (-+00).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pv[0]; // P
              types[1] = (int) TOUCHEDGE;
              pos[2] = pu[1]; // [B, C]
              pos[3] = pv[1]; // Q
            } else { // s4 < 0
              // Q = l, [P, Q] in [A, l] (-+0-).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pv[0]; // P
              types[1] = (int) ACROSSEDGE;
              pos[2] = pu[1]; // [B, C]
              pos[3] = pv[0]; // [P, Q]
            }
          }
        } else { // s3 < 0
          if (s2 > 0) {
            if (s4 > 0) {
              // [P, Q] in [A, l] (-+-+).
              types[0] = (int) TOUCHFACE;
              pos[0] = 3;     // [A, B, C]
              pos[1] = pv[0]; // P
              types[0] = (int) TOUCHFACE;
              pos[0] = 3;     // [A, B, C]
              pos[1] = pv[1]; // Q
            } else {
              if (s4 == 0) {
                // Q = l, [P, Q] in [A, l] (-+-0).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // [A, B, C]
                pos[1] = pv[0]; // P
                types[0] = (int) TOUCHEDGE;
                pos[0] = pu[1]; // [B, C]
                pos[1] = pv[1]; // Q
              } else { // s4 < 0
                // [P, Q] overlaps [A, l] (-+--).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // [A, B, C]
                pos[1] = pv[0]; // P
                types[0] = (int) ACROSSEDGE;
                pos[0] = pu[1]; // [B, C]
                pos[1] = pv[0]; // [P, Q]
              }
            }
          } else { // s2 == 0
            // P = l (#0##).
            types[0] = (int) TOUCHEDGE;
            pos[0] = pu[1]; // [B, C]
            pos[1] = pv[0]; // P
            types[1] = (int) DISJOINT; // No 2nd intersection
          }
        }
      }
    } else { // s1 == 0
      // Q = A (0###).
      types[0] = (int) SHAREVERT;
      pos[0] = pu[0]; // A
      pos[1] = pv[1]; // Q
      types[1] = (int) DISJOINT; // No 2nd intersection
    }
  } else if (z1 == 3) {  // (tritri-33)
    if (s1 < 0) {
      if (s3 > 0) {
        assert(s2 > 0); // SELF_CHECK
        if (s4 > 0) {
          // [P, Q] overlaps [A, B] (-+++).
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[0]; // A
          pos[1] = pv[0]; // [P, Q]
          types[1] = (int) TOUCHEDGE;
          pos[2] = pu[0]; // [A, B]
          pos[3] = pv[1]; // Q
        } else {
          if (s4 == 0) {
            // Q = B, [P, Q] contains [A, B] (-++0).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) SHAREVERT;
            pos[2] = pu[1]; // B
            pos[3] = pv[1]; // Q
          } else { // s4 < 0
            // [P, Q] contains [A, B] (-++-).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // [P, Q]
            types[1] = (int) ACROSSVERT;
            pos[2] = pu[1]; // B
            pos[3] = pv[0]; // [P, Q]
          }
        }
      } else {
        if (s3 == 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // P = A, [P, Q] in [A, B] (-+0+).
            types[0] = (int) SHAREVERT;
            pos[0] = pu[0]; // A
            pos[1] = pv[0]; // P
            types[1] = (int) TOUCHEDGE;
            pos[2] = pu[0]; // [A, B]
            pos[3] = pv[1]; // Q
          } else {
            if (s4 == 0) {
              // [P, Q] = [A, B] (-+00).
              types[0] = (int) SHAREEDGE;
              pos[0] = pu[0]; // [A, B]
              pos[1] = pv[0]; // [P, Q]
              types[1] = (int) SHAREEDGE;
              pos[2] = pu[0]; // [A, B]
              pos[3] = pv[0]; // [P, Q]
            } else { // s4 < 0
              // P= A, [P, Q] in [A, B] (-+0-).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pv[0]; // P
              types[1] = (int) ACROSSVERT;
              pos[2] = pu[1]; // B
              pos[3] = pv[0]; // [P, Q]
            }
          }
        } else { // s3 < 0
          if (s2 > 0) {
            if (s4 > 0) {
              // [P, Q] in [A, B] (-+-+).
              types[0] = (int) TOUCHEDGE;
              pos[0] = pu[0]; // [A, B]
              pos[1] = pv[0]; // P
              types[1] = (int) TOUCHEDGE;
              pos[2] = pu[0]; // [A, B]
              pos[3] = pv[1]; // Q
            } else {
              if (s4 == 0) {
                // Q = B, [P, Q] in [A, B] (-+-0).
                types[0] = (int) TOUCHEDGE;
                pos[0] = pu[0]; // [A, B]
                pos[1] = pv[0]; // P
                types[1] = (int) SHAREVERT;
                pos[2] = pu[1]; // B
                pos[3] = pv[1]; // Q
              } else { // s4 < 0
                // [P, Q] overlaps [A, B] (-+--).
                types[0] = (int) TOUCHEDGE;
                pos[0] = pu[0]; // [A, B]
                pos[1] = pv[0]; // P
                types[1] = (int) ACROSSVERT;
                pos[2] = pu[1]; // B
                pos[3] = pv[0]; // [P, Q]
              }
            }
          } else { // s2 == 0
            // P = B (#0##).
            types[0] = (int) SHAREVERT;
            pos[0] = pu[1]; // B
            pos[1] = pv[0]; // P
            types[1] = (int) DISJOINT; // No 2nd intersection
          }
        }
      }
    } else { // s1 == 0
      // Q = A (0###).
      types[0] = (int) SHAREVERT;
      pos[0] = pu[0]; // A
      pos[1] = pv[1]; // Q
      types[1] = (int) DISJOINT; // No 2nd intersection
    }
  }

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_edge_test()    Triangle-edge intersection test.                       //
//                                                                           //
// This routine takes a triangle T (with vertices A, B, C) and an edge E (P, //
// Q) in 3D, and tests if they intersect each other.  Return 1 if they are   //
// intersected, i.e., T \cap E is not empty, otherwise, return 0.            //
//                                                                           //
// If the point 'R' is not NULL, it lies strictly above the plane defined by //
// A, B, C. It is used in test when T and E are coplanar.                    //
//                                                                           //
// If T1 and T2 intersect each other (return 1), they may intersect in diff- //
// erent ways. If 'level' > 0, their intersection type will be reported in   //
// combinations of 'types' and 'pos'.                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::tri_edge_test(point A, point B, point C, point P, point Q, 
  point R, int level, int *types, int *pos)
{
  point U[3], V[3], Ptmp;
  int pu[3], pv[3], itmp;
  REAL sP, sQ, s1, s2, s3;
  int z1;

  // Test the locations of P and Q with respect to ABC.
  sP = orient3d(A, B, C, P);
  sQ = orient3d(A, B, C, Q);
  orient3dcount+=2;

  if (b->epsilon > 0) {
    if ((sP != 0) && iscoplanar(A, B, C, P, sP)) sP = 0;
    if ((sQ != 0) && iscoplanar(A, B, C, Q, sQ)) sQ = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-edge (%d %d %d)-(%d %d) (%c%c).\n", pointmark(A),
      pointmark(B), pointmark(C), pointmark(P), pointmark(Q),
      sP>0 ? '+' : (sP<0 ? '-' : '0'), sQ>0 ? '+' : (sQ<0 ? '-' : '0'));
  }
  triedgcount++;

  if (sP < 0) {
    if (sQ < 0) { // (--) disjoint
      return 0;
    } else {
      if (sQ > 0) { // (-+)
        SETVECTOR3(U, A, B, C);
        SETVECTOR3(V, P, Q, R);
        SETVECTOR3(pu, 0, 1, 2);
        SETVECTOR3(pv, 0, 1, 2);
        z1 = 0;
      } else { // (-0)
        SETVECTOR3(U, A, B, C);
        SETVECTOR3(V, P, Q, R);
        SETVECTOR3(pu, 0, 1, 2);
        SETVECTOR3(pv, 0, 1, 2);
        z1 = 1;
      }
    }
  } else {
    if (sP > 0) { // (+-)
      if (sQ < 0) {
        SETVECTOR3(U, A, B, C);
        SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
        SETVECTOR3(pu, 0, 1, 2);
        SETVECTOR3(pv, 1, 0, 2);
        z1 = 0;
      } else {
        if (sQ > 0) { // (++) disjoint
          return 0;
        } else { // (+0)
          SETVECTOR3(U, B, A, C); // A and B are flipped.
          SETVECTOR3(V, P, Q, R);
          SETVECTOR3(pu, 1, 0, 2);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        }
      }
    } else { // sP == 0
      if (sQ < 0) { // (0-)
        SETVECTOR3(U, A, B, C);
        SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
        SETVECTOR3(pu, 0, 1, 2);
        SETVECTOR3(pv, 1, 0, 2);
        z1 = 1;
      } else {
        if (sQ > 0) { // (0+)
          SETVECTOR3(U, B, A, C);  // A and B are flipped.
          SETVECTOR3(V, Q, P, R);  // P and Q are flipped.
          SETVECTOR3(pu, 1, 0, 2);
          SETVECTOR3(pv, 1, 0, 2);
          z1 = 1;
        } else { // (00)
          // A, B, C, P, and Q are coplanar.
          z1 = 2;
        }
      }
    }
  }

  if (z1 == 2) {
    // The triangle and the edge are coplanar.
    return tri_edge_2d(A, B, C, P, Q, R, level, types, pos);
  }

  s1 = orient3d(U[0], U[1], V[0], V[1]); orient3dcount++;
  if (b->epsilon) {
    if ((s1 != 0) && iscoplanar(U[0], U[1], V[0], V[1], s1)) s1 = 0;
  }
  if (s1 < 0) {
    return 0;
  }

  s2 = orient3d(U[1], U[2], V[0], V[1]); orient3dcount++;
  if (b->epsilon) {
    if ((s2 != 0) && iscoplanar(U[1], U[2], V[0], V[1], s2)) s2 = 0;
  }
  if (s2 < 0) {
    return 0;
  }

  s3 = orient3d(U[2], U[0], V[0], V[1]); orient3dcount++;
  if (b->epsilon) {
    if ((s3 != 0) && iscoplanar(U[2], U[0], V[0], V[1], s3)) s3 = 0;
  }
  if (s3 < 0) {
    return 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-edge (%d %d %d)-(%d %d) (%c%c%c).\n", pointmark(U[0]),
      pointmark(U[1]), pointmark(U[2]), pointmark(V[0]), pointmark(V[1]),
      s1>0 ? '+' : (s1<0 ? '-' : '0'), s2>0 ? '+' : (s2<0 ? '-' : '0'),
      s3>0 ? '+' : (s3<0 ? '-' : '0'));
  }

  if (level == 0) {
    return 1;  // The are intersected.
  }

  types[1] = (int) DISJOINT; // No second intersection point.

  if (z1 == 0) {
    if (s1 > 0) {
      if (s2 > 0) {
        if (s3 > 0) { // (+++)
          // [P, Q] passes interior of [A, B, C].
          types[0] = (int) ACROSSFACE;
          pos[0] = 3;  // interior of [A, B, C]
          pos[1] = 0;  // [P, Q]
        } else { // s3 == 0 (++0)
          // [P, Q] intersects [C, A].
          types[0] = (int) ACROSSEDGE;
          pos[0] = pu[2];  // [C, A]
          pos[1] = 0;  // [P, Q]
        }
      } else { // s2 == 0
        if (s3 > 0) { // (+0+)
          // [P, Q] intersects [B, C].
          types[0] = (int) ACROSSEDGE;
          pos[0] = pu[1];  // [B, C]
          pos[1] = 0;  // [P, Q]
        } else { // s3 == 0 (+00)
          // [P, Q] passes C.
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[2];  // C
          pos[1] = 0;  // [P, Q]
        }
      }
    } else { // s1 == 0
      if (s2 > 0) {
        if (s3 > 0) { // (0++)
          // [P, Q] intersects [A, B].
          types[0] = (int) ACROSSEDGE;
          pos[0] = pu[0];  // [A, B]
          pos[1] = 0;  // [P, Q]
        } else { // s3 == 0 (0+0)
          // [P, Q] passes A.
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[0];  // A
          pos[1] = 0;  // [P, Q]
        }
      } else { // s2 == 0
        if (s3 > 0) { // (00+)
          // [P, Q] passes B.
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[1];  // B
          pos[1] = 0;  // [P, Q]
        } else { // s3 == 0 (000)
          // Impossible.
          assert(0);
        }
      }
    }
  } else { // z1 == 1
    if (s1 > 0) {
      if (s2 > 0) {
        if (s3 > 0) { // (+++)
          // Q lies in [A, B, C].
          types[0] = (int) TOUCHFACE;
          pos[0] = 0; // [A, B, C]
          pos[1] = pv[1]; // Q
        } else { // s3 == 0 (++0)
          // Q lies on [C, A].
          types[0] = (int) TOUCHEDGE;
          pos[0] = pu[2]; // [C, A]
          pos[1] = pv[1]; // Q
        }
      } else { // s2 == 0
        if (s3 > 0) { // (+0+)
          // Q lies on [B, C].
          types[0] = (int) TOUCHEDGE;
          pos[0] = pu[1]; // [B, C]
          pos[1] = pv[1]; // Q
        } else { // s3 == 0 (+00)
          // Q = C.
          types[0] = (int) SHAREVERT;
          pos[0] = pu[2]; // C
          pos[1] = pv[1]; // Q
        }
      }
    } else { // s1 == 0
      if (s2 > 0) {
        if (s3 > 0) { // (0++)
          // Q lies on [A, B].
          types[0] = (int) TOUCHEDGE;
          pos[0] = pu[0]; // [A, B]
          pos[1] = pv[1]; // Q
        } else { // s3 == 0 (0+0)
          // Q = A.
          types[0] = (int) SHAREVERT;
          pos[0] = pu[0]; // A
          pos[1] = pv[1]; // Q
        }
      } else { // s2 == 0
        if (s3 > 0) { // (00+)
          // Q = B.
          types[0] = (int) SHAREVERT;
          pos[0] = pu[1]; // B
          pos[1] = pv[1]; // Q
        } else { // s3 == 0 (000)
          // Impossible.
          assert(0);
        }
      }
    }
  }

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_tri_2d()    Triangle-triangle coplanar intersection test.             //
//                                                                           //
// This routine takes two triangles T1 (with vertices A, B, C) and T2 (with  //
// vertices P, Q, R) in 3D, and T1 and T2 are coplanar, and tests if they    //
// intersect each other.  Return 1 if they intersect, ie., T1 \cap T2 is not //
// empty, otherwise, return 0.                                               //
//                                                                           //
// If 'O' is not NULL, it lies strictly above A, B, C.                       //
//                                                                           //
// If T1 and T2 intersect (return 1) and if 'level' > 0,                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::tri_tri_2d(point A, point B, point C, point P, point Q, 
  point R, point O, int level, int *types, int *pos)
{
  point U[3], V[3];
  int pu[3], pv[3], iv;
  REAL s1, s2, s3, s4;
  REAL s5, s6 ,s7, s8, s9, s10;
  int z1, z2, z3;

  if (O == NULL) {
    REAL n[3], len;
    // Calculate a lift point, saved in dummypoint.
    facenormal(A, B, C, n, 1);
    len = sqrt(DOT(n, n));
    n[0] /= len;
    n[1] /= len;
    n[2] /= len;
    len = DIST(A, B);
    len += DIST(B, C);
    len += DIST(C, A);
    len /= 3.0;
    O = dummypoint;
    O[0] = A[0] + len * n[0];
    O[1] = A[1] + len * n[1];
    O[2] = A[2] + len * n[2];
  }

  s1 = orient3d(A, B, O, P);
  s2 = orient3d(B, C, O, P);
  s3 = orient3d(C, A, O, P);
  orient3dcount+=3;

  if (b->epsilon) {
    if ((s1 != 0) && iscoplanar(A, B, O, P, s1)) s1 = 0;
    if ((s2 != 0) && iscoplanar(B, C, O, P, s2)) s2 = 0;
    if ((s3 != 0) && iscoplanar(C, A, O, P, s3)) s3 = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-tri-2d (%d %d %d)-(%d %d %d) (%c%c%c)\n", pointmark(A),
      pointmark(B), pointmark(C), pointmark(P), pointmark(Q), pointmark(R),
      s1>0 ? '+' : (s1 < 0 ? '-' : '0'), s2>0 ? '+' : (s2<0 ? '-' : '0'),
      s3>0 ? '+' : (s3<0 ? '-' : '0'));
  }

  if (s1 < 0) {
    if (s2 < 0) {
      if (s3 < 0) { // (---)
        assert(0); // Not possible.
      } else {
        if (s3 > 0) { // (--+)
          SETVECTOR3(U, B, C, A);  // PT = ST x ST
          SETVECTOR3(pu, 1, 2, 0);
          z1 = 2;
        } else { // (--0)
          assert(0); // Not possible.
        }
      }
    } else {
      if (s2 > 0) {
        if (s3 < 0) { // (-+-)
          SETVECTOR3(U, A, B, C);  // I3
          SETVECTOR3(pu, 0, 1, 2);
          z1 = 2;
        } else {
          if (s3 > 0) { // (-++)
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(pu, 0, 1, 2);
            z1 = 1;
          } else { // (-+0)
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(pu, 0, 1, 2);
            z1 = 4;
          }
        }
      } else { // s2 == 0
        if (s3 < 0) { //(-0-)
          assert(0); // Not possible
        } else {
          if (s3 > 0) { // (-0+)
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(pu, 1, 2, 0);
            z1 = 3;
          } else { // (-00)
            assert(0); // Not possible
          }
        }
      }
    }
  } else {
    if (s1 > 0) {
      if (s2 < 0) {
        if (s3 < 0) { // (+--)
          SETVECTOR3(U, C, A, B);  // PT = ST
          SETVECTOR3(pu, 2, 0, 1);
          z1 = 2;
        } else {
          if (s3 > 0) { // (+-+)
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(pu, 1, 2, 0);
            z1 = 1;
          } else { // (+-0)
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(pu, 2, 0, 1);
            z1 = 3;
          }
        }
      } else { 
        if (s2 > 0) {
          if (s3 < 0) { // (++-)
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(pu, 2, 0, 1);
            z1 = 1;
          } else {
            if (s3 > 0) { // (+++)
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(pu, 0, 1, 2);
              z1 = 7;
            } else { // (++0)
              SETVECTOR3(U, C, A, B);  // PT = ST
              SETVECTOR3(pu, 2, 0, 1);
              z1 = 5;
            }
          }
        } else { // s2 == 0
          if (s3 < 0) { // (+0-)
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(pu, 1, 2, 0);
            z1 = 4;
          } else {
            if (s3 > 0) { // (+0+)
              SETVECTOR3(U, B, C, A);  // PT = ST x ST
              SETVECTOR3(pu, 1, 2, 0);
              z1 = 5;
            } else { // (+00)
              SETVECTOR3(U, C, A, B);  // PT = ST
              SETVECTOR3(pu, 2, 0, 1);
              z1 = 6;
            }
          }
        }
      }
    } else { // s1 == 0
      if (s2 < 0) {
        if (s3 < 0) { // (0--)
          assert(0); // Not possible
        } else {
          if (s3 > 0) { // (0-+)
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(pu, 1, 2, 0);
            z1 = 4;
          } else { // (0-0)
            assert(0); // Not possible
          }
        }
      } else {
        if (s2 > 0) {
          if (s3 < 0) { // (0+-)
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(pu, 0, 1, 2);
            z1 = 3;
          } else {
            if (s3 > 0) { // (0++)
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(pu, 0, 1, 2);
              z1 = 5;
            } else { // (0+0)
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(pu, 0, 1, 2);
              z1 = 6;
            }
          }
        } else { // s2 == 0
          if (s3 < 0) { // (00-)
            assert(0); // Not possible
          } else {
            if (s3 > 0) { // (00+)
              SETVECTOR3(U, B, C, A);  // PT = ST x ST
              SETVECTOR3(pu, 1, 2, 0);
              z1 = 6;
            } else { // (000)
              assert(0); // Not possible
            }
          }
        }
      }
    }
  }

  if (b->verbose > 2) {
    printf("      Tri-tri-2d (%d %d %d)-(%d) z1(%d)\n", pointmark(U[0]), 
      pointmark(U[1]), pointmark(U[2]), pointmark(P), z1);
  }

  if (level == 0) {
    if (z1 > 4) {
      return 1;  // Intersect.
    }
  }

  // Orient P, Q, R to be CCW wrt O, keep P as origin
  s4 = orient3d(P, Q, R, O);
  orient3dcount++;
  if (b->epsilon) {
    if ((s4 != 0) && iscoplanar(P, Q, R, O, s4)) s4 = 0;
  }

  assert(s4 != 0);
  if (s4 < 0) {
    SETVECTOR3(V, P, Q, R);
    SETVECTOR3(pv, 0, 1, 2);
    iv = 0;
  } else {
    SETVECTOR3(V, P, R, Q);
    SETVECTOR3(pv, 0, 2, 1);
    iv = 1;
  }

  // Decide whether they intersect of not.

  if (z1 == 1) {

    s5 = orient3d(V[0], U[0], O, V[1]); // P, A, O, Q
    s6 = orient3d(V[0], U[1], O, V[1]); // P, B, O, Q
    s7 = orient3d(U[0], U[1], O, V[1]); // A, B, O, Q
    orient3dcount+=3;
    if (b->epsilon) {
      if ((s5 != 0) && iscoplanar(V[0], U[0], O, V[1], s5)) s5 = 0;
      if ((s6 != 0) && iscoplanar(V[0], U[1], O, V[1], s6)) s6 = 0;
      if ((s7 != 0) && iscoplanar(U[0], U[1], O, V[1], s7)) s7 = 0;
    }

    if (s5 < 0) {
      if (s6 < 0) {
        if (s7 < 0) { // (---)
          z2 = 6;
        } else {
          if (s7 > 0) { // (--+)
            z2 = 5;
          } else { // (--0)
            z2 = 12;
          }
        }
      } else {
        if (s6 > 0) {
          if (s7 < 0) { // (-+-)
            // z2 = 1; [P, Q] intersects [A, B, C].
            if (level > 0) {
              types[0] = (int) TRIEDGEINT;
              pos[0] = 3; // [A, B, C]
              pos[1] = pv[0]; // [P, Q]
              types[1] = (int) DISJOINT; // Skipped.
            }
            return 1;
          } else {
            if (s7 > 0) { // (-++)
              z2 = 7;
            } else { // (-+0)
              z2 = 16; assert(0);
            }
          }
        } else { // s6 == 0
          if (s7 < 0) { //(-0-)
            z2 = 13;
          } else {
            if (s7 > 0) { // (-0+)
              z2 = 15;
            } else { // (-00)
              if (level == 0) {
                return 1;
              }
              z2 = 18; // Q = B.
            }
          }
        }
      }
    } else {
      if (s5 > 0) {
        if (s6 < 0) {
          if (s7 < 0) { // (+--)
            assert(0); // Not possible.
          } else {
            if (s7 > 0) { // (+-+)
              z2 = 4;
            } else { // (+-0)
              assert(0); // Not possible.
            }
          }
        } else { 
          if (s6 > 0) {
            if (s7 < 0) { // (++-)
              return 0; // z2 = 2
            } else {
              if (s7 > 0) { // (+++)
                return 0; // z2 = 3
              } else { // (++0)
                return 0; // z2 = 9
              }
            }
          } else { // s6 == 0
            if (s7 < 0) { // (+0-)
              assert(0); // Not possible.
            } else {
              if (s7 > 0) { // (+0+)
                return 0; // z2 = 10
              } else { // (+00)
                assert(0); // Not possible.
              }
            }
          }
        }
      } else { // s5 == 0
        if (s6 < 0) {
          if (s7 < 0) { // (0--)
            assert(0); // Not possible.
          } else {
            if (s7 > 0) { // (0-+)
              z2 = 11;
            } else { // (0-0)
              assert(0); // Not possible.
            }
          }
        } else {
          if (s6 > 0) {
            if (s7 < 0) { // (0+-)
              z2 = 8;
            } else {
              if (s7 > 0) { // (0++)
                return 0; // z2 = 14
              } else { // (0+0)
                // Q = A; // z2 = 17;
                if (level > 0) {
                  types[0] = (int) SHAREVERT;
                  pos[0] = pu[0]; // A
                  pos[1] = pv[1]; // Q
                  types[1] = (int) DISJOINT;
                }
                return 1;
              }
            }
          } else { // s6 == 0
            if (s7 < 0) { // (00-)
              assert(0); // Not possible.
            } else {
              if (s7 > 0) { // (00+)
                z2 = 19; assert(0); // Q = P.
              } else { // (000)
                assert(0); // Not possible.
              }
            }
          }
        }
      }
    }

    s8 = orient3d(V[0], U[0], O, V[2]); // P, A, O, R
    s9 = orient3d(V[0], U[1], O, V[2]); // P, B, O, R
    s10 = orient3d(U[0], U[1], O, V[2]); // A, B, O, R
    orient3dcount+=3;
    if (b->epsilon) {
      if ((s8 != 0) && iscoplanar(V[0], U[0], O, V[2], s8)) s8 = 0;
      if ((s9 != 0) && iscoplanar(V[0], U[1], O, V[2], s9)) s9 = 0;
      if ((s10 != 0) && iscoplanar(U[0], U[1], O, V[2], s10)) s10 = 0;
    }

    if (s8 < 0) {
      if (s9 < 0) {
        if (s10 < 0) { // (---)
          return 0; // z3 = 6;
        } else {
          if (s10 > 0) { // (--+)
            return 0; // z3 = 5;
          } else { // (--0)
            return 0; // z3 = 12;
          }
        }
      } else {
        if (s9 > 0) {
          if (s10 < 0) { // (-+-)
            // [P, R] intersects [A, B, C], z3 = 1;
            if (level > 0) {
              types[0] = (int) TRIEDGEINT;
              pos[0] = 3; // [A, B, C]
              pos[1] = pv[2]; // [R, P]
              types[1] = (int) DISJOINT; // Skipped.
            }
            return 1;
          } else {
            if (s10 > 0) { // (-++)
              z3 = 7;
            } else { // (-+0)
              z3 = 16; assert(0);
            }
          }
        } else { // s9 == 0
          if (s10 < 0) { //(-0-)
            z3 = 13;
          } else {
            if (s10 > 0) { // (-0+)
              return 0; // z3 = 15;
            } else { // (-00)
              // z3 = 18; // R = B.
              if (level > 0) {
                types[0] = (int) SHAREVERT;
                pos[0] = pu[1]; // B
                pos[1] = pv[2]; // R
                types[1] = (int) DISJOINT;
              }
              return 1;
            }
          }
        }
      }
    } else {
      if (s8 > 0) {
        if (s9 < 0) {
          if (s10 < 0) { // (+--)
            assert(0); // Not possible.
          } else {
            if (s10 > 0) { // (+-+)
              return 0; // z3 = 4;
            } else { // (+-0)
              assert(0); // Not possible.
            }
          }
        } else { 
          if (s9 > 0) {
            if (s10 < 0) { // (++-)
              z3 = 2;
            } else {
              if (s10 > 0) { // (+++)
                z3 = 3;
              } else { // (++0)
                z3 = 9;
              }
            }
          } else { // s9 == 0
            if (s10 < 0) { // (+0-)
              assert(0); // Not possible.
            } else {
              if (s10 > 0) { // (+0+)
                z3 = 10;
              } else { // (+00)
                assert(0); // Not possible.
              }
            }
          }
        }
      } else { // s8 == 0
        if (s9 < 0) {
          if (s10 < 0) { // (0--)
            assert(0); // Not possible.
          } else {
            if (s10 > 0) { // (0-+)
              return 0; // z3 = 11;
            } else { // (0-0)
              assert(0); // Not possible.
            }
          }
        } else {
          if (s9 > 0) {
            if (s10 < 0) { // (0+-)
              z3 = 8;
            } else {
              if (s10 > 0) { // (0++)
                z3 = 14;
              } else { // (0+0)
                if (level == 0) {
                  return 1;
                }
                z3 = 17; // R = A;
              }
            }
          } else { // s9 == 0
            if (s10 < 0) { // (00-)
              assert(0); // Not possible.
            } else {
              if (s10 > 0) { // (00+)
                z3 = 19; assert(0); // R = P.
              } else { // (000)
                assert(0); // Not possible.
              }
            }
          }
        }
      }
    }

    if (z2 == 4) {
      return 0;
    }

    if (z2 == 5) {
      if (z3 == 2) {
        s11 = orient3d(V[1], V[2], O, U[0]); // Q, R, O, A
        orient3dcount++;
        if (b->epsilon) {
          if ((s11 != 0) && iscoplanar(V[1], V[2], O, U[0], s11)) s11 = 0;
        }
        if (s11 < 0) {
          return 0;
        } else {
          if (level > 0) {
            if (s11 > 0) {
              // [A, B] intersects [P, Q, R]
              types[0] = (int) EDGETRIINT;
              pos[0] = pu[0]; // [A, B]
              pos[1] = 3; // [P, Q, R]
              types[1] = (int) DISJOINT; // Test skipped.
            } else { // s11 == 0
              // [Q, R] passes A.
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = pv[1]; // [Q, R]
              types[1] = (int) DISJOINT;
            }
          }
          return 1;
        }
      }
      if (z3 == 8) {
        if (level > 0) {
          // [A, B] intersects [P, Q, R]
          types[0] = (int) EDGETRIINT;
          pos[0] = pu[0]; // [A, B]
          pos[1] = 3; // [P, Q, R]
          types[1] = (int) DISJOINT; // Test skipped.
        }
        return 1;
      }
      // Other cases of z3 are all DISJOINT.
      return 0;
    }

    if (z2 == 6) {
    
    }

  } // z1 == 1

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_tri_test()    Triangle-triangle intersection test.                    //
//                                                                           //
// This routine takes two triangles T1 (with vertices A, B, C) and T2 (P, Q  //
// R) in 3D, and tests if they intersect each other.  Return 1 if they are   //
// intersected, i.e., T1 \cap T2 is not empty, otherwise, return 0.          //
//                                                                           //
// If the point 'O' is not NULL, it lies strictly above the plane defined by //
// A, B, C. It is used in test when T1 and T2 are coplanar.                  //
//                                                                           //
// If T1 and T2 intersect each other (return 1), they may intersect in diff- //
// erent ways. If 'level' > 0, their intersection type will be reported in   //
// combinations of 'types' and 'pos'.                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::tri_tri_test(point A, point B, point C, point P, point Q, 
  point R, point O, int level, int *types, int *pos)
{
  point U[3], V[3], W[3], Ptmp;  // The permuted vectors of points.
  int pu[3], pv[3], pw[3], iu, iv, itmp; 
  REAL sA, sB, sC, sP, sQ, sR;
  REAL s1, s2, s3, s4;
  int z1, z2;

  // Test A's, B's, and C's orientations wrt plane PQR. 
  sA = orient3d(P, Q, R, A);
  sB = orient3d(P, Q, R, B);
  sC = orient3d(P, Q, R, C);
  orient3dcount+=3;

  if (b->epsilon) {
    if ((sA != 0) && iscoplanar(P, Q, R, A, sA)) sA = 0;
    if ((sB != 0) && iscoplanar(P, Q, R, B, sB)) sB = 0;
    if ((sC != 0) && iscoplanar(P, Q, R, C, sC)) sC = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-tri (%d %d %d)-(%d %d %d) (%c%c%c)\n", pointmark(A),
      pointmark(B), pointmark(C), pointmark(P), pointmark(Q), pointmark(R),
      sA > 0 ? '+' : (sA < 0 ? '-' : '0'), sB>0 ? '+' : (sB<0 ? '-' : '0'),
      sC>0 ? '+' : (sC<0 ? '-' : '0'));
  }
  // tritricount++;
  iu = iv = 0;

  if (sA < 0) {
    if (sB < 0) {
      if (sC < 0) { // (---).
        return DISJOINT; 
      } else {
        if (sC > 0) { // (--+).
          // All points are in the right positions.
          SETVECTOR3(U, A, B, C);  // I3
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 0, 1, 2);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else { // (--0).
          SETVECTOR3(U, A, B, C);  // I3
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 0, 1, 2);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        }
      }
    } else { 
      if (sB > 0) {
        if (sC < 0) { // (-+-).
          SETVECTOR3(U, C, A, B);  // PT = ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 2, 0, 1);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else {
          if (sC > 0) { // (-++).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 0;
          } else { // (-+0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          }
        }
      } else {
        if (sC < 0) { // (-0-).
          SETVECTOR3(U, C, A, B);  // PT = ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 2, 0, 1);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        } else {
          if (sC > 0) { // (-0+).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 2;
          } else { // (-00).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 3;
          }
        }
      }
    }
  } else {
    if (sA > 0) {
      if (sB < 0) {
        if (sC < 0) { // (+--).
          SETVECTOR3(U, B, C, A);  // PT = ST x ST
          SETVECTOR3(V, P, Q, R);  // I2
          SETVECTOR3(pu, 1, 2, 0);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 0;
        } else {
          if (sC > 0) { // (+-+).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 0;
          } else { // (+-0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 2;
          }
        }
      } else { 
        if (sB > 0) {
          if (sC < 0) { // (++-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 0;
          } else {
            if (sC > 0) { // (+++).
              return DISJOINT; 
            } else { // (++0).
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 0, 1, 2);
              SETVECTOR3(pv, 1, 0, 2); iv = 1;
              z1 = 1; 
            }
          }
        } else {
          if (sC < 0) { // (+0-).
            SETVECTOR3(U, B, C, A);  // PT = ST x ST
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 1, 2, 0);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          } else {
            if (sC > 0) { // (+0+).
              SETVECTOR3(U, C, A, B);  // PT = ST
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 2, 0, 1);
              SETVECTOR3(pv, 1, 0, 2); iv = 1;
              z1 = 1;
            } else { // (+00).
              SETVECTOR3(U, B, C, A);  // PT = ST x ST
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 1, 2, 0);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            }
          }
        }
      }
    } else { 
      if (sB < 0) {
        if (sC < 0) { // (0--).
          SETVECTOR3(U, B, C, A);  // PT = ST x ST
          SETVECTOR3(V, P, Q, R);  // PL = I2
          SETVECTOR3(pu, 1, 2, 0);
          SETVECTOR3(pv, 0, 1, 2);
          z1 = 1;
        } else {
          if (sC > 0) { // (0-+).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, P, Q, R);  // I2
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 0, 1, 2);
            z1 = 2;
          } else { // (0-0).
            SETVECTOR3(U, C, A, B);  // PT = ST
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 2, 0, 1);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 3; 
          }
        }
      } else { 
        if (sB > 0) {
          if (sC < 0) { // (0+-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 2;
          } else {
            if (sC > 0) { // (0++).
              SETVECTOR3(U, B, C, A);  // PT = ST x ST
              SETVECTOR3(V, Q, P, R);  // PL = SL
              SETVECTOR3(pu, 1, 2, 0);
              SETVECTOR3(pv, 1, 0, 2); iv = 1;
              z1 = 1;
            } else { // (0+0).
              SETVECTOR3(U, C, A, B);  // PT = ST
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 2, 0, 1);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            }
          }
        } else {
          if (sC < 0) { // (00-).
            SETVECTOR3(U, A, B, C);  // I3
            SETVECTOR3(V, Q, P, R);  // PL = SL
            SETVECTOR3(pu, 0, 1, 2);
            SETVECTOR3(pv, 1, 0, 2); iv = 1;
            z1 = 3; 
          } else {
            if (sC > 0) { // (00+).
              SETVECTOR3(U, A, B, C);  // I3
              SETVECTOR3(V, P, Q, R);  // I2
              SETVECTOR3(pu, 0, 1, 2);
              SETVECTOR3(pv, 0, 1, 2);
              z1 = 3; 
            } else { // (000)
              // (A, B, C) is coplanar with (P, Q, R).
              z1 = 4;
            }
          }
        }
      }
    }
  }

  sP = orient3d(U[0], U[1], U[2], V[0]);
  sQ = orient3d(U[0], U[1], U[2], V[1]);
  sR = orient3d(U[0], U[1], U[2], V[2]);
  orient3dcount+=3;

  if (b->epsilon) {
    if ((sP != 0) && iscoplanar(U[0], U[1], U[2], V[0], sP)) sP = 0;
    if ((sQ != 0) && iscoplanar(U[0], U[1], U[2], V[1], sQ)) sQ = 0;
    if ((sR != 0) && iscoplanar(U[0], U[1], U[2], V[2], sR)) sR = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-tri (%d %d %d)-(%d %d %d) (%c%c%c)\n", pointmark(U[0]),
      pointmark(U[1]), pointmark(U[2]), pointmark(V[0]), pointmark(V[1]), 
      pointmark(V[2]), sP>0 ? '+' : (sP<0 ? '-' : '0'), 
      sQ>0 ? '+' : (sQ<0 ? '-' : '0'), sR>0 ? '+' : (sR<0 ? '-' : '0'));
  }

  if (sP < 0) {
    if (sQ < 0) {
      if (sR < 0) {  // (---)
        return DISJOINT;
      } else {
        if (sR > 0) { // (--+)
          // P1->Q1 is opposite to A1->B1. Swicth A1 and B1.
          SETVECTOR3(W, V[0], V[1], V[2]);  // I3
          SWAP2(U[0], U[1], Ptmp);  // PL = SL
          SETVECTOR3(pw, pv[0], pv[1], pv[2]);
          SWAP2(pu[0], pu[1], itmp); iu = 1;
          z2 = 0;
        } else {  // (--0)
          SETVECTOR3(W, V[0], V[1], V[2]);  // I3
          SWAP2(U[0], U[1], Ptmp);  // PL = SL
          SETVECTOR3(pw, pv[0], pv[1], pv[2]);
          SWAP2(pu[0], pu[1], itmp); iu = 1;
          z2 = 1;
        }
      }
    } else {
      if (sQ > 0) {
        if (sR < 0) {  // (-+-)
          SETVECTOR3(W, V[2], V[0], V[1]);  // PT = ST
          SWAP2(U[0], U[1], Ptmp);  // PL = SL
          SETVECTOR3(pw, pv[2], pv[0], pv[1]);
          SWAP2(pu[0], pu[1], itmp); iu = 1;
          z2 = 0;
        } else {
          if (sR > 0) { // (-++)
            SETVECTOR3(W, V[1], V[2], V[0]);  // PT = ST x ST
            SETVECTOR3(pw, pv[1], pv[2], pv[0]);
            z2 = 0;
          } else { // (-+0)
            SETVECTOR3(W, V[2], V[0], V[1]);  // PT = ST
            SWAP2(U[0], U[1], Ptmp); // PL = SL
            SETVECTOR3(pw, pv[2], pv[0], pv[1]);
            SWAP2(pu[0], pu[1], itmp); iu = 1;
            z2 = 2;
          }
        }
      } else {
        if (sR < 0) {  // (-0-)
          SETVECTOR3(W, V[2], V[0], V[1]);  // PT = ST
          SWAP2(U[0], U[1], Ptmp); // PL = SL
          SETVECTOR3(pw, pv[2], pv[0], pv[1]);
          SWAP2(pu[0], pu[1], itmp); iu = 1;
          z2 = 1;
        } else {
          if (sR > 0) {  // (-0+)
            SETVECTOR3(W, V[1], V[2], V[0]);  // PT = ST x ST
            SETVECTOR3(pw, pv[1], pv[2], pv[0]);
            z2 = 2;
          } else {  // (-00)
            SETVECTOR3(W, V[1], V[2], V[0]);  // PT = ST x ST
            SETVECTOR3(pw, pv[1], pv[2], pv[0]);
            z2 = 3;
          }
        }
      }
    }
  } else {
    if (sP > 0) {
      if (sQ < 0) {
        if (sR < 0) {  // (+--)
          SETVECTOR3(W, V[1], V[2], V[0]);  // PT = ST x ST
          SWAP2(U[0], U[1], Ptmp); // PL = SL
          SETVECTOR3(pw, pv[1], pv[2], pv[0]);
          SWAP2(pu[0], pu[1], itmp); iu = 1;
          z2 = 0;
        } else {
          if (sR > 0) {  // (+-+)
            SETVECTOR3(W, V[2], V[0], V[1]);  // PT = ST
            SETVECTOR3(pw, pv[2], pv[0], pv[1]);
            z2 = 0;
          } else {  // (+-0)
            SETVECTOR3(W, V[2], V[0], V[1]);  // PT = ST
            SETVECTOR3(pw, pv[2], pv[0], pv[1]);
            z2 = 2;
          }
        }
      } else {
        if (sQ > 0) {
          if (sR < 0) {  // (++-)
            SETVECTOR3(W, V[0], V[1], V[2]);  // I3
            SETVECTOR3(pw, pv[0], pv[1], pv[2]);
            z2 = 0;
          } else {
            if (sR > 0) {  // (+++)
              return DISJOINT;
            } else {  // (++0)
              SETVECTOR3(W, V[0], V[1], V[2]);  // I3
              SETVECTOR3(pw, pv[0], pv[1], pv[2]);
              z2 = 1;
            }
          }
        } else { // sQ == 0
          if (sR < 0) {  // (+0-)
            SETVECTOR3(W, V[1], V[2], V[0]);  // PT = ST x ST
            SWAP2(U[0], U[1], Ptmp); // PL = SL
            SETVECTOR3(pw, pv[1], pv[2], pv[0]);
            SWAP2(pu[0], pu[1], itmp); iu = 1;
            z2 = 2;
          } else {
            if (sR > 0) {  // (+0+)
              SETVECTOR3(W, V[2], V[0], V[1]);  // PT = ST
              SETVECTOR3(pw, pv[2], pv[0], pv[1]);
              z2 = 1;
            } else {  // (+00)
              SETVECTOR3(W, V[1], V[2], V[0]);  // PT = ST x ST
              SWAP2(U[0], U[1], Ptmp); // PL = SL
              SETVECTOR3(pw, pv[1], pv[2], pv[0]);
              SWAP2(pu[0], pu[1], itmp); iu = 1;
              z2 = 3;
            }
          }
        }
      }
    } else { // sP == 0
      if (sQ < 0) {
        if (sR < 0) {  // (0--)
          SETVECTOR3(W, V[1], V[2], V[0]);  // PT = ST x ST
          SWAP2(U[0], U[1], Ptmp); // PL = SL
          SETVECTOR3(pw, pv[1], pv[2], pv[0]);
          SWAP2(pu[0], pu[1], itmp); iu = 1;
          z2 = 1;
        } else {
          if (sR > 0) {  // (0-+)
            SETVECTOR3(W, V[0], V[1], V[2]);  // I3
            SWAP2(U[0], U[1], Ptmp); // PL = SL
            SETVECTOR3(pw, pv[0], pv[1], pv[2]);
            SWAP2(pu[0], pu[1], itmp); iu = 1;
            z2 = 2;
          } else {  // (0-0)
            SETVECTOR3(W, V[2], V[0], V[1]);  // PT = ST
            SETVECTOR3(pw, pv[2], pv[0], pv[1]);
            z2 = 3;
          }
        }
      } else {
        if (sQ > 0) {
          if (sR < 0) {  // (0+-)
            SETVECTOR3(W, V[0], V[1], V[2]);  // I3
            SETVECTOR3(pw, pv[0], pv[1], pv[2]);
            z2 = 2;
          } else {
            if (sR > 0) {  // (0++)
              SETVECTOR3(W, V[1], V[2], V[0]);  // PT = ST x ST
              SETVECTOR3(pw, pv[1], pv[2], pv[0]);
              z2 = 1;
            } else {  // (0+0)
              SETVECTOR3(W, V[2], V[0], V[1]);  // PT = ST
              SWAP2(U[0], U[1], Ptmp); // PL = SL
              SETVECTOR3(pw, pv[2], pv[0], pv[1]);
              SWAP2(pu[0], pu[1], itmp); iu = 1;
              z2 = 3;
            }
          }
        } else { // sQ == 0
          if (sR < 0) {  // (00-)
            SETVECTOR3(W, V[0], V[1], V[2]);  // I3
            SETVECTOR3(pw, pv[0], pv[1], pv[2]);
            z2 = 3;
          } else {
            if (sR > 0) {  // (00+)
              SETVECTOR3(W, V[0], V[1], V[2]);  // I3
              SWAP2(U[0], U[1], Ptmp); // PL = SL
              SETVECTOR3(pw, pv[0], pv[1], pv[2]);
              SWAP2(pu[0], pu[1], itmp); iu = 1;
              z2 = 3;
            } else {  // (000)
              z2 = 4;
            }
          }
        }
      }
    }
  }

  if (z1 == 4) {
    assert(z2 == 4);  // SELF_CHECK
    assert(0); // Not implement yet.
    // return tri_tri_2d();
  }

  if (z2 == 1) {
    s1 = orient3d(U[0], U[2], W[0], W[2]);  // A, C, P, R
    s2 = orient3d(U[1], U[2], W[1], W[2]);  // B, C, Q, R
    orient3dcount+=2;
    if (b->epsilon > 0) {
      if ((s1 != 0) && iscoplanar(U[0], U[2], W[0], W[2], s1)) s1 = 0;
      if ((s2 != 0) && iscoplanar(U[1], U[2], W[1], W[2], s2)) s2 = 0;
    }
  } else {
    s1 = orient3d(U[0], U[2], W[2], W[1]);  // A, C, R, Q
    s2 = orient3d(U[1], U[2], W[2], W[0]);  // B, C, R, P
    orient3dcount+=2;
    if (b->epsilon > 0) {
      if ((s1 != 0) && iscoplanar(U[0], U[2], W[2], W[1], s1)) s1 = 0;
      if ((s2 != 0) && iscoplanar(U[1], U[2], W[2], W[0], s2)) s2 = 0;
    }
  }

  if (b->verbose > 2) {
    printf("      Tri-tri (%d %d %d)-(%d %d %d) (%d-%d) (%c%c)\n",
      pointmark(U[0]), pointmark(U[1]), pointmark(U[2]), pointmark(W[0]), 
      pointmark(W[1]), pointmark(W[2]), z1, z2,
      s1>0 ? '+' : (s1<0 ? '-' : '0'), s2>0 ? '+' : (s2<0 ? '-' : '0'));
  }

  if (z2 == 1) {
    if (s1 < 0) {
      return 0;
    }
    if (s2 > 0) {
      return 0;
    }
  } else {
    if (s1 > 0) {
      return 0;
    }
    if (s2 < 0) {
      return 0;
    }
  }

  // Do we need to classify the intersection type of T1 and T2?
  if (level == 0) {
    return 1;
  }

  if (z2 == 1) {

    if (z1 == 0) {  // (01)
      if (s1 == 0) {
        // R = k in [A, C] (tritri-010###).
        types[0] = (int) TOUCHEDGE;
        pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // [C, A]
        pos[1] = pw[2]; // R
        types[1] = (int) DISJOINT; // No 2nd inter-point.
      } else {
        if (s2 == 0) {
          // R = l in [B, C] (tritri-01#0##).
          types[0] = (int) TOUCHEDGE;
          pos[0] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // [B, C]
          pos[1] = pw[2]; // R
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        } else {
          // R in [k, l] in [A, B, C] (tritri-01+-##).
          types[0] = (int) TOUCHFACE;
          pos[0] = 3; // Interior of [A, B, C]
          pos[1] = pw[2]; // R
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        }
      }
    } else if (z1 == 1) {  // (11)
      assert(s1 == 0); // SELF_CHECK
      // C = R (tritri-110###).
      types[0] = (int) SHAREVERT;
      pos[0] = pu[2]; // C
      pos[1] = pw[2]; // R
      types[1] = (int) DISJOINT; // No 2nd inter-point.
    } else if (z1 == 2) {  // (21)
      if (s1 == 0) {
        // R = A (tritri-210###).
        types[0] = (int) SHAREVERT;
        pos[0] = pu[0]; // A
        pos[1] = pw[2]; // R
        types[1] = (int) DISJOINT; // No 2nd inter-point.
      } else {
        if (s2 == 0) {
          // R = l, R in [B, C] (tritri-21#0##).
          types[0] = (int) TOUCHEDGE;
          pos[0] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // [B, C]
          pos[1] = pw[2]; // R
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        } else {
          // R in [k, l] in [A, B, C] (tritri-21+-##).
          types[0] = (int) TOUCHFACE;
          pos[0] = 3; // Interior of [A, B, C]
          pos[1] = pw[2]; // R
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        }
      }
    } else if (z1 == 3) {  // (31)
      if (s1 == 0) {
        // R = A (tritri-310###).
        types[0] = (int) SHAREVERT;
        pos[0] = pu[0]; // A
        pos[1] = pw[2]; // R
        types[1] = (int) DISJOINT; // No 2nd inter-point.
      } else {
        if (s2 == 0) {
          // R = B (tritri-31#0##).
          types[0] = (int) SHAREVERT;
          pos[0] = pu[1]; // B
          pos[1] = pw[2]; // R
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        } else {
          // R in [A, B] (tritri-31+-##).
          types[0] = (int) TOUCHEDGE;
          pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // [A, B]
          pos[1] = pw[2]; // R
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        }
      }
    }

    return 1;
  } // z2 == 1

  if (z1 == 1) {
  
    if (z2 == 0) { // (10)
      if (s1 == 0) {
        // C = j, C in [Q, R] (tritri-100###).
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[2]; // C
        pos[1] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // [Q, R]
        types[1] = (int) DISJOINT; // No 2nd inter-point.
      } else {
        if (s2 == 0) {
          // C = i, C in [P, R] (tritri-10#0##).
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[2]; // C
          pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // [R, P]
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        } else {
          // C in [i, j] in [P, Q, R] (tritri-10-+##).
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[2]; // C
          pos[1] = 3; // Interior of [P, Q, R]
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        }
      }
    } else if (z2 == 2) {  // (12)
      if (s1 == 0) { // 
        // C = j, C in [Q, R] (tritri-120###).
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[2]; // C
        pos[1] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // [Q, R]
        types[1] = (int) DISJOINT; // No 2nd inter-point.
      } else {
        if (s2 == 0) {
          // C = P (tritri-12#0##).
          types[0] = (int) SHAREVERT;
          pos[0] = pu[2]; // C
          pos[1] = pw[0]; // P
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        } else {
          // C in [i, j] in [P, Q, R] (tritri-12-+##).
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[2]; // C
          pos[1] = 3; // Interior of [P, Q, R]
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        }
      }
    } else if (z2 == 3) {  // (13)
      if (s1 == 0) {
        // C = Q (tritri-130###).
        types[0] = (int) SHAREVERT;
        pos[0] = pu[2]; // C
        pos[1] = pw[1]; // Q
        types[1] = (int) DISJOINT; // No 2nd inter-point.
      } else {
        if (s2 == 0) {
          // C = P (tritri-13#0##).
          types[0] = (int) SHAREVERT;
          pos[0] = pu[2]; // C
          pos[1] = pw[0]; // P
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        } else {
          // C in [P, Q] (tritri-13-+##).
          types[0] = (int) ACROSSVERT;
          pos[0] = pu[2]; // C
          pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // [P, Q]
          types[1] = (int) DISJOINT; // No 2nd inter-point.
        }
      }
    }

    return 1;
  } // if (z1 == 1)

  // Two additional orientation tests are needed.
  s3 = orient3d(U[0], U[2], W[2], W[0]);  // A, C, R, P
  s4 = orient3d(U[1], U[2], W[2], W[1]);  // B, C, R, Q
  orient3dcount+=2;
  if (b->epsilon > 0) {
    if ((s3 != 0) && iscoplanar(U[0], U[2], W[2], W[0], s3)) s3 = 0;
    if ((s4 != 0) && iscoplanar(U[1], U[2], W[2], W[1], s4)) s4 = 0;
  }

  if (b->verbose > 2) {
    printf("      (tritri-%d%d%c%c%c%c)\n", z1, z2, 
      s1>0 ? '+' : (s1<0 ? '-' : '0'), s2>0 ? '+' : (s2<0 ? '-' : '0'), 
      s3>0 ? '+' : (s3<0 ? '-' : '0'), s4>0 ? '+' : (s4<0 ? '-' : '0'));
  }

  if (z1 == 0) {
    
    if (z2 == 0) {  // (00)
      if (s1 < 0) {
        if (s3 > 0) {
          assert(s2 > 0); // SELF_CHECK.
          if (s4 > 0) {
            // [i, j] overlaps [k, l] (tritri-00-+++).
            types[0] = (int) ACROSSEDGE;
            pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
            pos[1] = 3;     // Int([P, Q, R])
            types[1] = (int) ACROSSFACE;
            pos[2] = 3;     // Int([A, B, C])
            pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
          } else {
            if (s4 == 0) {
              // j = l, [i, j] contains [k, l] (tritri-00-++0).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else {
              // [i, j] contains [k, l] (tritri-00-++-).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = 3;     // Int([P, Q, R])
            }
          }
        } else {
          if (s3 == 0) { // (00-+0#)
            assert(s2 > 0); // SELF_CHECK.
            if (s4 > 0) {
              // i = k, [i, j] overlaps [k, l] (tritri-00-+0+).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
              types[1] = (int) ACROSSFACE;
              pos[2] = 3;     // Int([A, B, C])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else {
              if (s4 == 0) {
                // [i, j] = [k, l] (tritri-00-+00).
                types[0] = (int) ACROSSEDGE;
                pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else {
                // i = k, [i, j] contains [k, l] (tritri-00-+0-).
                types[0] = (int) ACROSSEDGE;
                pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = 3;     // Int([P, Q, R]) 
              }
            }
          } else { // s3 < 0
            if (s2 > 0) {
              if (s4 > 0) {
                // [i, j] in [k, l] in [A, B, C] (tritri-00-+-+).
                types[0] = (int) ACROSSFACE;
                pos[0] = 3;     // Int([A, B, C])
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSFACE;
                pos[2] = 3;     // Int([A, B, C])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else {
                if (s4 == 0) {
                  // j = l, [i, j] in [k, l] (tritri-00-+-0)
                  types[0] = (int) ACROSSFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
                } else { // s4 < 0
                  // [i, j] overlaps [k, l] (tritri-00-+--)
                  types[0] = (int) ACROSSFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = 3;     // Int([P, Q, R])
                }
              }
            } else { // s2 == 0
              assert(s4 < 0); // SELF_CHECK
              // i = l, [P, R] intersects [B, C] (tritri-00#0##).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
              types[1] = (int) DISJOINT; // No 2nd intersection.
            }
          }
        }
      } else { // s1 == 0
        // j = k, [Q, R] intersects [A, B] (tritri-000###).
        types[0] = (int) ACROSSEDGE;
        pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
        pos[1] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
        types[1] = (int) DISJOINT; // No 2nd intersection.
      }
    } else if (z2 == 2) {  // (02)
      if (s1 < 0) {
        if (s3 > 0) { // (02-#+#)
          assert(s2 > 0); // SELF_CHECK;
          if (s4 > 0) {
            // [P, j] overlaps [k, l] (tritri-02-+++).
            types[0] = (int) ACROSSEDGE;
            pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
            pos[1] = 3;     // Int([P, Q, R])
            types[1] = (int) ACROSSFACE;
            pos[2] = 3;     // Int([A, B, C])
            pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
          } else {
            if (s4 == 0) {
              // j = k, [i, j] contains [k, l] (tritri-02-++0).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else { // s4 < 0
              // [i, j] contains [k, l] (tritri-02-++-).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = 3;     // Int([P, Q, R])
            }
          }
        } else {
          if (s3 == 0) { // (02-#0#)
            assert(s2 > 0); // SELF_CHECK
            if (s4 > 0) {
              // P = k, [P, j] in [k, l] (tritri-02-+0+).
              types[0] = (int) TOUCHEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = pv[0]; // P
              types[1] = (int) ACROSSFACE;
              pos[2] = 3;     // Int([A, B, C])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else {
              if (s4 == 0) {
                // [P, j] = [k, l] (tritri-02-+00).
                types[0] = (int) TOUCHEDGE;
                pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
                pos[1] = pv[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else {  // s4 < 0
                // P = k, [P, j] contains [k, l] (tritri-02-+0-).
                types[0] = (int) TOUCHEDGE;
                pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
                pos[1] = pv[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = 3;     // Int([P, Q, R])
              }
            }
          } else { // s3 < 0.
            if (s2 > 0) {
              if (s4 > 0) {
                // [P, j] in [k, l] in [A, B, C] (tritri-02-+-+).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // Int([A, B, C])
                pos[1] = pv[0]; // P
                types[1] = (int) ACROSSFACE;
                pos[2] = 3;     // Int([A, B, C])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else {
                if (s4 == 0) {
                  // j = l, [P, j] overlaps [k, l] (tritri-02-+-0).
                  types[0] = (int) TOUCHFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = pv[0]; // P
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
                } else { // s4 < 0
                  // [P, j] overlaps [k, l] (tritri-02-+--).
                  types[0] = (int) TOUCHFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = pv[0]; // P
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = 3;     // Int([P, Q, R])
                }
              }
            } else {  // s2 == 0
              assert(s4 < 0); // SELF_CHECK
              // P = k, P in [B, C] (tritri-02#0##).
              types[0] = (int) TOUCHEDGE;
              pos[0] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[1] = pw[0]; // P
              types[1] = (int) DISJOINT;  // No 2nd intersection
            }
          }
        }
      } else { // s1 == 0
        // j = k, [Q, R] intersects [A, C] (tritri-020###).
        types[0] = (int) ACROSSEDGE;
        pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
        pos[1] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
        types[1] = (int) DISJOINT;  // No 2nd intersection
      }
    } else if (z2 == 3) {  // (03)
      if (s1 < 0) {
        if (s3 > 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // [P, Q] overlaps [k, l] (tritri-03-+++).
            types[0] = (int) ACROSSEDGE;
            pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
            pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
            types[1] = (int) TOUCHFACE;
            pos[2] = 3;     // Int([A, B, C])
            pos[3] = pw[1]; // Q
          } else {
            if (s4 == 0) {
              // Q = l, [P, Q] contains [k, l] (tritri-03-++0).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              types[1] = (int) TOUCHEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = pw[1]; // Q
            } else { // s4 < 0
              // [P, Q] contains [k, l] (tritri-03-++-).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
            }
          }
        } else {
          if (s3 == 0) {
            assert(s2 > 0);  // SELF_CHECK
            if (s4 > 0) {
              // P = k, [P, Q] in [k, l] (tritri-03-+0+).
              types[0] = (int) TOUCHEDGE;
              pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
              pos[1] = pw[0]; // P
              types[1] = (int) TOUCHFACE;
              pos[2] = 3;     // Int([A, B, C])
              pos[3] = pw[1]; // Q
            } else {
              if (s4 == 0) {
                // [P, Q] = [k, l] (tritri-03-+00).
                types[0] = (int) TOUCHEDGE;
                pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
                pos[1] = pw[0]; // P
                types[1] = (int) TOUCHEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = pw[1]; // Q
              } else {  // s4 < 0
                // P = k, [P, Q] contains [k, l] (tritri-03-+0-).
                types[0] = (int) TOUCHEDGE;
                pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              }
            }
          } else { // s3 < 0
            if (s2 > 0) {
              if (s4 > 0) {
                // [P, Q] in [k, l] (tritri-03-+-+).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // Int([A, B, C])
                pos[1] = pw[0]; // P
                types[1] = (int) TOUCHFACE;
                pos[2] = 3;     // Int([A, B, C])
                pos[3] = pw[1]; // Q
              } else {
                if (s4 == 0) {
                  // Q = l, [P, Q] in [k, l] (tritri-03-+-0).
                  types[0] = (int) TOUCHFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = pw[0]; // P
                  types[1] = (int) TOUCHEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = pw[1]; // Q
                } else { // s4 < 0
                  // [P, Q] overlaps [k, l] (tritri-03-+--).
                  types[0] = (int) TOUCHFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = pw[0]; // P
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
                }
              }
            } else {  // s2 == 0
              // P = l in [B, C] (tritri-03#0##).
              types[0] = (int) TOUCHEDGE;
              pos[0] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[1] = pw[0]; // P
              types[1] = (int) DISJOINT; // No 2nd intersection
            }
          }
        }
      } else {
        // Q = k in [A, C] (tritri-030###).
        types[0] = (int) TOUCHEDGE;
        pos[0] = (iu == 0 ? pu[2] : mi1mo3[pu[2]]); // Int([C, A])
        pos[1] = pw[1]; // Q
        types[1] = (int) DISJOINT; // No 2nd intersection
      }
    }

    return 1;
  } // if (z1 == 0)

  if (z1 == 2) {

    if (z2 == 0) {  // (20)
      if (s1 < 0) {
        if (s3 > 0) {
          assert(s2 > 0);  // SELF_CHECK
          if (s4 > 0) {
            // [i, j] overlaps [A, l] (tritri-20-+++).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = 3;     // Int([P, Q, R])
            types[1] = (int) ACROSSFACE;
            pos[2] = 3;     // Int([A, B, C])
            pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
          } else {
            if (s4 == 0) {
              // j = l, [i, j] contains [A, l] (tritri-20-++0).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else {  // s4 < 0
              // [i, j] contains [A, l] (tritri-20-++-).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = 3;     // Int([P, Q, R])
            }
          }
        } else {
          if (s3 == 0) {
            assert(s2 > 0); // SELF_CHECK
            if (s4 > 0) {
              // i = A, [i, j] overlaps [A, l] (tritri-20-+0+).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
              types[1] = (int) ACROSSFACE;
              pos[2] = 3;     // Int([A, B, C])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else {
              if (s4 == 0) {
                // [i, j] = [A, l] (tritri-20-+00).
                types[0] = (int) ACROSSVERT;
                pos[0] = pu[0]; // A
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else { // s4 < 0
                // i = A, [i, j] contains [A, l] (tritri-20-+0-).
                types[0] = (int) ACROSSVERT;
                pos[0] = pu[0]; // A
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = 3;     // Int([P, Q, R])
              }
            }
          } else { // s3 < 0
            if (s2 > 0) {
              if (s4 > 0) {
                // [i, j] in [A, l] in [A, B, C] (tritri-20-+-+).
                types[0] = (int) ACROSSFACE;
                pos[0] = 3;     // Int([A, B, C])
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSFACE;
                pos[2] = 3;     // Int([A, B, C])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else {
                if (s4 == 0) {
                  // j = l, [i, j] in [A, l] (tritri-20-+-0).
                  types[0] = (int) ACROSSFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
                } else { // s4 < 0
                  // [i, j] overlaps [A, l] (tritri-20-+--).
                  types[0] = (int) ACROSSFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = 3;     // Int([P, Q, R])
                }
              }
            } else { // s2 == 0
              // i = l, [P, R] intersects [B, C] (tritri-20#0##).
              types[0] = (int) ACROSSEDGE;
              pos[0] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
              types[1] = (int) DISJOINT; // No 2nd intersection
            }
          }
        }
      } else { // s1 == 0
        // j = A in [Q, R] (tritri-200###).
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[0]; // A
        pos[1] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
        types[1] = (int) DISJOINT; // No 2nd intersection
      }
    } else if (z2 == 2) {  // (22)
      if (s1 < 0) {
        if (s3 > 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // [P, j] overlaps [A, l] (tritri-22-+++).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = 3;     // Int([P, Q, R])
            types[1] = (int) ACROSSFACE;
            pos[2] = 3;     // Int([A, B, C])
            pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
          } else {
            if (s4 == 0) {
              // j = l, [P, j] contains [A, l] (tritri-22-++0).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else { // s4 < 0
              // [P, j] contains [A, l] (tritri-22-++-).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = 3;     // Int([P, Q, R])
            }
          }
        } else {
          if (s3 == 0) {
            assert(s2 > 0); // SELF_CHECK
            if (s4 > 0) {
              // P = A, [P, j] in [A, l] (tritri-22-+0+).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pw[0]; // P
              types[1] = (int) ACROSSFACE;
              pos[2] = 3;     // Int([A, B, C])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else {
              if (s4 == 0) {
                // [P, j] = [A, l] (tritri-22-+00).
                types[0] = (int) SHAREVERT;
                pos[0] = pu[0]; // A
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else { // s4 < 0
                // P = A, [P, j] contains [A, l] (tritri-22-+0-).
                types[0] = (int) SHAREVERT;
                pos[0] = pu[0]; // A
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = 3;     // Int([P, Q, R])
              }
            }
          } else { // s3 < 0
            if (s2 > 0) {
              if (s4 > 0) {
                // [P, j] in [A, l] in [A, B, C] (tritri-22-+-+).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // Int([A, B, C])
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSFACE;
                pos[2] = 3;     // Int([A, B, C])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else {
                if (s4 == 0) {
                  // j = l, [P, j] in [A, l] (tritri-22-+-0).
                  types[0] = (int) TOUCHFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = pw[0]; // P
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
                } else { // s4 < 0
                  // P, j] overlaps [A, l] (tritri-22-+--).
                  types[0] = (int) TOUCHFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = pw[0]; // P
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = 3;     // Int([P, Q, R])
                }
              }
            } else { // s2 == 0
              // P = l, P in [B, C] (tritri-22#0##).
              types[0] = (int) TOUCHEDGE;
              pos[0] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int({[B, C])
              pos[1] = pw[0]; // P
              types[1] = (int) DISJOINT; // No 2nd intersection
            }
          }
        }
      } else {
        // j = A in [Q, R] (tritri-220###).
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[0]; // A
        pos[1] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
        types[1] = (int) DISJOINT; // No 2nd intersection
      }
    } else if (z2 == 3) {  // (23)
      if (s1 < 0) {
        if (s3 > 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // [P, Q] overlaps [A, l] (tritri-23-+++).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
            types[1] = (int) TOUCHFACE;
            pos[2] = 3;     // Int([A, B, C])
            pos[3] = pw[1]; // Q
          } else {
            if (s4 == 0) {
              // Q = l, [P, Q] contains [A, l] (tritri-23-++0).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              types[1] = (int) TOUCHEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = pw[1]; // Q
            } else { // s4 < 0
              // [P, Q] contains [A, l] (tritri-23-++-).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
            }
          }
        } else {
          if (s3 == 0) {
            assert(s2 > 0); // SELF_CHECK
            if (s4 > 0) {
              // P = A, [P, Q] in [A, l] (tritri-23-+0+).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pw[0]; // P
              types[1] = (int) TOUCHFACE;
              pos[2] = 3;     // Int([A, B, C])
              pos[3] = pw[1]; // Q
            } else {
              if (s4 == 0) {
                // [P, Q] = [A, l] (tritri-23-+00).
                types[0] = (int) SHAREVERT;
                pos[0] = pu[0]; // A
                pos[1] = pw[0]; // P
                types[1] = (int) TOUCHEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = pw[1]; // Q
              } else {  // s4 < 0
                // [P, Q] contains [A, l] (tritri-23-+0-).
                types[0] = (int) SHAREVERT;
                pos[0] = pu[0]; // A
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              }
            }
          } else { // s3 < 0
            if (s2 > 0) {
              if (s4 > 0) {
                // [P, Q] in [A, l] (tritri-23-+-+).
                types[0] = (int) TOUCHFACE;
                pos[0] = 3;     // Int([A, B, C])
                pos[1] = pw[0]; // P
                types[1] = (int) TOUCHFACE;
                pos[2] = 3;     // Int([A, B, C])
                pos[3] = pw[1]; // Q
              } else {
                if (s4 == 0) {
                  // Q = l, [P, Q] in [A, l] (tritri-23-+-0).
                  types[0] = (int) TOUCHFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = pw[0]; // P
                  types[1] = (int) TOUCHEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = pw[1]; // Q
                } else { // s4 < 0
                  // [P, Q] overlaps [A, l] (tritri-23-+--).
                  types[0] = (int) TOUCHFACE;
                  pos[0] = 3;     // Int([A, B, C])
                  pos[1] = pw[0]; // P
                  types[1] = (int) ACROSSEDGE;
                  pos[2] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
                  pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
                }
              }
            } else { // s2 == 0
              // P = l in [B, C] (tritri-23#0##).
              types[0] = (int) TOUCHEDGE;
              pos[0] = (iu == 0 ? pu[1] : mi1mo3[pu[1]]); // Int([B, C])
              pos[1] = pw[0]; // P
              types[1] = (int) DISJOINT; // No 2nd intersection
            }
          }
        }
      } else {
        // Q = A (tritri-230###).
        types[0] = (int) SHAREVERT;
        pos[0] = pu[0]; // A
        pos[1] = pw[1]; // Q
        types[1] = (int) DISJOINT; // No 2nd intersection
      }
    }

    return 1;
  } // if (z1 == 2)

  if (z1 == 3) {

    if (z2 == 0) {  // (30)
      if (s1 < 0) {
        if (s3 > 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // [i, j] overlaps [A, B] (tritri-30-+++).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = 3;     // Int([P, Q, R])
            types[1] = (int) ACROSSEDGE;
            pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
            pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
          } else {
            if (s4 == 0) {
              // j = B, [i, j] contains [A, B] (tritri-30-++0).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSVERT;
              pos[2] = pu[1]; // B
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else { // s4 < 0
              // [i, j] contains [A, B] (tritri-30-++-).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSVERT;
              pos[2] = pu[1]; // B
              pos[3] = 3;     // Int([P, Q, R])
            }
          }
        } else {
          if (s3 == 0) {
            assert(s2 > 0); // SELF_CHECK
            if (s4 > 0) {
              // i = A, [i, j] in [A, B] (tritri-30-+0+).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else {
              if (s4 == 0) {
                // [i, j] = [A, B] (tritri-30-+00).
                types[0] = (int) ACROSSVERT;
                pos[0] = pu[0]; // A
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSVERT;
                pos[2] = pu[1]; // B
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else { // s4 < 0
                // i = A, [i, j] contains [A, B] (tritri-30-+0-).
                types[0] = (int) ACROSSVERT;
                pos[0] = pu[0]; // A
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSVERT;
                pos[2] = pu[1]; // B
                pos[3] = 3;     // Int([P, Q, R])
              }
            }
          } else { // s3 < 0
            if (s2 > 0) {
              if (s4 > 0) {
                // [i, j] in [A, B] (tritri-30-+-+).
                types[0] = (int) ACROSSEDGE;
                pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else {
                if (s4 == 0) {
                  // j = B, [i, j] in [A, B] (tritri-30-+-0).
                  types[0] = (int) ACROSSEDGE;
                  pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                  pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                  types[1] = (int) ACROSSVERT;
                  pos[2] = pu[1]; // B
                  pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
                } else { // s4 < 0
                  // [i, j] overlaps [A, B] (tritri-30-+--).
                  types[0] = (int) ACROSSEDGE;
                  pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                  pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
                  types[1] = (int) ACROSSVERT;
                  pos[2] = pu[1]; // B
                  pos[3] = 3;     // Int([P, Q, R])
                }
              }
            } else { // s2 == 0
              // i = B in [P, R] (tritri-30#0##).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[1]; // B
              pos[1] = (iv == 0 ? pw[2] : mi1mo3[pw[2]]); // Int([R, P])
              types[1] = (int) DISJOINT; // No 2nd intersection
            }
          }
        }
      } else {
        // j = A in [Q, R] (tritri-300###).
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[0]; // A
        pos[1] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
        types[1] = (int) DISJOINT; // No 2nd intersection
      }
    } else if (z2 == 2) {  // (32)
      if (s1 < 0) {
        if (s3 > 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // [P, j] overlaps [A, B] (tritri-32-+++).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = 3;     // Int([P, Q, R])
            types[1] = (int) ACROSSEDGE;
            pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
            pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
          } else {
            if (s4 == 0) {
              // j = B, [P, j] contains [A, B] (tritri-32-++0).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSVERT;
              pos[2] = pu[1]; // B
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else { // s4 < 0
              // [P, j] contains [A, B] (tritri-32-++-).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = 3;     // Int([P, Q, R])
              types[1] = (int) ACROSSVERT;
              pos[2] = pu[1]; // B
              pos[3] = 3;     // Int([P, Q, R])
            }
          }
        } else {
          if (s3 == 0) {
            assert(s2 > 0); // SELF_CHECK
            if (s4 > 0) {
              // P = A, [P, j] in [A, B] (tritri-32-+0+).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pw[0]; // P
              types[1] = (int) ACROSSEDGE;
              pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
              pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
            } else {
              if (s4 == 0) {
                // [P, j] = [A, B] (tritri-32-+00).
                types[0] = (int) SHAREVERT;
                pos[0] = pu[0]; // A
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSVERT;
                pos[2] = pu[1]; // B
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else { // s4 < 0
                // P = A, [P, j] contains [A, B] (tritri-32-+0-).
                types[0] = (int) SHAREVERT;
                pos[0] = pu[0]; // A
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSVERT;
                pos[2] = pu[1]; // B
                pos[3] = 3;     // Int([P, Q, R])
              }
            }
          } else { // s3 < 0
            if (s2 > 0) {
              if (s4 > 0) {
                // [P, j] in [A, B] (tritri-32-+-+).
                types[0] = (int) TOUCHEDGE;
                pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSEDGE;
                pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
              } else {
                if (s4 == 0) {
                  // j = B, [P, j] in [A, B] (tritri-32-+-0).
                  types[0] = (int) TOUCHEDGE;
                  pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                  pos[1] = pw[0]; // P
                  types[1] = (int) ACROSSVERT;
                  pos[2] = pu[1]; // B
                  pos[3] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
                } else { // s4 < 0
                  // [P, j] overlaps [A, B] (tritri-32-+--).
                  types[0] = (int) TOUCHEDGE;
                  pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                  pos[1] = pw[0]; // P
                  types[1] = (int) ACROSSVERT;
                  pos[2] = pu[1]; // B
                  pos[3] = 3;     // Int([P, Q, R])
                }
              }
            } else { // s2 == 0
              // P = B (tritri-32#0##).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[1]; // B
              pos[1] = pw[0]; // P
              types[1] = (int) DISJOINT; // No 2nd intersection
            }
          }
        }
      } else {
        // j = A in [Q, R] (tritri-320###).
        types[0] = (int) ACROSSVERT;
        pos[0] = pu[0]; // A
        pos[1] = (iv == 0 ? pw[1] : mi1mo3[pw[1]]); // Int([Q, R])
        types[1] = (int) DISJOINT; // No 2nd intersection
      }
    } else if (z2 == 3) {  // (33)
      if (s1 < 0) {
        if (s3 > 0) {
          assert(s2 > 0); // SELF_CHECK
          if (s4 > 0) {
            // [P, Q] overlaps [A, B] (tritri-33-+++).
            types[0] = (int) ACROSSVERT;
            pos[0] = pu[0]; // A
            pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
            types[1] = (int) TOUCHEDGE;
            pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
            pos[3] = pw[1]; // Q
          } else {
            if (s4 == 0) {
              // Q = B, [P, Q] contains [A, B] (tritri-33-++0).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              types[1] = (int) SHAREVERT;
              pos[2] = pu[1]; // B
              pos[3] = pw[1]; // Q
            } else { // s4 < 0
              // [P, Q] contains [A, B] (tritri-33-++-).
              types[0] = (int) ACROSSVERT;
              pos[0] = pu[0]; // A
              pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              types[1] = (int) ACROSSVERT;
              pos[2] = pu[1]; // B
              pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
            }
          }
        } else {
          if (s3 == 0) {
            assert(s2 > 0); // SELF_CHECK
            if (s4 > 0) {
              // P = A, [P, Q] in [A, B] (tritri-33-+0+).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[0]; // A
              pos[1] = pw[0]; // P
              types[1] = (int) TOUCHEDGE;
              pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
              pos[3] = pw[1]; // Q
            } else {
              if (s4 == 0) {
                // [P, Q] = [A, B] (tritri-33-+00).
                types[0] = (int) SHAREEDGE;
                pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
                types[1] = (int) SHAREEDGE;
                pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                pos[1] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              } else { // s4 < 0
                // P = A, [P, Q] contains [A, B] (tritri-33-+0-).
                types[0] = (int) SHAREVERT;
                pos[0] = pu[0]; // A
                pos[1] = pw[0]; // P
                types[1] = (int) ACROSSVERT;
                pos[2] = pu[0]; // A
                pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
              }
            }
          } else { // s3 < 0
            if (s2 > 0) {
              if (s4 > 0) {
                // [P, Q] in [A, B] (tritri-33-+-+).
                types[0] = (int) TOUCHEDGE;
                pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                pos[1] = pw[0]; // P
                types[1] = (int) TOUCHEDGE;
                pos[2] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                pos[3] = pw[1]; // Q
              } else {
                if (s4 == 0) {
                  // Q = B, [P, Q] overlaps [A, B] (tritri-33-+-0).
                  types[0] = (int) TOUCHEDGE;
                  pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                  pos[1] = pw[0]; // P
                  types[1] = (int) SHAREVERT;
                  pos[2] = pu[1]; // B
                  pos[3] = pw[1]; // Q
                } else { // s4 < 0
                  // [P, Q] overlaps [A, B] (tritri-33-+--).
                  types[0] = (int) TOUCHEDGE;
                  pos[0] = (iu == 0 ? pu[0] : mi1mo3[pu[0]]); // Int([A, B])
                  pos[1] = pw[0]; // P
                  types[1] = (int) ACROSSVERT;
                  pos[2] = pu[1]; // B
                  pos[3] = (iv == 0 ? pw[0] : mi1mo3[pw[0]]); // Int([P, Q])
                }
              }
            } else { // s2 == 0
              // P = B (tritri-33#0##).
              types[0] = (int) SHAREVERT;
              pos[0] = pu[1]; // B
              pos[1] = pw[0]; // P
              types[1] = (int) DISJOINT;  // No 2nd intersection
            }
          }
        }
      } else {
        // Q = A (tritri-330###).
        types[0] = (int) SHAREVERT;
        pos[0] = pu[0]; // A
        pos[1] = pw[1]; // Q
        types[1] = (int) DISJOINT;  // No 2nd intersection
      }
    }

    return 1;
  } // if (z1 == 3)
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insphere_sos()    Insphere test with symbolic perturbation.               //
//                                                                           //
// Given four points pa, pb, pc, and pd, test if the point pe lies inside or //
// outside the circumscirbed sphere of the four points.  Here we assume that //
// the orientation of the sequence {pa, pb, pc, pd} is negative (NOT zero),  //
// i.e., pd lies at the negative side of the plane defined by pa, pb, and pc.//
//                                                                           //
// Return a positive value (> 0) if pe lies outside, a negative value (< 0)  //
// if pe lies inside the sphere, the returned value will not be zero.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::insphere_sos(point pa, point pb, point pc, point pd, point pe)
{
  REAL sign;

  inspherecount++;

  sign = insphere(pa, pb, pc, pd, pe);
  if (sign != 0.0) {
    return sign;
  }

  insphere_sos_count++;

  // Symbolic perturbation.
  point pt[5], swappt;
  REAL oriA, oriB;
  int swaps, count;
  int n, i;

  pt[0] = pa;
  pt[1] = pb;
  pt[2] = pc;
  pt[3] = pd;
  pt[4] = pe;
  
  // Sort the five points such that their indices are in the increasing
  //   order. An optimized bubble sort algorithm is used, i.e., it has
  //   the worst case O(n^2) runtime, but it is usually much faster.
  swaps = 0; // Record the total number of swaps.
  n = 5;
  do {
    count = 0;
    n = n - 1;
    for (i = 0; i < n; i++) {
      if (pointmark(pt[i]) > pointmark(pt[i+1])) {
        swappt = pt[i]; pt[i] = pt[i+1]; pt[i+1] = swappt;
        count++;
      }
    }
    swaps += count;
  } while (count > 0); // Continue if some points are swapped.

  oriA = orient3d(pt[1], pt[2], pt[3], pt[4]);
  if (oriA != 0.0) {
    // Flip the sign if there are odd number of swaps.
    if ((swaps % 2) != 0) oriA = -oriA;
    return oriA;
  }
  
  oriB = -orient3d(pt[0], pt[2], pt[3], pt[4]);
  assert(oriB != 0.0); // SELF_CHECK
  // Flip the sign if there are odd number of swaps.
  if ((swaps % 2) != 0) oriB = -oriB;
  return oriB;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// isincircle()    Test if d lies inside the circumcircle of abc.            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::incircle3d(point pa, point pb, point pc, point pd)
{
  point pk, pl, pm, pn;
  REAL area2[4], n[3], c[3];
  REAL amax, sign, r, l, q;
  int imax;

  // Calculate the areas of the four triangles in a, b, c, and d.
  //   Get the base triangle which has the largest area.
  facenormal(pa, pb, pc, n, 1);
  area2[0] = DOT(n, n); 
  facenormal(pb, pa, pd, n, 1);
  area2[1] = DOT(n, n);
  if (area2[0] < area2[1]) {
    amax = area2[1]; imax = 1;
  } else {
    amax = area2[0]; imax = 0;
  }
  facenormal(pc, pd, pb, n, 1);
  area2[2] = DOT(n, n);
  if (amax < area2[2]) {
    amax = area2[2]; imax = 2;
  }
  facenormal(pd, pc, pa, n, 1);
  area2[3] = DOT(n, n);
  if (amax < area2[3]) {
    amax = area2[3]; imax = 3;
  }

  // Permute the vertices s. t. the base triangle is (pk, pl, pm).
  if (imax == 0) {
    pk = pa; pl = pb; pm = pc; pn = pd; sign = 1.0;
  } else if (imax == 1) {
    pk = pb; pl = pa; pm = pd; pn = pc; sign = 1.0;
  } else if (imax == 2) {
    pk = pc; pl = pd; pm = pb; pn = pa; sign = -1.0;
  } else {
    pk = pd; pl = pc; pm = pa; pn = pb; sign = -1.0;
  }

  // Make sure that the base triangle is not degenerate.
  l = DIST(pk, pl);
  l += DIST(pl, pm);
  l += DIST(pm, pk);
  l /= 3.0;

  if (sqrt(amax) > (l * l * b->epsilon)) {
    // Calculate the circumcenter and radius.
    circumsphere(pk, pl, pm, NULL, c, &r);
    l = DIST(c, pn);
    q = fabs((l - r) / r);
  } else {
    // A (nearly) degenerate base triangle.
    assert(0);  // Not handle yet.
    q = 0;
  }

  if (q > b->epsilon) {
    return (l - r) * sign;  // Adjust the sign.
  } else {
    return 0;  // Round to zero.
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// iscoplanar()    Check if four points are approximately coplanar.          //
//                                                                           //
// 'tol' is the relative error tolerance.  The coplanarity is determined by  //
// the equation q < tol, where q = fabs(6 * vol) / L^3, vol is the volume of //
// the tet klmn, and L is the average edge length of the tet.                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::iscoplanar(point k, point l, point m, point n, REAL ori)
{
  REAL L, q;

  L = DIST(k, l);
  L += DIST(l, m);
  L += DIST(m, k);
  L += DIST(k, n);
  L += DIST(l, n);
  L += DIST(m, n);
  assert(L > 0.0);  // SELF_CHECK
  L /= 6.0;

  q = fabs(ori) / (L * L * L);

  return q <= b->epsilon;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// interiorangle()    Return the interior angle (0 - 2 * PI) between vectors //
//                    o->p1 and o->p2.                                       //
//                                                                           //
// 'n' is the normal of the plane containing face (o, p1, p2).  The interior //
// angle is the total angle rotating from o->p1 around n to o->p2.  Exchange //
// the position of p1 and p2 will get the complement angle of the other one. //
// i.e., interiorangle(o, p1, p2) = 2 * PI - interiorangle(o, p2, p1).  Set  //
// 'n' be NULL if you only want the interior angle between 0 - PI.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::interiorangle(point o, point p1, point p2, REAL* n)
{
  REAL v1[3], v2[3], np[3];
  REAL theta, costheta, lenlen;
  REAL ori, len1, len2;

  // Get the interior angle (0 - PI) between o->p1, and o->p2.
  v1[0] = p1[0] - o[0];
  v1[1] = p1[1] - o[1];
  v1[2] = p1[2] - o[2];
  v2[0] = p2[0] - o[0];
  v2[1] = p2[1] - o[1];
  v2[2] = p2[2] - o[2];
  len1 = sqrt(DOT(v1, v1));
  len2 = sqrt(DOT(v2, v2));
  lenlen = len1 * len2;
  assert(lenlen != 0.0);  // SELF_CHECK
  costheta = DOT(v1, v2) / lenlen;
  if (costheta > 1.0) {
    costheta = 1.0; // Roundoff. 
  } else if (costheta < -1.0) {
    costheta = -1.0; // Roundoff. 
  }
  theta = acos(costheta);
  if (n != NULL) {
    // Get a point above the face (o, p1, p2);
    np[0] = o[0] + n[0];
    np[1] = o[1] + n[1];
    np[2] = o[2] + n[2];
    // Adjust theta (0 - 2 * PI).
    ori = orient3d(p1, o, np, p2);
    if (ori > 0.0) {
      theta = 2 * PI - theta;
    }
  }

  return theta;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// facenormal()    Calculate the normal of the face.                         //
//                                                                           //
// The normal of the face abc can be calculated by the cross product of 2 of //
// its 3 edge vectors.  A better choice of two edge vectors will reduce the  //
// numerical error during the calculation.  Burdakov proved that the optimal //
// basis problem is equivalent to the minimum spanning tree problem with the //
// edge length be the functional, see Burdakov, "A greedy algorithm for the  //
// optimal basis problem", BIT 37:3 (1997), 591-599. If 'pivot' > 0, the two //
// short edges in abc are chosen for the calculation.                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::facenormal(point pa, point pb, point pc, REAL *n, int pivot)
{
  REAL v1[3], v2[3], v3[3], *pv1, *pv2;
  REAL L1, L2, L3;

  v1[0] = pb[0] - pa[0];  // edge vector v1: a->b
  v1[1] = pb[1] - pa[1];
  v1[2] = pb[2] - pa[2];
  v2[0] = pa[0] - pc[0];  // edge vector v2: c->a
  v2[1] = pa[1] - pc[1];
  v2[2] = pa[2] - pc[2];

  // Default, normal is calculated by: v1 x (-v2) (see Fig. fnormal).
  if (pivot > 0) {
    // Choose edge vectors by Burdakov's algorithm.
    v3[0] = pc[0] - pb[0];  // edge vector v3: b->c
    v3[1] = pc[1] - pb[1];
    v3[2] = pc[2] - pb[2];
    L1 = DOT(v1, v1);
    L2 = DOT(v2, v2);
    L3 = DOT(v3, v3);
    // Sort the three edge lengths.
    if (L1 < L2) {
      if (L2 < L3) {
        pv1 = v1; pv2 = v2; // n = v1 x (-v2).
      } else {
        pv1 = v3; pv2 = v1; // n = v3 x (-v1).
      }
    } else {
      if (L1 < L3) {
        pv1 = v1; pv2 = v2; // n = v1 x (-v2).
      } else {
        pv1 = v2; pv2 = v3; // n = v2 x (-v3).
      }
    }
  } else {
    pv1 = v1; pv2 = v2; // n = v1 x (-v2).
  }

  // Calculate the face normal.
  CROSS(pv1, pv2, n);
  // Inverse the direction;
  n[0] = -n[0];
  n[1] = -n[1];
  n[2] = -n[2];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// circumsphere()    Calculate the smallest circumsphere (center and radius) //
//                   of the given three or four points.                      //
//                                                                           //
// The circumsphere of four points (a tetrahedron) is unique if they are not //
// degenerate. If 'pd == NULL', the smallest circumsphere of three points is //
// the diametral sphere of the triangle (pa, pb, pc).                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::circumsphere(point pa, point pb, point pc, point pd,
  REAL* cent, REAL* radius)
{
  REAL A[4][4], rhs[4], D;
  int indx[4];

  // Compute the coefficient matrix A (3x3).
  A[0][0] = pb[0] - pa[0];
  A[0][1] = pb[1] - pa[1];
  A[0][2] = pb[2] - pa[2];
  A[1][0] = pc[0] - pa[0];
  A[1][1] = pc[1] - pa[1];
  A[1][2] = pc[2] - pa[2];
  if (pd != NULL) {
    A[2][0] = pd[0] - pa[0];
    A[2][1] = pd[1] - pa[1]; 
    A[2][2] = pd[2] - pa[2];
  } else {
    CROSS(A[0], A[1], A[2]);
  }

  // Compute the right hand side vector b (3x1).
  rhs[0] = 0.5 * DOT(A[0], A[0]);
  rhs[1] = 0.5 * DOT(A[1], A[1]);
  if (pd != NULL) {
    rhs[2] = 0.5 * DOT(A[2], A[2]);
  } else {
    rhs[2] = 0.0;
  }

  // Solve the 3 by 3 equations use LU decomposition with partial pivoting
  //   and backward and forward substitute..
  if (!lu_decmp(A, 3, indx, &D, 0)) {
    assert(0);  // No solution.
  }
  lu_solve(A, 3, indx, rhs, 0);
  if (cent != (REAL *) NULL) {
    cent[0] = pa[0] + rhs[0];
    cent[1] = pa[1] + rhs[1];
    cent[2] = pa[2] + rhs[2];
  }
  if (radius != (REAL *) NULL) {
    *radius = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lu_decmp()    Compute the LU decomposition of a matrix.                   //
//                                                                           //
// Compute the LU decomposition of a (non-singular) square matrix A using    //
// partial pivoting and implicit row exchanges.  The result is:              //
//     A = P * L * U,                                                        //
// where P is a permutation matrix, L is unit lower triangular, and U is     //
// upper triangular.  The factored form of A is used in combination with     //
// 'lu_solve()' to solve linear equations: Ax = b, or invert a matrix.       //
//                                                                           //
// The inputs are a square matrix 'lu[N..n+N-1][N..n+N-1]', it's size is 'n'.//
// On output, 'lu' is replaced by the LU decomposition of a rowwise permuta- //
// tion of itself, 'ps[N..n+N-1]' is an output vector that records the row   //
// permutation effected by the partial pivoting, effectively,  'ps' array    //
// tells the user what the permutation matrix P is; 'd' is output as +1/-1   //
// depending on whether the number of row interchanges was even or odd,      //
// respectively.                                                             //
//                                                                           //
// Return true if the LU decomposition is successfully computed, otherwise,  //
// return false in case that A is a singular matrix.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::lu_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N)
{
  REAL scales[4];
  REAL pivot, biggest, mult, tempf;
  int pivotindex = 0;
  int i, j, k;

  *d = 1.0;                                      // No row interchanges yet.

  for (i = N; i < n + N; i++) {                             // For each row.
    // Find the largest element in each row for row equilibration
    biggest = 0.0;
    for (j = N; j < n + N; j++)
      if (biggest < (tempf = fabs(lu[i][j])))
        biggest  = tempf;
    if (biggest != 0.0)
      scales[i] = 1.0 / biggest;
    else {
      scales[i] = 0.0;
      return false;                            // Zero row: singular matrix.
    }
    ps[i] = i;                                 // Initialize pivot sequence.
  }

  for (k = N; k < n + N - 1; k++) {                      // For each column.
    // Find the largest element in each column to pivot around.
    biggest = 0.0;
    for (i = k; i < n + N; i++) {
      if (biggest < (tempf = fabs(lu[ps[i]][k]) * scales[ps[i]])) {
        biggest = tempf;
        pivotindex = i;
      }
    }
    if (biggest == 0.0) {
      return false;                         // Zero column: singular matrix.
    }
    if (pivotindex != k) {                         // Update pivot sequence.
      j = ps[k];
      ps[k] = ps[pivotindex];
      ps[pivotindex] = j;
      *d = -(*d);                          // ...and change the parity of d.
    }

    // Pivot, eliminating an extra variable  each time
    pivot = lu[ps[k]][k];
    for (i = k + 1; i < n + N; i++) {
      lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;
      if (mult != 0.0) {
        for (j = k + 1; j < n + N; j++)
          lu[ps[i]][j] -= mult * lu[ps[k]][j];
      }
    }
  }

  // (lu[ps[n + N - 1]][n + N - 1] == 0.0) ==> A is singular.
  return lu[ps[n + N - 1]][n + N - 1] != 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lu_solve()    Solves the linear equation:  Ax = b,  after the matrix A    //
//               has been decomposed into the lower and upper triangular     //
//               matrices L and U, where A = LU.                             //
//                                                                           //
// 'lu[N..n+N-1][N..n+N-1]' is input, not as the matrix 'A' but rather as    //
// its LU decomposition, computed by the routine 'lu_decmp'; 'ps[N..n+N-1]'  //
// is input as the permutation vector returned by 'lu_decmp';  'b[N..n+N-1]' //
// is input as the right-hand side vector, and returns with the solution     //
// vector. 'lu', 'n', and 'ps' are not modified by this routine and can be   //
// left in place for successive calls with different right-hand sides 'b'.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::lu_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N)
{
  int i, j;
  REAL X[4], dot;

  for (i = N; i < n + N; i++) X[i] = 0.0;

  // Vector reduction using U triangular matrix.
  for (i = N; i < n + N; i++) {
    dot = 0.0;
    for (j = N; j < i + N; j++)
      dot += lu[ps[i]][j] * X[j];
    X[i] = b[ps[i]] - dot;
  }

  // Back substitution, in L triangular matrix.
  for (i = n + N - 1; i >= N; i--) {
    dot = 0.0;
    for (j = i + 1; j < n + N; j++)
      dot += lu[ps[i]][j] * X[j];
    X[i] = (X[i] - dot) / lu[ps[i]][i];
  }

  for (i = N; i < n + N; i++) b[i] = X[i];
}

#endif // #ifndef geomCXX