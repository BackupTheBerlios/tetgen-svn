#ifndef geomCXX
#define geomCXX

#include "tetgen.h"

REAL tetgenmesh::PI = 3.14159265358979323846264338327950288419716939937510582;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_vert_inter()    Test whether a triangle (ABC) and a coplanar vertex   //
//                     (P) intersect or not.                                 //
//                                                                           //
// ABC is not degenerate.  We know that R lies above ABC, i.e., ABCR is a    //
// valid tetrahedron.                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::tri_vert_inter(point A, point B,
  point C, point P, point R, int* pos)
{
  REAL s3, s4, s5;

  // Test P's orientations wrt edges AB, BC, CA. 
  s3 = orient3d(A, B, R, P);
  s4 = orient3d(B, C, R, P);
  s5 = orient3d(C, A, R, P);
  orient3dcount+=3;

  if (b->epsilon) {
    // Re-evaluate the sign with respect to the tolerance.
    if ((s3 != 0) && iscoplanar(A, B, R, P, s3)) s3 = 0;
    if ((s4 != 0) && iscoplanar(B, C, R, P, s4)) s4 = 0;
    if ((s5 != 0) && iscoplanar(C, A, R, P, s5)) s5 = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-vert cop (%d %d %d)-(%d)-%d, s3=%g, s4=%g, s5=%g.\n",
      pointmark(A), pointmark(B), pointmark(C), pointmark(P), pointmark(R),
      s3, s4, s5);
  }
  trivercopcount++;

  if (s3 < 0) {
    return DISJOINT;
  } else {
    if (s4 < 0) {
      return DISJOINT;
    } else {
      if (s5 < 0) {
        return DISJOINT;
      }
    }
  }

  // P lies eother inside or on the boundary of ABC.
  if (s3 == 0) {
    if (s4 == 0) {
      if (s5 == 0) {
        assert(0); // Not possible.
      }
      // P = B.
      *pos = 1;
      return ACROSSVERT;
    }
    if (s5 == 0) {
      // P = A.
      *pos = 0;
      return ACROSSVERT;
    }
    // P lies on AB.
    *pos = 0;
    return ACROSSEDGE;
  }
  if (s4 == 0) {
    if (s5 == 0) {
      // P = C.
      *pos = 2;
      return ACROSSVERT;
    }
    // P lies on BC.
    *pos = 1;
    return ACROSSEDGE;
  }
  if (s5 == 0) {
    // P lies on CA.
    *pos = 1;
    return ACROSSEDGE;
  }

  // P lies inside ABC.
  return ACROSSFACE;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_edge_inter_cop()    Test whether a triangle (ABC) and a coplanar edge //
//                         (PQ) intersect or not.                            //
//                                                                           //
// ABC and PQ are not degenerate.  We know that R lies above ABC, i.e., ABCR //
// is a valid tetrahedron.                                                   //
//                                                                           //
// In addition, we assume that both P and Q locate outside the triangle ABC. //
// This is the case which we search intersecting faces along a missing seg.  //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::tri_edge_inter_cop(point A, point B,
  point C, point P, point Q, point R, int *pos)
{
  REAL s3, s4, s5;
  REAL s3s4, s4s5, s5s3;
  REAL p_ab, p_bc, p_ca;
  REAL q_ab, q_bc, q_ca;

  // Test A's, B's, and C's orientations wrt plane PQR. 
  s3 = orient3d(P, Q, R, A);
  s4 = orient3d(P, Q, R, B);
  s5 = orient3d(P, Q, R, C);
  orient3dcount+=3;

  if (b->epsilon) {
    // Re-evaluate the sign with respect to the tolerance.
    if ((s3 != 0) && iscoplanar(P, Q, R, A, s3)) s3 = 0;
    if ((s4 != 0) && iscoplanar(P, Q, R, B, s4)) s4 = 0;
    if ((s5 != 0) && iscoplanar(P, Q, R, C, s5)) s5 = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-edge cop (%d %d %d)-(%d, %d)-%d, s3=%g, s4=%g, s5=%g.\n",
      pointmark(A), pointmark(B), pointmark(C), pointmark(P), pointmark(Q),
      pointmark(R), s3, s4, s5);
  }
  triedgcopcount++;

  s3s4 = s3 * s4;
  s4s5 = s4 * s5;
  s5s3 = s5 * s3;

  // Check DISJOINT cases (see Fig. 2).
  if (s3s4 > 0) {
    if (s4s5 > 0) {
      if (!(s5s3 > 0)) {
        assert(0); // Impossible case.
      }
      return DISJOINT;  // (+++) (---).
    }
  }

  // Test the orientations of P wrt ABR, BCR, and CAR.
  p_ab = orient3d(A, B, R, P);
  p_bc = orient3d(B, C, R, P);
  p_ca = orient3d(C, A, R, P);
  orient3dcount+=3;
  // Test the orientations of Q wrt ABR, BCR, and CAR.
  q_ab = orient3d(A, B, R, Q);
  q_bc = orient3d(B, C, R, Q);
  q_ca = orient3d(C, A, R, Q);
  orient3dcount+=3;

  if (b->epsilon) {
    if ((p_ab != 0) && iscoplanar(A, B, R, P, p_ab)) p_ab = 0;
    if ((p_bc != 0) && iscoplanar(B, C, R, P, p_bc)) p_bc = 0;
    if ((p_ca != 0) && iscoplanar(C, A, R, P, p_ca)) p_ca = 0;
    if ((q_ab != 0) && iscoplanar(A, B, R, Q, q_ab)) q_ab = 0;
    if ((q_bc != 0) && iscoplanar(B, C, R, Q, q_bc)) q_bc = 0;
    if ((q_ca != 0) && iscoplanar(C, A, R, Q, q_ca)) q_ca = 0;
  }

  if (b->verbose > 2) {
    printf("        p_ab=%g, p_bc=%g, p_ca=%g.\n", p_ab, p_bc, p_ca);
    printf("        q_ab=%g, q_bc=%g, q_ca=%g.\n", q_ab, q_bc, q_ca);
  }

  // Check the DISJOINT cases (see Fig. 1).
  if (p_ab < 0) {
    if (q_ab < 0) {
      // PQ lies above ABR.
      return DISJOINT;  // (--####)
    }
  }
  if (p_bc < 0) {
    if (q_bc < 0) {
      // PQ lies above BCR.
      return DISJOINT; // (##--##)
    }
  }
  if (p_ca < 0) {
    if (q_ca < 0) {
      // PQ lies above CAR.
      return DISJOINT; // (####--)
    }
  }

  if (p_ab < 0) {
    if (q_ab > 0) {
      if (s3 == 0) {
        if (s4 == 0) {
          // PQ is collinear with AB. Impossible for p_ab < 0.
          assert(0);
        }
        if (s5 == 0) {
          // PQ is collinear with CA.
        }
        // PQ intersects A (see Fig. 4 top).
        *pos = 0;
        return ACROSSVERT;
      }
      if (s4 == 0) {
        if (s5 == 0) {
          // PQ is collinear with BC.
        }
        // PQ intersects B (see Fig. 4 top).
        *pos = 1;
        return ACROSSVERT;
      }
      // (-+####) PQ intersect AB (see Fig. 4 top).
      *pos = 0;
      return ACROSSEDGE;
    } else {
      // q_ab = 0. This case has been ruled out by our assumption.
      assert(0);
    }
  } else {
    if (p_ab == 0) {
      // This case has been ruled out by our assumption.
      assert(0);
    }
  }

  if (p_bc < 0) {
    if (q_bc > 0) {
      if (s4 == 0) {
        if (s5 == 0) {
          // PQ is collinear with BC. Impossible by p_bc < 0.
          assert(0);
        }
        if (s3 == 0) {
          // PQ is collinear with AB.
        }
        // PQ intersects B.
        *pos = 1;
        return ACROSSVERT;
      }
      if (s5 == 0) {
        if (s3 == 0) {
          // PQ is collinear with CA.
        }
        // PQ intersects C.
        *pos = 2;
        return ACROSSVERT;
      }
      // PQ intersects BC.
      *pos = 1;
      return ACROSSEDGE;
    } else { 
      // q_bc == 0; Not possible by our assumption.
      assert(0);
    }
  } else {
    if (p_bc == 0) {
      // Not possible by our assumption.
      assert(0);
    }
  }

  assert(p_ca < 0);
  assert(q_ca > 0);

  if (s5 == 0) {
    if (s3 == 0) {
      // PQ is collinear with CA. Impossible by p_ca < 0.
      assert(0);
    }
    if (s4 == 0) {
      // PQ is collinear with BC.
    }
    // PQ intersects C.
    *pos = 2;
    return ACROSSVERT;
  }
  if (s3 == 0) {
    if (s4 == 0) {
      // PQ is collinear with AB.
    }
    // PQ intersects A.
    *pos = 0;
    return ACROSSVERT;
  }
  // PQ intersects CA.
  *pos = 2;
  return ACROSSEDGE;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_edge_inter_tail()    Test whether a triangle (ABC) and a non-coplanar //
//                          edge (PQ) intersect or not.                      //
//                                                                           //
// ABC and PQ are not degenerate.  We know that P lies above ABC, i.e., ABCP //
// is a valid tetrahedron, and Q lies below ABC.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::tri_edge_inter_tail(point A, point B,
  point C, point P, point Q, int *pos)
{
  REAL s3, s4, s5;

  s3 = orient3d(A, B, P, Q);
  s4 = orient3d(B, C, P, Q);
  s5 = orient3d(C, A, P, Q);
  orient3dcount+=3;

  if (b->epsilon) {
    // Re-evaluate the sign with respect to the tolerance.
    if ((s3 != 0) && iscoplanar(A, B, P, Q, s3)) s3 = 0;
    if ((s4 != 0) && iscoplanar(B, C, P, Q, s4)) s4 = 0;
    if ((s5 != 0) && iscoplanar(C, A, P, Q, s5)) s5 = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-edge tail (%d %d %d)-(%d, %d), s3=%g, s4=%g, s5=%g.\n",
      pointmark(A), pointmark(B), pointmark(C), pointmark(P), pointmark(Q),
      s3, s4, s5);
  }

  // Use the signs to decide whether PQ intersects ABC (see Fig. tri_edge).
  //   There are total 3^3 = 27 cases, some are impossible (by assumptions).
  if (s3 > 0) {
    if (s4 > 0) {
      if (s5 > 0) {
        // (+++) PQ intersects ABC in its interior.
        return ACROSSFACE;
      } else {
        if (s5 < 0) {
          // (++-) PQ lies outside of CAP.
          return DISJOINT;
        } else {
          // (++0) PQ intersects wih edge CA.
          *pos = 2;
          return ACROSSEDGE;
        }
      }
    } else {
      if (s4 < 0) {
        if (s5 > 0) {
          // (+-+) PQ lies outside of BCP.
          return DISJOINT;
        } else {
          if (s5 < 0) {
            // (+--) PQ lies both outsides of BCP and CAP.
            return DISJOINT;
          } else {
            // (+-0) PQ lies both outsides of BCP, and coplanar with CAP.
            return DISJOINT;
          }
        }
      } else {
        if (s5 > 0) {
          // (+0+) PQ intersects edge BC.
          *pos = 1;
          return ACROSSEDGE;
        } else {
          if (s5 < 0) {
            // (+0-) PQ lies both outsides of BCP and CAP.
            return DISJOINT;
          } else {
            // (+00) PQ is collinear with PC.
            *pos = 2;
            return ACROSSVERT;
          }
        }
      }
    }
  } else {
    if (s3 < 0) {
      if (s4 > 0) {
        if (s5 > 0) {
          // (-++) PQ lies outside ABP.
          return DISJOINT;
        } else {
          if (s5 < 0) {
            // (-+-) PQ lies both outsides of ABP and CAP.
            return DISJOINT;
          } else {
            // (-+0) PQ lies outside of ABP, and coplanar with CAP.
            return DISJOINT;
          }
        }
      } else {
        if (s4 < 0) {
          if (s5 > 0) {
            // (--+) PQ lies both outsides of ABP and BCP.
            return DISJOINT;
          } else {
            if (s5 < 0) {
              // (---) PQ lies both outsides of ABP, BCP, and CAP.
              assert(0);  // Impossible by assumption.
            } else {
              // (--0) PQ lies both outsides of ABP, BCP, and cop. with CAP.
              assert(0);  // Impossible by assumption.
            }
          }
        } else {
          if (s5 > 0) {
            // (-0+) PQ lies outsides of ABP, and coplanar with BCP.
            return DISJOINT;
          } else {
            if (s5 < 0) {
              // (-0-) PQ lies both outsides of ABP, CAP, and cop with BCP.
              assert(0);  // Impossible by assumption.
            } else {
              // (-00) PQ lies collinear with CP, and Q lies above ABC.
              assert(0);  // Impossible by assumption.
            }
          }
        }
      }
    } else {
      if (s4 > 0) {
        if (s5 > 0) {
          // (0++) PQ intersects edge AB.
          *pos = 0;
          return ACROSSEDGE;
        } else {
          if (s5 < 0) {
            // (0+-) PQ lies outside ACP.
            return DISJOINT;
          } else {
            // (0+0) PQ is coincident with edge PA.
            *pos = 0;
            return ACROSSVERT;
          }
        }
      } else {
        if (s4 < 0) {
          if (s5 > 0) {
            // (0-+) PQ lies outside of BCP.
            return DISJOINT;
          } else {
            if (s5 < 0) {
              // (0--) Q must lie above ABC.
              assert(0);  // Impossible by assumption.
            } else {
              // (0-0) Q must lie above ABC. Q is collinaer with AP.
              assert(0);  // Impossible by assumption.
            }
          }
        } else {
          if (s5 > 0) {
            // (00+) PQ is collinear with edge PB.
            *pos = 1;
            return ACROSSVERT;
          } else {
            if (s5 < 0) {
              // (00-) Q must lies above ABC.
              assert(0);  // Impossible by assumption.
            } else {
              // (000) Q = P.
              assert(0);  // Impossible by assumption.
            }
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tri_edge_inter()    Test whether a triangle (ABC) and an edge (PQ)        //
//                     intersect or not.                                     //
//                                                                           //
// R is a reference point which lies (strictly) above ABC.                   //
//                                                                           //
// Returned DISJOINT if they don't intersect. Otherwise, they intersect, and //
// 'pos' indicates the intersection position.                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::tri_edge_inter(point A, point B,
  point C, point P, point Q, point R, int *pos)
{
  REAL s1, s2, s1s2;

  // Test the locations of P and Q with respect to ABC.
  s1 = orient3d(A, B, C, P);
  s2 = orient3d(A, B, C, Q);
  orient3dcount+=2;

  if (b->epsilon > 0) {
    if ((s1 != 0) && iscoplanar(A, B, C, P, s1)) s1 = 0;
    if ((s2 != 0) && iscoplanar(A, B, C, Q, s2)) s2 = 0;
  }

  if (b->verbose > 2) {
    printf("      Tri-edge test (%d %d %d)-(%d, %d), s1=%g, s2=%g.\n",
      pointmark(A), pointmark(B), pointmark(C), pointmark(P), pointmark(Q),
      s1, s2);
  }
  triedgcount++;

  s1s2 = s1 * s2;

  // Classify the 3^2 = 9 cases.
  if (s1s2 > 0) {
    // Both P and Q lie at the same side of ABC.
    return DISJOINT;  // (++) or (--)
  } else {
    if (s1s2 < 0) {
      if (s1 < 0) {
        return tri_edge_inter_tail(A, B, C, P, Q, pos);  // (-+)
      } else {
        return tri_edge_inter_tail(A, B, C, Q, P, pos);  // (+-)
      }
    } else {
      if (s1 == 0) {
        if (s2 == 0) {
          return tri_edge_inter_cop(A, B, C, P, Q, R, pos); // (00)
        } else {
          return tri_vert_inter(A, B, C, P, R, pos); // (0+) or (0-) 
        }
      }
      if (s2 == 0) {
        return tri_vert_inter(A, B, C, Q, R, pos); // (+0) or (-0)
      }
    }
  }
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