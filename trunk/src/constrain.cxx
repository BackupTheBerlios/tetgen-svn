#ifndef constrainCXX
#define constrainCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// finddirection()    Find the tet on the path from one point to another.    //
//                                                                           //
// The path starts from 'searchtet''s origin and ends at 'endpt'. On finish, //
// 'searchtet' contains a tet on the path, its origin does not change.       //
//                                                                           //
// The return value indicates one of the following cases (let 'searchtet' be //
// abcd, a is the origin of the path):                                       //
//   - ACROSSVERT, edge ab is collinear with the path;                       //
//   - ACROSSEDGE, edge bc intersects with the path;                         //
//   - ACROSSFACE, face bcd intersects with the path.                        //
//                                                                           //
// WARNING: This routine is designed for convex triangulations, and will not //
// generally work after the holes and concavities have been carved.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::finddirection(triface* searchtet, 
  point endpt)
{
  triface neightet;
  point pa, pb, pc, pd, pn;
  enum {HMOVE, RMOVE, LMOVE} nextmove;
  enum {HCOPLANE, RCOPLANE, LCOPLANE, NCOPLANE} cop;
  REAL hori, rori, lori;
  REAL dmin, dist;

  tetrahedron ptr;
  int *iptr, tver;

  // The origin is fixed.
  pa = org(*searchtet);
  if ((point) searchtet->tet[7] == dummypoint) {
    // A hull tet. Choose the neighbor of its base face.
    searchtet->loc = 0;
    symself(*searchtet);
    // Reset the origin to be pa.
    if ((point) searchtet->tet[4] == pa) {
      searchtet->loc = 0; searchtet->ver = 0;
    } else if ((point) searchtet->tet[5] == pa) {
      searchtet->loc = 0; searchtet->ver = 2;
    } else if ((point) searchtet->tet[6] == pa) {
      searchtet->loc = 0; searchtet->ver = 4;
    } else {
      assert((point) searchtet->tet[7] == pa); // SELF_CHECK
      searchtet->loc = 1; searchtet->ver = 2;
    }
  }
  if (searchtet->ver & 01) {
    // Switch to the 0th edge ring.
    esymself(*searchtet);
    enextself(*searchtet);
  }
  pb = dest(*searchtet);
  pc = apex(*searchtet);

  // Check whether the destination or apex is 'endpt'.
  if (pb == endpt) {
    // pa->pb is the search edge.
    return ACROSSVERT;
  }
  if (pc == endpt) {
    // pa->pc is the search edge.
    enext2self(*searchtet);
    esymself(*searchtet);
    return ACROSSVERT;
  }

  // Walk through tets at pa until the right one is found.
  while (1) {

    pd = oppo(*searchtet);

    if (b->verbose > 2) {
      printf("      From tet (%d, %d, %d, %d) to %d.\n", pointmark(pa),
        pointmark(pb), pointmark(pc), pointmark(pd), pointmark(endpt));
    }

    // Check whether the opposite vertex is 'endpt'.
    if (pd == endpt) {
      // pa->pd is the search edge.
      enext0fnextself(*searchtet);
      enext2self(*searchtet);
      esymself(*searchtet);
      return ACROSSVERT;
    }
    // Check if we have entered outside of the domain.
    if (pd == dummypoint) {
      assert(0);
    }

    // Now assume that the base face abc coincides with the horizon plane,
    //   and d lies above the horizon.  The search point 'endpt' may lie
    //   above or below the horizon.  We test the orientations of 'endpt'
    //   with respect to three planes: abc (horizon), bad (right plane),
    //   and acd (left plane). 
    hori = orient3d(pa, pb, pc, endpt);
    rori = orient3d(pb, pa, pd, endpt);
    lori = orient3d(pa, pc, pd, endpt);
    orient3dcount += 3;

    // Now decide the tet to move.  It is possible there are more than one
    //   tet are viable moves. Use the opposite points of thier neighbors
    //   to discriminate, i.e., we choose the tet whose opposite point has
    //   the shortest distance to 'endpt'.
    if (hori > 0) {
      if (rori > 0) {
        if (lori > 0) {
          // Any of the three neighbors is a viable move.
          nextmove = HMOVE;
          sym(*searchtet, neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dmin = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dmin = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
          }
          enext0fnext(*searchtet, neightet);
          symself(neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dist = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dist = dmin;
          }
          if (dist < dmin) {
            nextmove = RMOVE;
            dmin = dist;
          }
          enext2fnext(*searchtet, neightet);
          symself(neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dist = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dist = dmin;
          }
          if (dist < dmin) {
            nextmove = LMOVE;
            dmin = dist;
          }
        } else {
          // Two tets, below horizon and below right, are viable.
          nextmove = HMOVE;
          sym(*searchtet, neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dmin = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dmin = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
          }
          enext0fnext(*searchtet, neightet);
          symself(neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dist = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dist = dmin;
          }
          if (dist < dmin) {
            nextmove = RMOVE;
            dmin = dist;
          }
        }
      } else {
        if (lori > 0) {
          // Two tets, below horizon and below left, are viable.
          nextmove = HMOVE;
          sym(*searchtet, neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dmin = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dmin = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
          }
          enext2fnext(*searchtet, neightet);
          symself(neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dist = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dist = dmin;
          }
          if (dist < dmin) {
            nextmove = LMOVE;
            dmin = dist;
          }
        } else {
          // The tet below horizon is chosen.
          nextmove = HMOVE;
        }
      }
    } else {
      if (rori > 0) {
        if (lori > 0) {
          // Two tets, below right and below left, are viable.
          nextmove = RMOVE;
          enext0fnext(*searchtet, neightet);
          symself(neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dmin = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dmin = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
          }
          enext2fnext(*searchtet, neightet);
          symself(neightet);
          pn = oppo(neightet);
          if (pn != dummypoint) {
            dist = NORM2(endpt[0] - pn[0], endpt[1] - pn[1], endpt[2] - pn[2]);
          } else {
            dist = dmin;
          }
          if (dist < dmin) {
            nextmove = LMOVE;
            dmin = dist;
          }
        } else {
          // The tet below right is chosen.
          nextmove = RMOVE;
        }
      } else {
        if (lori > 0) {
          // The tet below left is chosen.
          nextmove = LMOVE;
        } else {
          // 'endpt' lies either on the plane(s) or across face bcd.
          if (hori == 0) {
            if (rori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pb.
              return ACROSSVERT;
            }
            if (lori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pc.
              enext2self(*searchtet);
              esymself(*searchtet);
              return ACROSSVERT;
            }
            // pa->'endpt' crosses the edge pb->pc.
            // enextself(*searchtet);
            // return ACROSSEDGE;
            cop = HCOPLANE;
            break;
          }
          if (rori == 0) {
            if (lori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pd.
              enext0fnextself(*searchtet); // face abd.
              enext2self(*searchtet);
              esymself(*searchtet);
              return ACROSSVERT;
            }
            // pa->'endpt' crosses the edge pb->pd.
            // enext0fnextself(*searchtet); // face abd.
            // enextself(*searchtet);
            // return ACROSSEDGE;
            cop = RCOPLANE;
            break;
          }
          if (lori == 0) {
            // pa->'endpt' crosses the edge pc->pd.
            // enext2fnextself(*searchtet);  // face cad
            // enext2self(*searchtet);
            // return ACROSSEDGE;
            cop = LCOPLANE;
            break;
          }
          // pa->'endpt' crosses the face bcd.
          // enextfnextself(*searchtet);
          // return ACROSSFACE;
          cop = NCOPLANE;
          break;
        }
      }
    }

    // Move to the next tet, fix pa as its origin.
    if (nextmove == RMOVE) {
      fnextself(*searchtet);
    } else if (nextmove == LMOVE) {
      enext2self(*searchtet);
      fnextself(*searchtet);
      enextself(*searchtet);
    } else { // HMOVE
      symedgeself(*searchtet);
      enextself(*searchtet);
    }
    assert(org(*searchtet) == pa); // SELF_CHECK
    pb = dest(*searchtet);
    pc = apex(*searchtet);

  } // while (1)

  // Either case ACROSSEDGE or ACROSSFACE.
  if (b->epsilon > 0) {
    // Use tolerance to re-evaluate the orientations.
    if (cop != HCOPLANE) {
      if (iscoplanar(pa, pb, pc, endpt, hori)) hori = 0;
    }
    if (cop != RCOPLANE) {
      if (iscoplanar(pb, pa, pd, endpt, rori)) rori = 0;
    }
    if (cop != LCOPLANE) {
      if (iscoplanar(pa, pc, pd, endpt, lori)) lori = 0;
    }
    // It is not possible that all orientations are zero.
    assert(!((hori == 0) && (rori == 0) && (lori == 0))); // SELF_CHECK
  }
  
  // Now decide the degenerate cases.
  if (hori == 0) {
    if (rori == 0) {
      // pa->'endpt' is COLLINEAR with pa->pb.
      return ACROSSVERT;
    }
    if (lori == 0) {
      // pa->'endpt' is COLLINEAR with pa->pc.
      enext2self(*searchtet);
      esymself(*searchtet);
      return ACROSSVERT;
    }
    // pa->'endpt' crosses the edge pb->pc.
    return ACROSSEDGE;
  }
  if (rori == 0) {
    if (lori == 0) {
      // pa->'endpt' is COLLINEAR with pa->pd.
      enext0fnextself(*searchtet); // face abd.
      enext2self(*searchtet);
      esymself(*searchtet);
      return ACROSSVERT;
    }
    // pa->'endpt' crosses the edge pb->pd.
    enext0fnextself(*searchtet); // face abd.
    esymself(*searchtet);
    enextself(*searchtet);
    return ACROSSEDGE;
  }
  if (lori == 0) {
    // pa->'endpt' crosses the edge pc->pd.
    enext2fnextself(*searchtet);  // face cad
    esymself(*searchtet);
    return ACROSSEDGE;
  }
  // pa->'endpt' crosses the face bcd.
  return ACROSSFACE;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutsegment()    Look for a given segment in the tetrahedralization T.   //
//                                                                           //
// Search an edge in the tetrahedralization that matches the given segmment. //
// If such an edge is found, the segment is 'locked' at the edge.            //
//                                                                           //
// If 'searchtet' != NULL, it's origin must be the origin of 'sseg'.  It is  //
// used as the starting tet for searching the edge.                          //
//                                                                           //
// The returned value indicates one of the following cases:                  //
//   - SHAREVERT, the segment exists and is inserted in T;                   //
//   - ACROSSVERT, a vertex ('refpt') lies on the segment;                   //
//   - ACROSSEDGE, the segment is missing;                                   //
//   - ACROSSFACE, the segment is missing;                                   //
//                                                                           //
// If the returned value is ACROSSEDGE or ACROSSFACE, i.e., the segment is   //
// missing, 'refpt' returns the reference point for splitting thus segment,  //
// 'searchtet' returns a tet containing the 'refpt'.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::scoutsegment(face* sseg,
  triface* searchtet, point* refpt)
{
  triface neightet, reftet;
  face splitsh, checkseg;
  point startpt, endpt;
  point pa, pb, pc, pd;
  enum location loc;
  enum intersection dir;
  REAL angmax, ang;
  long facecount;
  bool orgflag;
  int types[2], poss[4];
  int shver, pos, i;

  tetrahedron ptr;
  shellface sptr;
  int *iptr, tver;

  // Is 'searchtet' a valid handle?
  if (searchtet->tet == NULL) {
    orgflag = false;
    // Search a tet whose origin is one of the endpoints of 'sseg'.
    for (shver = 0; shver < 2 && !orgflag; shver++) {
      startpt = (point) sseg->sh[shver + 3];
      decode(point2tet(startpt), *searchtet);
      if ((searchtet->tet != NULL) && (searchtet->tet[4] != NULL)) {
        // Check if this tet contains pa.
        for (i = 4; i < 8 && !orgflag; i++) {
          if ((point) searchtet->tet[i] == startpt) {
            // Found. Set pa as its origin.
            switch (i) {
              case 4: searchtet->loc = 0; searchtet->ver = 0; break;
              case 5: searchtet->loc = 0; searchtet->ver = 2; break;
              case 6: searchtet->loc = 0; searchtet->ver = 4; break;
              case 7: searchtet->loc = 1; searchtet->ver = 2; break;
            }
            sseg->shver = shver;
            orgflag = true;
          }
        }
      }
    }
    if (!orgflag) {
      // Locate startpt in tetrahedralization.
      startpt = sorg(*sseg);
      randomsample(startpt, searchtet);
      loc = locate(startpt, searchtet);
      assert(loc == ONVERTEX);  // SELF_CHECK
      force_ptloc_count++;
    }
  } else {
    startpt = sorg(*sseg);
    assert(org(*searchtet) == startpt); // SELF_CHECK
  }
  endpt = sdest(*sseg);

  if (b->verbose > 1) {
    printf("    Scout seg (%d, %d).\n", pointmark(startpt), pointmark(endpt));
  }

  dir = finddirection(searchtet, endpt);

  if (dir == ACROSSVERT) {
    pd = dest(*searchtet);
    if (pd == endpt) {
      // Found! Insert the segment.
      tsspivot(*searchtet, checkseg);  // SELF_CHECK
      assert(checkseg.sh == NULL);  // SELF_CHECK
      neightet = *searchtet;
      do {
        tssbond1(neightet, *sseg);
        fnextself(neightet);
      } while (neightet.tet != searchtet->tet);
      // The job is done. 
      return SHAREVERT;
    } else {
      // A point is on the path.
      *refpt = pd;
      return ACROSSVERT;
    }
  }

  if (b->verbose > 1) {
    printf("    Scout ref point of seg (%d, %d).\n", pointmark(startpt), 
      pointmark(endpt));
  }
  facecount = across_face_count;

  enextfnextself(*searchtet); // Go to the opposite face.
  symedgeself(*searchtet); // Enter the adjacent tet.

  pa = org(*searchtet);
  angmax = interiorangle(pa, startpt, endpt, NULL);
  *refpt = pa;
  pb = dest(*searchtet);
  ang = interiorangle(pb, startpt, endpt, NULL);
  if (ang > angmax) {
    angmax = ang;
    *refpt = pb;
  }

  // Check whether two segments are intersecting.
  if (dir == ACROSSEDGE) {
    tsspivot(*searchtet, checkseg);
    if (checkseg.sh != NULL) {
      printf("  Invalid PLC!  Two segments intersect each other.\n");
      startpt = farsorg(*sseg);
      endpt = farsdest(*sseg);
      pa = farsorg(checkseg);
      pb = farsdest(checkseg);
      printf("    1st: (%d, %d), 2nd: (%d, %d).\n", pointmark(startpt), 
        pointmark(endpt), pointmark(pa), pointmark(pb));
      terminatetetgen(1);
    }
    across_edge_count++;
  }
  
  pc = apex(*searchtet);
  ang = interiorangle(pc, startpt, endpt, NULL);
  if (ang > angmax) {
    angmax = ang;
    *refpt = pc;
  }
  reftet = *searchtet; // Save the tet containing the refpt.

  // Search intersecting faces along the segment.
  while (1) {

    pd = oppo(*searchtet);
    assert(pd != dummypoint);  // SELF_CHECK

    if (b->verbose > 2) {
      printf("      Passing face (%d, %d, %d, %d), dir(%d).\n", pointmark(pa),
        pointmark(pb), pointmark(pc), pointmark(pd), (int) dir);
    }
    across_face_count++;

    // Stop if we meet 'endpt'.
    if (pd == endpt) break;

    ang = interiorangle(pd, startpt, endpt, NULL);
    if (ang > angmax) {
      angmax = ang;
      *refpt = pd;
      reftet = *searchtet;
    }

    // Find a face intersecting the segment.
    if (dir == ACROSSFACE) {
      // One of the three oppo faces in 'searchtet' intersects the segment.
      neightet.tet = searchtet->tet;
      neightet.ver = 0;
      for (i = 0; i < 3; i++) {
        neightet.loc = locpivot[searchtet->loc][i];
        pa = org(neightet);
        pb = dest(neightet);
        pc = apex(neightet);
        pd = oppo(neightet); // The above point.
        if (tri_edge_test(pa, pb, pc, startpt, endpt, pd, 1, types, poss)) {
          dir = (enum intersection) types[0];
          pos = poss[0];
          break;
        } else {
          dir = DISJOINT;
          pos = 0;
        }
      }
      assert(dir != DISJOINT);  // SELF_CHECK
    } else { // dir == ACROSSEDGE
      // Check the two opposite faces (of the edge) in 'searchtet'.
      neightet = *searchtet;
      neightet.ver = 0;
      for (i = 0; i < 2; i++) {
        neightet.loc = locverpivot[searchtet->loc][searchtet->ver][i];
        pa = org(neightet);
        pb = dest(neightet);
        pc = apex(neightet);
        pd = oppo(neightet); // The above point.
        if (tri_edge_test(pa, pb, pc, startpt, endpt, pd, 1, types, poss)) {
          dir = (enum intersection) types[0];
          pos = poss[0];
          break;
        } else {
          dir = DISJOINT;
          pos = 0;
        }
      }
      if (dir == DISJOINT) {
        // No intersection. Go to the next tet.
        dir = ACROSSEDGE;
        fnextself(*searchtet);
        continue;
      }
    }

    if (dir == ACROSSVERT) {
      // This segment passing a vertex. Choose it and return.
      for (i = 0; i < pos; i++) {
        enextself(neightet);
      }
      pd = org(neightet);
      if (b->verbose > 2) {
        angmax = interiorangle(pd, startpt, endpt, NULL);
      }
      *refpt = pd;
      break;
    }
    if (dir == ACROSSEDGE) {
      // Get the edge intersects with the segment.
      for (i = 0; i < pos; i++) {
        enextself(neightet);
      }
    }
    // Go to the next tet.
    symedge(neightet, *searchtet);

    if (dir == ACROSSEDGE) {
      // Check whether two segments are intersecting.
      tsspivot(*searchtet, checkseg);
      if (checkseg.sh != NULL) {
        printf("  Invalid PLC!  Two segments intersect each other.\n");
        startpt = farsorg(*sseg);
        endpt = farsdest(*sseg);
        pa = farsorg(checkseg);
        pb = farsdest(checkseg);
        printf("    1st: (%d, %d), 2nd: (%d, %d).\n", pointmark(startpt), 
          pointmark(endpt), pointmark(pa), pointmark(pb));
        terminatetetgen(1);
      }
      across_edge_count++;
    }

  } // while (1)

  // dir is either ACROSSVERT, or ACROSSEDGE, or ACROSSFACE.
  if (b->verbose > 2) {
    printf("      Refpt %d (%g), visited %ld faces.\n", pointmark(*refpt),
      angmax / PI * 180.0, (int) dir, across_face_count - facecount);
  }
  if (across_face_count - facecount > across_max_count) {
    across_max_count = across_face_count - facecount;
  }

  *searchtet = reftet;
  return dir;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// getsegmentsplitpoint()    Calculate a split point in the given segment.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::getsegmentsplitpoint(face* sseg, point refpt, REAL* vt)
{
  point ei, ej, ek;
  REAL split, L, d, d1, d2, d3;
  int stype, sign;
  int i;

  // Decide the type of this segment.
  sign = 1;
  ei = sorg(*sseg);
  ej = sdest(*sseg);

  if (getpointtype(ei) == ACUTEVERTEX) {
    if (getpointtype(ej) == ACUTEVERTEX) {
      // Both ei and ej are ACUTEVERTEX.
      stype = 0;
    } else {
      // ej is either a RIDGEVERTEX or a STEINERVERTEX.
      stype = 1;
    }
  } else {
    if (getpointtype(ei) == RIDGEVERTEX) {
      if (getpointtype(ej) == ACUTEVERTEX) {
        stype = 1; sign = -1;
      } else {
        if (getpointtype(ej) == RIDGEVERTEX) {
          // Both ei and ej are non-acute.
          stype = 0;
        } else {
          // ej is a STEINERVETEX.
          ek = farsdest(*sseg);
          if (getpointtype(ek) == ACUTEVERTEX) {
            stype = 1; sign = -1;
          } else {
            stype = 0;
          }
        }
      }
    } else {
      // ei is a STEINERVERTEX.
      if (getpointtype(ej) == ACUTEVERTEX) {
        stype = 1; sign = -1;
      } else {
        ek = farsorg(*sseg);
        if (getpointtype(ej) == RIDGEVERTEX) {
          if (getpointtype(ek) == ACUTEVERTEX) {
            stype = 1;
          } else {
            stype = 0;
          }
        } else {
          // Both ei and ej are STEINERVETEXs. ei has priority.
          if (getpointtype(ek) == ACUTEVERTEX) {
            stype = 1;
          } else {
            ek = farsdest(*sseg);
            if (getpointtype(ek) == ACUTEVERTEX) {
              stype = 1; sign = -1;
            } else {
              stype = 0;
            }
          }
        }
      }
    }
  }

  // Adjust the endpoints: ei, ej.
  if (sign == -1) {
    sesymself(*sseg);
    ei = sorg(*sseg);
    ej = sdest(*sseg);
  }

  if (b->verbose > 1) {
    printf("    Split a type-%d seg(%d, %d) ref(%d)", stype,
      pointmark(ei), pointmark(ej), pointmark(refpt));
    if (stype) {
      ek = farsorg(*sseg);
      printf(" ek(%d)", pointmark(ek));
    }
    printf(".\n");
  }

  // Calculate the split point.
  if (stype == 0) {
    // Use rule-1.
    L = DIST(ei, ej);
    d1 = DIST(ei, refpt);
    d2 = DIST(ej, refpt);
    if (d1 < d2) {
      // Choose ei as center.
      if (d1 < 0.5 * L) {
        split = d1 / L;
      } else {
        split = 0.5;
      }
      for (i = 0; i < 3; i++) {
        vt[i] = ei[i] + split * (ej[i] - ei[i]);
      }
    } else {
      // Choose ej as center.
      if (d2 < 0.5 * L) {
        split = d2 / L;
      } else {
        split = 0.5;
      }
      for (i = 0; i < 3; i++) {
        vt[i] = ej[i] + split * (ei[i] - ej[i]);
      }
    }
    r1count++;
  } else {
    // Use rule-2.
    ek = farsorg(*sseg);
    L = DIST(ek, ej);
    d = DIST(ek, refpt);
    split = d / L;
    for (i = 0; i < 3; i++) {
      vt[i] = ek[i] + split * (ej[i] - ek[i]);
    }
    d1 = DIST(vt, refpt);
    d2 = DIST(vt, ej);
    if (d1 > d2) {
      // Use rule-3.
      d3 = DIST(ei, refpt);
      if (d1 < 0.5 * d3) {
        split = (d - d1) / L;
      } else {
        split = (d - 0.5 * d3) / L;
      }
      for (i = 0; i < 3; i++) {
        vt[i] = ek[i] + split * (ej[i] - ek[i]);
      }
    }
    d1 > d2 ? r3count++ : r2count++;
  }

  if (b->verbose > 1) {
    printf("    split (%g), vt (%g, %g, %g).\n", split, vt[0], vt[1], vt[2]);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// markacutevertices()    Classify vertices as ACUTEVERTEXs or RIDGEVERTEXs. //
//                                                                           //
// Initially, all segment vertices are marked as VOLVERTEX (after calling    //
// incrementaldelaunay()).                                                   //
//                                                                           //
// A segment vertex is ACUTEVERTEX if it two segments incident it form an    //
// interior angle less than 60 degree, otherwise, it is a RIDGEVERTEX.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::markacutevertices()
{
  point pa, pb, pc;
  REAL anglimit, ang;
  bool acuteflag;
  int acutecount;
  int idx, i, j;

  face* segperverlist;
  int* idx2seglist;

  if (b->verbose) {
    printf("  Marking acute vertices.\n");
  }

  // Construct a map from points to segments.
  makepoint2submap(subsegpool, idx2seglist, segperverlist);

  anglimit = PI / 3.0;  // 60 degree.
  acutecount = 0;

  // Loop over the set of vertices.
  pointpool->traversalinit();
  pa = pointtraverse();
  while (pa != NULL) {
    idx = pointmark(pa) - in->firstnumber;
    // Mark it if it is an endpoint of some segments.
    if (idx2seglist[idx + 1] > idx2seglist[idx]) {
      acuteflag = false;
      // Do a brute-force pair-pair check.
      for (i = idx2seglist[idx]; i < idx2seglist[idx + 1] && !acuteflag; i++) {
        pb = sdest(segperverlist[i]);
        for (j = i + 1; j < idx2seglist[idx + 1] && !acuteflag; j++) {
          pc = sdest(segperverlist[j]);
          ang = interiorangle(pa, pb, pc, NULL);
          acuteflag = ang < anglimit;
        }
      }
      // Now mark the vertex.
      if (b->verbose > 1) {
        printf("    Mark %d as %s.\n", pointmark(pa), acuteflag ?
          "ACUTEVERTEX" : "RIDGEVERTEX");
      }
      setpointtype(pa, acuteflag ? ACUTEVERTEX : RIDGEVERTEX);
      acutecount += (acuteflag ? 1 : 0);
    }
    pa = pointtraverse();
  }

  if (b->verbose) {
    printf("  %d acute vertices.\n", acutecount);
  }

  delete [] idx2seglist;
  delete [] segperverlist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// delaunizesegments()    Recover segments in a Delaunay tetrahedralization. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::delaunizesegments()
{
  triface searchtet;
  face *psseg, sseg, nsseg, splitshs[2];
  point refpt, newpt;
  enum intersection dir;
  int s;

  shellface sptr;

  if (b->verbose) {
    printf("  Delaunizing segments.\n");
  }

  // Loop until 'subsegstack' is empty.
  while (subsegstack->objects > 0l) {
    // seglist is used as a stack.
    subsegstack->objects--;
    psseg = (face *) fastlookup(subsegstack, subsegstack->objects);
    sseg = *psseg;

    if (!sinfected(sseg)) continue; // Not a missing segment.
    suninfect(sseg);

    // Insert the segment.
    searchtet.tet = NULL;
    dir = scoutsegment(&sseg, &searchtet, &refpt);

    if (dir != SHAREVERT) {
      // The segment is missing, split it.
      spivot(sseg, splitshs[0]);
      if (dir != ACROSSVERT) {
        // Create the new point.
        makepoint(&newpt);
        getsegmentsplitpoint(&sseg, refpt, newpt);
        // Split the segment by newpt.
        flipn2nf(newpt, splitshs, 1);
        setpointtype(newpt, STEINERVERTEX);
        sspivot(splitshs[0], sseg);
        sspivot(splitshs[1], nsseg);
        // Some subfaces may be non-Delaunay.
        lawsonflip();
        // Insert newpt into the DT.
        insertvertex(newpt, &searchtet, true);
      } else {
        // Split the segment by refpt.
        flipn2nf(refpt, splitshs, 1);
        if (getpointtype(refpt) != ACUTEVERTEX) {
          setpointtype(refpt, RIDGEVERTEX);
        }
        sspivot(splitshs[0], sseg);
        sspivot(splitshs[1], nsseg);
        lawsonflip();
      }
      // Add two subsegments into the stack.
      if (b->order == 4) {  // '-o4' option (for debug)
        sinfect(sseg);
        subsegstack->newindex((void **) &psseg);
        *psseg = sseg;
        sinfect(nsseg);
        subsegstack->newindex((void **) &psseg);
        *psseg = nsseg;
      } else {
        s = randomnation(subsegstack->objects);
        subsegstack->newindex((void **) &psseg);
        *psseg = * (face *) fastlookup(subsegstack, s);
        sinfect(sseg); 
        psseg = (face *) fastlookup(subsegstack, s);
        *psseg = sseg;
        s = randomnation(subsegstack->objects);
        subsegstack->newindex((void **) &psseg);
        *psseg = * (face *) fastlookup(subsegstack, s);
        sinfect(nsseg);
        psseg = (face *) fastlookup(subsegstack, s);
        *psseg = nsseg;
      }
    }
  }

  if (b->verbose) {
    printf("  %d protecting points.\n", r1count + r2count + r3count);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutsubface()    Look for a given subface in the tetrahedralization T.   //
//                                                                           //
// 'ssub' is the subface, denoted as abc. If abc exists in T, it is 'locked' //
// at the place where the two tets sharing at it.                            //
//                                                                           //
// The returned value indicates one of the following cases:                  //
//   - SHAREFACE, abc exists and is inserted;                                //
//   - TOUCHEDGE, a vertex (the origin of 'searchtet') lies on ab.           //
//   - EDGETRIINT, all three edges of abc are missing.                       //
//   - ACROSSTET, a tet (in 'searchtet') crosses the facet containg abc.     //
//                                                                           //
// If the retunred value is ACROSSTET, the subface is missing.  'searchtet'  //
// returns a tet which shares the same edge as 'pssub'.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::scoutsubface(face* pssub,
  triface* searchtet)
{
  triface spintet;
  face checksh;
  point pa, pb, pc, pd;
  enum intersection dir;
  int i;

  tetrahedron ptr;
  int tver;

  if (searchtet->tet == NULL) {
    // Search an edge of 'ssub' in tetrahedralization.
    pssub->shver = 0;
    for (i = 0; i < 3; i++) {
      pa = sorg(*pssub);
      pb = sdest(*pssub);
      // Do not search dummypoint.
      assert(pa != dummypoint); // SELF_CHECK
      assert(pb != dummypoint); // SELF_CHECK
      // Get a tet whose origin is pa.
      decode(point2tet(pa), *searchtet);
      assert(searchtet->tet != NULL); // SELF_CHECK
      if ((point) searchtet->tet[4] == pa) {
        searchtet->loc = 0; searchtet->ver = 0;
      } else if ((point) searchtet->tet[5] == pa) {
        searchtet->loc = 0; searchtet->ver = 2;
      } else if ((point) searchtet->tet[6] == pa) {
        searchtet->loc = 0; searchtet->ver = 4;
      } else {
        assert((point) searchtet->tet[7] == pa); // SELF_CHECK
        searchtet->loc = 1; searchtet->ver = 2;
      }
      // Search the edge from pa->pb.
      dir = finddirection(searchtet, pb);
      if (dir == ACROSSVERT) {
        if (dest(*searchtet) == pb) {
          // Found the edge. Break the loop.
          break;
        } else {
          // A vertex lies on the search edge. Return it.
          enextself(*searchtet);
          return TOUCHEDGE;
        }
      }
      senextself(*pssub);
    }
    if (i == 3) {
      // None of the three edges exists.
      return EDGETRIINT; // ab intersects the face in 'searchtet'.
    }
  } else {
    // 'searchtet' holds the current edge of 'pssub'.
    pa = org(*searchtet);
    pb = dest(*searchtet);
  }

  pc = sapex(*pssub);

  if (b->verbose > 1) {
    printf("    Scout subface (%d, %d, %d) (%ld).\n", pointmark(pa),
      pointmark(pb), pointmark(pc), subfacstack->objects);
  }

  // Searchtet holds edge pa->pb. Search a face with apex pc.
  spintet = *searchtet;
  while (1) {
    fnextself(spintet);
    pd = apex(spintet);  // pd may be dummypoint. Search the face anyway.
    if (pd == pc) {
      // Found! Insert the subface.
      tspivot(spintet, checksh); // SELF_CHECK
      assert(checksh.sh == NULL); // SELF_CHECK
      tsbond(spintet, *pssub);
      symedgeself(spintet);
      tspivot(spintet, checksh); // SELF_CHECK
      assert(checksh.sh == NULL); // SELF_CHECK
      tsbond(spintet, *pssub);
      return SHAREFACE;
    }
    if (pd == apex(*searchtet)) break;
  }

  return ACROSSTET;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutcrosstet()    Scout a tetrahedron across a facet.                    //
//                                                                           //
// A subface (abc) of the facet (F) is given in 'misregion[0]', abc is not a //
// Delaunay face and is intersected by some Delaunay tets. 'searchtet' holds //
// the edge ab, it is the tet starting the search.                           //
//                                                                           //
// The subface (abc) was produced by a 2D CDT algorithm under the Assumption //
// that F is flat. In real data, however, F may not be strictly flat.  Hence //
// a tet (abde) that crosses abc may be in one of the two cases: (i) abde    //
// intersects F in its interior, or (ii) abde intersects F on its boundary.  //
// In case (i) F (or part of it) is missing in DT and needs to be recovered. //
// In (ii) F is not missing, the surface mesh of F needs to be adjusted.     //
//                                                                           //
// This routine distinguishes the two cases by the returned value, which is  //
//   - ACROSSTET, if it is case (i), 'searchtet' is abde, d and e lies below //
//     and above abc, respectively, neither d nor e is dummypoint; or        //
//   - ACROSSFACE, if it is case (ii), 'searchtet' is abde, where the face   //
//     abd intersects abc, e lies above abc, e may be dummypoint.            //
//                                                                           //
// On return, 'misregion' contains all subfaces of F, and all subfaces of F  //
// are marked as tested. All vertuces of F are infected, they will be unmark //
// or uninfect in later called functions.                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::scoutcrosstet(arraypool* misregion, 
  triface* searchtet)
{
  triface spintet;
  face *pssub, worksh, neighsh;
  face checkseg;
  point pa, pb, pc, pd, pe;
  REAL ori;
  int i, j;

  // Get the first missing subface abc.
  pssub = (face *) fastlookup(misregion, 0);
  pa = sorg(*pssub);
  pb = sdest(*pssub);
  pc = sapex(*pssub);

  // Infect the vertices of this face.
  pinfect(pa);
  pinfect(pb);
  pinfect(pc);
  // Mark this face as tested.
  smarktest(*pssub);

  // Infect all vertices of the facet.
  for (i = 0; i < misregion->objects; i++) {
    worksh = * (face *) fastlookup(misregion, i);
    for (j = 0; j < 2; j++) {
      senextself(worksh);
      sspivot(worksh, checkseg);
      if (checkseg.sh == NULL) {
        spivot(worksh, neighsh);
        assert(neighsh.sh != NULL); // SELF_CHECK
        if (!smarktested(neighsh)) {
          pd = sapex(neighsh);
          if (!pinfected(pd)) {
            pinfect(pd);
          }
          smarktest(neighsh);
          misregion->newindex((void **) &pssub);
          *pssub = neighsh;
        }
      }
    }
  }

  // Search an edge crossing the facet containing abc.
  if (searchtet->ver & 01) {
    // Adjust to 0th edge ring.
    esymself(*searchtet);
    pssub = (face *) fastlookup(misregion, 0);
    sesymself(*pssub);
    pa = org(*searchtet);
    pb = dest(*searchtet);
  }

  // Get a face containing ab and its apex lies below abc.
  pd = apex(*searchtet);
  spintet = *searchtet;
  while (1) {
    if (pd == dummypoint) break;
    ori = orient3d(pa, pb, pc, pd);
    if ((ori != 0) && pinfected(pd)) {
      ori = 0; // F is not flat. Force d be coplanar with abc.
    }
    if (ori > 0) break;
    fnextself(spintet);
    pd = apex(spintet);
    assert(pd != apex(*searchtet)); // SELF_CHECK
  }
  // Search a tet whose apex->oppo crosses the facet containig abc.
  while (1) {
    pe = oppo(spintet);
    ori = orient3d(pa, pb, pc, pe);
    if ((ori != 0) && pinfected(pe)) {
      ori = 0; // Force pe be coplanar with abc.
    }
    if (ori <= 0) break;  // stop at pd->pe.    
    fnextself(spintet);
  }

  if (ori == 0) {
    fnext(spintet, *searchtet);
    pd = apex(*searchtet);
    pe = oppo(*searchtet);
    if (b->verbose > 1) {
      printf("    Found a coplanar face (%d, %d, %d) op (%d).\n", 
        pointmark(pa), pointmark(pb), pointmark(pd), pointmark(pe));
    }
    if (getpointtype(pd) == VOLVERTEX) {
      // A vertex (pd) lies on the facet.
      enext2self(*searchtet); // org(*searchtet) == pd
      return TOUCHFACE;
    }
    return ACROSSFACE;
  } else {
    // Return a crossing tet.
    *searchtet = spintet;
    if (b->verbose > 1) {
      printf("    Found a crossing tet (%d, %d, %d, %d).\n", pointmark(pa),
        pointmark(pb), pointmark(apex(spintet)), pointmark(pe));
    }
    return ACROSSTET; // abc intersects the volume of 'searchtet'.
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// adjustsurfmesh()    Adjust surface mesh to match the volume mesh.         //
//                                                                           //
// An edge flip (flip22) is performed to adjust the surface mesh.            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::adjustsurfmesh(arraypool* misregion, triface* crosstet)
{
  face *pssub, worksh, flipfaces[2];
  face checkseg;
  point pa, pb, pc, pd, pe;
  REAL ori1, len, n[3];
  int i;

  // Get the first missing subface [a, b, c].
  pssub = (face *) fastlookup(misregion, 0);
  pa = sorg(*pssub);
  pb = sdest(*pssub);
  pc = sapex(*pssub);

  // The crosstet is abde.
  pd = apex(*crosstet);
  pe = oppo(*crosstet);

  if (pe == dummypoint) {
    // Calculate a point above the faces.
    facenormal(pa, pb, pd, n, 1);
    len = sqrt(DOT(n, n));
    n[0] /= len;
    n[1] /= len;
    n[2] /= len;
    len = DIST(pa, pb);
    len += DIST(pb, pd);
    len += DIST(pd, pa);
    len /= 3.0;
    pe[0] = pa[0] + len * n[0];
    pe[1] = pa[1] + len * n[1];
    pe[2] = pa[2] + len * n[2];
  }

  while (1) {

    ori1 = orient3d(pb, pc, pe, pd);
    assert(ori1 != 0); // SELF_CHECK

    if (ori1 < 0) { // Flip edge [b, c]
      senext(*pssub, flipfaces[0]);
    } else { // Flip edge [c, a]
      senext2(*pssub, flipfaces[0]);
    }

    sspivot(flipfaces[0], checkseg); // SELF_CHECK
    assert(checkseg.sh == NULL); // SELF_CHECK
    spivot(flipfaces[0], flipfaces[1]);
    assert(sinfected(flipfaces[1])); // SELF_CHECK

    // Temporarily uninfect them.
    suninfect(flipfaces[0]);
    suninfect(flipfaces[1]);

    flip22(flipfaces, 0);

    // Infect them back (to be recovered).
    sinfect(flipfaces[0]);
    sinfect(flipfaces[1]);

    // Find the edge [a, b].
    if (ori1 < 0) { // Flip edge [b, c]
      senext(flipfaces[1], *pssub);
    } else { // Flip edge [c, a]
      senext2(flipfaces[0], *pssub);
    }
    assert(sorg(*pssub) == pa); // SELF_CHECK
    assert(sdest(*pssub) == pb); // SELF_CHECK

    pc = sapex(*pssub);
    if (pc == pd) break;

  }

  if (pe == dummypoint) {
    pe[0] = pe[1] = pe[2] = 0;
  }

  // Uninfect (unmark) all vertices (subfaces) of the facet.
  for (i = 0; i < misregion->objects; i++) {
    worksh = * (face *) fastlookup(misregion, i);
    sunmarktest(worksh);
    pa = sorg(worksh);
    pb = sdest(worksh);
    pc = sapex(worksh);
    puninfect(pa);
    puninfect(pb);
    puninfect(pc);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// formcavity()    Form the cavity of a missing region.                      //
//                                                                           //
// On input, 'misregion' contains all subfaces of a facet F (collected in    //
// scoutcrosstet()). The first subface (misregion[0]) is the face abc.       //
// 'crosstets' contains only one crossing tet, abde, where d and e lie below //
// and above abc, respectively. Starting from this tet(abde), other crossing //
// tets are found one by one and added into 'crosstets'.                     //
//                                                                           //
// 'topfaces' and 'botfaces' return the upper and lower boundary faces of C, //
// respectively. Moreover 'topfaces[0]' is face abe, 'botfaces[0]' is bad,   //
// i.e., topfaces[0] and botfaces[0] share edge ab.                          //
//                                                                           //
// NOTE: The set of crossing tets may contain hull tets(see fig/dump-cavity- //
// case2.lua for an example).                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::formcavity(arraypool* misregion, arraypool* crosstets,
  arraypool* topfaces, arraypool* botfaces)
{
  arraypool *crossedges;
  triface *parytet, crosstet, spintet, neightet;
  face *pssub, worksh;
  face checkseg;
  point pa, pb, pc, pf, pg;
  REAL ori;
  int i, j;

  int *iptr;

  // Get the first missing subface abc.
  pssub = (face *) fastlookup(misregion, 0);
  pa = sorg(*pssub);
  pb = sdest(*pssub);
  pc = sapex(*pssub);

  // Get a crossing tet abde. 
  parytet = (triface *) fastlookup(crosstets, 0); // face abd.
  // The edge de crosses the facet. d lies below abc.
  enext2fnext(*parytet, crosstet);
  enext2self(crosstet); 
  esymself(crosstet); // the edge d->e at face [d,e,a]
  infect(crosstet);
  *parytet = crosstet; // Save it in list.

  // Temporarily re-use 'topfaces'.
  crossedges = topfaces;
  crossedges->newindex((void **) &parytet);
  *parytet = crosstet;

  // Collect all crossing tets.  Each cross tet is saved in the standard
  //   form deab, where de is a corrsing edge, orient3d(d,e,a,b) < 0. 
  // NOTE: hull tets may be collected. See fig/dump-cavity-case2a(b).lua.
  //   Make sure that neither d nor e is dummypoint.
  for (i = 0; i < crossedges->objects; i++) {
    crosstet = * (triface *) fastlookup(crossedges, i);
    // It may already be tested.
    if (!edgemarked(crosstet)) {
      // Collect all tets sharing at the edge.
      pg = apex(crosstet);
      spintet = crosstet;
      while (1) {
        // Mark this edge as tested.
        markedge(spintet);
        if (!infected(spintet)) {
          infect(spintet);
          crosstets->newindex((void **) &parytet);
          *parytet = spintet;
        }
        // Go to the neighbor tet.
        fnextself(spintet);
        // Check the validity of the PLC.
        tspivot(spintet, worksh);
        if (worksh.sh != NULL) {
          printf("Error:  Invalid PLC.\n");
          terminatetetgen(1);
        }
        if (apex(spintet) == pg) break;
      }
      // Detect new cross edges.
      while (1) {
        // Remember: spintet is edge d->e, d lies below abc.
        pf = apex(spintet);
        if (pf != dummypoint) { // Do not grab a hull edge.
          if (!pinfected(pf)) {
            // There exist a crossing edge, either d->f, or f->e.
            ori = orient3d(pa, pb, pc, pf);
            assert(ori != 0);
            if (ori < 0) {
              // The edge d->f corsses the facet.
              enext2fnext(spintet, neightet);
              esymself(neightet); // d->f.
            } else {
              // The edge f->e crosses the face.
              enextfnext(spintet, neightet);
              esymself(neightet); // f->e.
            }
            if (!edgemarked(neightet)) {
              // Add a new cross edge.
              crossedges->newindex((void **) &parytet);
              *parytet = neightet;
            }
          }
        }
        fnextself(spintet);
        if (apex(spintet) == pg) break;
      }
    }
  }

  // All cross tets are found. Unmark cross edges.
  for (i = 0; i < crossedges->objects; i++) {
    crosstet = * (triface *) fastlookup(crossedges, i);
    if (edgemarked(crosstet)) {
      // Unmark this edge.
      pg = apex(crosstet);
      spintet = crosstet;
      while (1) {
        assert(edgemarked(spintet)); // SELF_CHECK
        unmarkedge(spintet);
        // Go to the neighbor tet.
        fnextself(spintet);
        if (apex(spintet) == pg) break;
      }
    }
  }

  if (b->verbose > 1) {
    printf("    Formed cavity: %ld (%ld) cross tets (edges).\n", 
      crosstets->objects, crossedges->objects);
  }
  crossedges->restart();

  // Collect the top and bottom faces. 
  //   Remember that each cross tet was saved in the standard form: deab,
  //   where de is a corrsing edge, orient3d(d,e,a,b) < 0.
  //   topfaces[0] is abe, and botfaces[0] is bad.
  // NOTE: Hull tets may be collected. Process them as normal one.
  for (i = 0; i < crosstets->objects; i++) {
    crosstet = * (triface *) fastlookup(crosstets, i);
    enextfnext(crosstet, spintet);
    enextself(spintet);
    symedge(spintet, neightet);
    if (!infected(neightet)) {
      // A top face.
      topfaces->newindex((void **) &parytet);
      *parytet = neightet;
    }
    enext2fnext(crosstet, spintet);
    enext2self(spintet);
    symedge(spintet, neightet);
    if (!infected(neightet)) {
      // A bottom face.
      botfaces->newindex((void **) &parytet);
      *parytet = neightet;
    }
  }

  // Unmark all facet vertices.
  puninfect(pa);
  puninfect(pb);
  puninfect(pc);
  // Unmark (uninfect) all vertices (subfaces) of the facet.
  for (i = 0; i < misregion->objects; i++) {
    worksh = * (face *) fastlookup(misregion, i);
    sunmarktest(worksh);
    pf = sapex(worksh);
    puninfect(pf);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// delaunizecavity()    Fill a cavity by Delaunay tetrahedra.                //
//                                                                           //
// The tetrahedralizing cavity is the half (top or bottom part) of the whole //
// cavity.  The boundary faces of the half cavity are given in 'cavfaces',   //
// the bounday faces of the internal facet are not given.  These faces will  //
// be recovered later in fillcavity().                                       //
//                                                                           //
// This routine first constructs the DT of the vertices by the Bowyer-Watson //
// algorithm.  Then it identifies the boundary faces of the cavity in DT.    //
// The DT is returned in 'newtets'.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::delaunizecavity(arraypool *cavfaces, arraypool *newtets)
{
  triface *parytet, searchtet, neightet, spintet;
  face checksh, tmpsh;
  point pa, pb, pc, pd, pt[3];
  enum intersection dir;
  REAL ori;
  int i, j;

  tetrahedron ptr;
  int tver;

  if (b->verbose > 1) {
    printf("    Delaunizing cavity: %ld faces.\n", cavfaces->objects);
  }

  if (cavfaces->objects > maxcavsize) {
    maxcavsize = cavfaces->objects;
  }

  // Get four non-coplanar points (no dummypoint).
  parytet = (triface *) fastlookup(cavfaces, 0);
  pa = org(*parytet);
  pb = dest(*parytet);
  pc = apex(*parytet);
  pinfect(pa);
  pinfect(pb);
  pinfect(pc);
  pd = NULL;
  for (i = 1; i < cavfaces->objects; i++) {
    parytet = (triface *) fastlookup(cavfaces, i);
    pt[0] = org(*parytet);
    pt[1] = dest(*parytet);
    pt[2] = apex(*parytet);
    for (j = 0; j < 3; j++) {
      if (pt[j] != dummypoint) { // Do not include a hull point.
        if (!pinfected(pt[j])) {
          ori = orient3d(pa, pb, pc, pt[j]);
          if (ori != 0) {
            pd = pt[j];
            if (ori > 0) {  // Swap pa and pb.
              pt[j] = pa; pa = pb; pb = pt[j]; 
            }
            break;
          }
        }
      }
    }
    if (pd != NULL) break;
  }
  assert(i < cavfaces->objects); // SELF_CHECK
  pinfect(pd);

  // Create an init DT.
  initialDT(pa, pb, pc, pd);

  // Incrementally insert other vertices.
  for (i = 0; i < cavfaces->objects; i++) {
    parytet = (triface *) fastlookup(cavfaces, i);
    pt[0] = org(*parytet);
    pt[1] = dest(*parytet);
    pt[2] = apex(*parytet);
    for (j = 0; j < 3; j++) {
      if (pt[j] != dummypoint) { // Do not include a hull point.
        if (!pinfected(pt[j])) {
          searchtet = recenttet;
          insertvertex(pt[j], &searchtet, true);
          pinfect(pt[j]);
        }
      }
    }
  }

  // Collect all tets of the DT.
  marktest(recenttet);
  newtets->newindex((void **) &parytet);
  *parytet = recenttet;
  for (i = 0; i < newtets->objects; i++) {
    searchtet = * (triface *) fastlookup(newtets, i);
    for (searchtet.loc = 0; searchtet.loc < 4; searchtet.loc++) {
      sym(searchtet, neightet);
      if (!marktested(neightet)) {
        marktest(neightet);
        newtets->newindex((void **) &parytet);
        *parytet = neightet;
      }
    }
  }

  // Indentify boundary faces. Mark interior tets.
  for (i = 0; i < cavfaces->objects; i++) {
    parytet = (triface *) fastlookup(cavfaces, i);
    pt[0] = org(*parytet);
    pt[1] = dest(*parytet);
    pt[2] = apex(*parytet);
    // The face orients in its 0th edge ring.
    assert((parytet->ver & 01) == 0); // SELF_CHECK
    // Uninfect the vertices.
    for (j = 0; j < 3; j++) {
      puninfect(pt[j]);
    }
    // Does this face contain dummypoint?
    for (j = 0; j < 3; j++) {
      if (pt[j] == dummypoint) {
        // Make dummypoint be its apex.
        parytet->ver = 4;
        pt[0] = org(*parytet);
        pt[1] = dest(*parytet);
        pt[2] = apex(*parytet);
        break;
      }
    }
    // Create a temp subface.
    makeshellface(subfacepool, &tmpsh);
    setshvertices(tmpsh, pt[0], pt[1], pt[2]);
    // Remember tmpsh (use the adjacent tet slot). 
    parytet->tet[parytet->loc] = (tetrahedron) sencode(tmpsh);
    // Insert tmpsh in DT.
    searchtet.tet = NULL; 
    dir = scoutsubface(&tmpsh, &searchtet);
    if (dir != SHAREFACE) {
      assert(0); // Face unmatched. Not process yet.
    }
    // Identify the inter and outer tets at tempsh.
    stpivot(tmpsh, neightet);
    // neightet and tmpsh refer to the same edge [pt[0], pt[1]].
    //   Morover, neightet is in 0th edge ring (see decode()).
    if (org(neightet) != pt[1]) {
      symedgeself(neightet);
      assert(org(neightet) == pt[1]); // SELF_CHECK
      // Make sure that tmpsh is connected with an interior tet. 
      tsbond(neightet, tmpsh);
    }
    assert(dest(neightet) == pt[0]); // SELF_CHECK
    // Mark neightet as interior.
    if (!infected(neightet)) {
      infect(neightet);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// fillcavity()    Fill new tets into the cavity.                            //
//                                                                           //
// The new tets are stored in two disjoint sets(which share the same facet). //
// 'topfaces' and 'botfaces' are the boundaries of these two sets, respect-  //
// ively. 'midfaces' is empty on input, and will store faces in the facet.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::fillcavity(arraypool* topfaces, arraypool* botfaces,
  arraypool* midfaces, arraypool* misregion)
{
  arraypool *cavfaces;
  triface *parytet, toptet, bottet, neightet, midface;
  face worksh, checksh, tmpsh;
  face checkseg;
  point pa, pb, pc, pf, pg;
  bool bflag;
  int i, j, k;

  tetrahedron ptr;
  int *iptr, tver;

  // Connect newtets to tets outside the cavity.
  for (k = 0; k < 2; k++) {
    cavfaces = (k == 0 ? topfaces : botfaces);
    for (i = 0; i < cavfaces->objects; i++) {
      parytet = (triface *) fastlookup(cavfaces, i);
      pa = org(*parytet);
      pb = dest(*parytet);
      pc = apex(*parytet);
      // Get the saved temp subface.
      sdecode(parytet->tet[parytet->loc], tmpsh);
      // Get the adjacent new tet.
      stpivot(tmpsh, neightet);
      assert(org(neightet) == pb); // SELF_CHECK
      assert(dest(neightet) == pa); // SELF_CHECK
      // Bond the two tets.
      bond(*parytet, neightet); // Also cleared the pointer to tmpsh.
      // Bond a subface (if it exists).
      tspivot(*parytet, checksh);
      if (checksh.sh != NULL) {
        tsbond(neightet, checksh); // Also cleared the pointer to tmpsh.
      } else {
        tsdissolve(neightet); // No subface, clear the pointer to tmpsh.
      }
      // Update the point-to-tets map.
      point2tet(pa) = encode(neightet);
      point2tet(pb) = encode(neightet);
      point2tet(pc) = encode(neightet);
      // Delete the temp subface.
      shellfacedealloc(subfacepool, tmpsh.sh);
    }
  }

  // Mark all facet vertices for detecting middle subfaces.
  worksh = * (face *) fastlookup(misregion, 0);
  pa = sorg(worksh);
  pb = sdest(worksh);
  pc = sapex(worksh);
  pinfect(pa);
  pinfect(pb);
  pinfect(pc);
  for (i = 0; i < misregion->objects; i++) {
    worksh = * (face *) fastlookup(misregion, i);
    pf = sapex(worksh);
    pinfect(pf);
  }

  // The first pair of top and bottom tets share the same edge [a, b].
  toptet = * (triface *) fastlookup(topfaces, 0);
  bottet = * (triface *) fastlookup(botfaces, 0);
  symedgeself(toptet);
  symedgeself(bottet);
  // Search a subface from the top mesh.
  while (1) {
    enext0fnextself(toptet); // The next face in the same tet.
    pc = apex(toptet);
    if (pinfected(pc)) break; // [a,b,c] is a subface.
    symedgeself(toptet); // Go to the same face in the adjacent tet.
  }
  // Search the subface [a,b,c] in the bottom mesh.
  while (1) {
    enext0fnextself(bottet); // The next face in the same tet.
    pf = apex(bottet);
    if (pf == pc) break; // Face matched.
    assert(!pinfected(pf)); // SELF_CHECK
    symedgeself(bottet);
  }
  // Connect the two tets together.
  bond(toptet, bottet);
  // Both are interior tets.
  infect(toptet);
  infect(bottet);
  // Add this face into search list.
  esymself(toptet); // Choose the 0th edge ring.
  midfaces->newindex((void **) &parytet);
  *parytet = toptet;

  // Match pairs of subfaces (middle faces), connect top and bottom tets.
  for (i = 0; i < midfaces->objects; i++) {
    // Get a matched middle face [a, b, c]
    midface = * (triface *) fastlookup(midfaces, i);
    // It is inside the cavity.
    assert(marktested(midface)); // SELF_CHECK 
    // Check the neighbors at edges [b, c] and [c, a].
    for (j = 0; j < 2; j++) {
      enextself(midface); // [b, c] or [c, a].
      pg = apex(midface);
      toptet = midface;
      bflag = false;
      while (1) {
        // Go to the next face in the same tet.
        enext0fnextself(toptet);
        pc = apex(toptet);
        if (pinfected(pc)) {
          break; // Find a subface.
        }
        if (pc == dummypoint) {
          break; // Find a subface.
        }
        assert(pc != pg); // SELF_CHECK
        // Go to the same face in the adjacent tet.
        symedgeself(toptet);
        // Do we walk outside the cavity?
        if (!marktested(toptet)) {
          // Yes, the adjacent face is not a middle face.
          bflag = true; break; 
        }
      }
      if (!bflag) {
        symedge(midface, bottet);
        while (1) {
          assert(marktested(bottet)); // SELF_CHECK
          enext0fnextself(bottet);
          pf = apex(bottet);
          if (pf == pc) break; // Face matched.
          assert(!pinfected(pf)); // SELF_CHECK
          symedgeself(bottet);
        }
        // Connect two tets together.
        bond(toptet, bottet);
        // Both are interior tets.
        infect(toptet);
        infect(bottet);
        // Add this face into list.
        esymself(toptet);
        midfaces->newindex((void **) &parytet);
        *parytet = toptet;
      }
    } // j
  } // i

  if (b->verbose > 1) {
    printf("    Found %ld middle subfaces.\n", midfaces->objects);
  }
  if (midfaces->objects > maxregionsize) {
    maxregionsize = midfaces->objects;
  }

  // Bond subsegments to new tets.
  for (k = 0; k < 2; k++) {
    cavfaces = (k == 0 ? topfaces : botfaces);
    for (i = 0; i < cavfaces->objects; i++) {
      parytet = (triface *) fastlookup(cavfaces, i);
      // Bond a subsegment (if it exists).
      for (j = 0; j < 3; j++) {
        tsspivot(*parytet, checkseg);
        if (checkseg.sh != NULL) {
          symedge(*parytet, neightet);
          assert(marktested(neightet)); // SELF_CHECK
          while (1) {
            tssbond1(neightet, checkseg);
            fnextself(neightet);
            if (!marktested(neightet)) break;
          }
        }
        enextself(*parytet);
      }
    }
  }

  // Unmark all facet vertices.
  worksh = * (face *) fastlookup(misregion, 0);
  pa = sorg(worksh);
  pb = sdest(worksh);
  pc = sapex(worksh);
  puninfect(pa);
  puninfect(pb);
  puninfect(pc);
  for (i = 0; i < misregion->objects; i++) {
    worksh = * (face *) fastlookup(misregion, i);
    pf = sapex(worksh);
    puninfect(pf);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// carvecavity()    Delete old tets and outer new tets of the cavity.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::carvecavity(arraypool *crosstets, arraypool *topnewtets,
  arraypool *botnewtets)
{
  arraypool *newtets;
  triface *parytet, *pnewtet, neightet;
  int i, j, k;

  // Delete the old tets in cavity.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    tetrahedrondealloc(parytet->tet);
  }
  crosstets->restart(); // crosstets will be re-used.

  // Collect infected new tets in cavity.
  for (k = 0; k < 2; k++) {
    newtets = (k == 0 ? topnewtets : botnewtets);
    for (i = 0; i < newtets->objects; i++) {
      parytet = (triface *) fastlookup(newtets, i);
      if (infected(*parytet)) {
        crosstets->newindex((void **) &pnewtet);
        *pnewtet = *parytet;
      }
    }
  }
  // Collect all new tets in cavity.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    for (j = 0; j < 4; j++) {
      decode(parytet->tet[j], neightet);
      if (marktested(neightet)) { // Is it a new tet?
        if (!infected(neightet)) {
          // Find an interior tet.
          assert((point) neightet.tet[7] != dummypoint); // SELF_CHECK
          infect(neightet);
          crosstets->newindex((void **) &pnewtet);
          *pnewtet = neightet;
        }
      }
    }
  }

  // Delete outer new tets.
  for (k = 0; k < 2; k++) {
    newtets = (k == 0 ? topnewtets : botnewtets);
    for (i = 0; i < newtets->objects; i++) {
      parytet = (triface *) fastlookup(newtets, i);
      if (infected(*parytet)) {
        // This is an interior tet.
        uninfect(*parytet);
        unmarktest(*parytet);
      } else {
        // An outer tet. Delete it.
        tetrahedrondealloc(parytet->tet);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// constrainedfacets()    Recover subfaces saved in 'subfacestack'.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::constrainedfacets()
{
  arraypool *crosstets, *topnewtets, *botnewtets;
  arraypool *topfaces, *botfaces, *midfaces;
  arraypool *misregion;
  triface *parytet, searchtet;
  face *pssub, *pssub1, ssub;
  enum intersection dir;
  long bakflipcount, cavitycount;
  int bakhullsize;
  int s;

  if (b->verbose) {
    printf("  Constraining facets.\n");
  }

  // Initialize arrays.
  crosstets = new arraypool(sizeof(triface), 10);
  topnewtets = new arraypool(sizeof(triface), 10);
  botnewtets = new arraypool(sizeof(triface), 10);
  topfaces = new arraypool(sizeof(triface), 10);
  botfaces = new arraypool(sizeof(triface), 10);
  midfaces = new arraypool(sizeof(triface), 10);
  misregion = new arraypool(sizeof(face), 8);

  bakflipcount = flip22count;
  cavitycount = 0l;

  // Loop until 'subfacstack' is empty.
  while (subfacstack->objects > 0l) {
    // The list is used as a stack.
    subfacstack->objects--;
    pssub = (face *) fastlookup(subfacstack, subfacstack->objects);
    ssub = *pssub;

    if (!sinfected(ssub)) continue; // Not a missing subface.
    suninfect(ssub);

    // Insert the subface.
    searchtet.tet = NULL;
    dir = scoutsubface(&ssub, &searchtet);
    if (dir == SHAREFACE) continue;

    // Push the face back into stack.
    sinfect(ssub);
    if (dir != EDGETRIINT) {
      subfacstack->newindex((void **) pssub);
      *pssub = ssub;
    } else {
      s = randomnation(subfacstack->objects - 1);
      pssub = (face *) fastlookup(subfacstack, s);
      subfacstack->newindex((void **) &pssub1);
      *pssub1 = *pssub;
      *pssub = ssub;
      continue;
    }

    // Search for a crossing tet.
    misregion->newindex((void **) &pssub);
    *pssub = ssub;
    dir = scoutcrosstet(misregion, &searchtet);

    if (dir == ACROSSTET) {
      // The subface is missing. Recover it by local retetrahedralization.
      cavitycount++;
      bakhullsize = hullsize;
      checksubsegs = 0;
      crosstets->newindex((void **) &parytet);
      *parytet = searchtet;
      // Form a cavity of crossing tets.
      formcavity(misregion, crosstets, topfaces, botfaces);
      // Tetrahedralize the top part.
      delaunizecavity(topfaces, topnewtets);
      // Tetrahedralize the bottom part.
      delaunizecavity(botfaces, botnewtets);
      // Fill the cavity with new tets.
      fillcavity(topfaces, botfaces, midfaces, misregion);
      // Delete old tets and outer new tets.
      carvecavity(crosstets, topnewtets, botnewtets);
      // Clear working lists.
      crosstets->restart();
      topnewtets->restart();
      botnewtets->restart();
      topfaces->restart();
      botfaces->restart();
      midfaces->restart();
      hullsize = bakhullsize;
      checksubsegs = 1;
    } else if (dir == ACROSSFACE) {
      // Recover subface(s) by edge flips.
      adjustsurfmesh(misregion, &searchtet);
    } else {
      assert(0); // Not handled yet.
    }
    misregion->restart();
  }

  if (b->verbose) {
    printf("  %ld flips in surface mesh.\n", flip22count - bakflipcount);
    printf("  %ld cavities remeshed.\n", cavitycount);
  }

  // Delete arrays.
  delete crosstets;
  delete topnewtets;
  delete botnewtets;
  delete topfaces;
  delete botfaces;
  delete midfaces;
  delete misregion;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// formskeleton()    Form a constrained tetrahedralization.                  //
//                                                                           //
// The segments and facets of a PLS will be recivered.                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::formskeleton()
{
  face *pssub, ssub;
  face *psseg, sseg;
  REAL bakeps;
  int s, i;

  if (!b->quiet) {
    printf("Recovering boundaries.\n");
  }

  // Bakup the epsilon.
  bakeps = b->epsilon;
  b->epsilon = 0;

  // Construct a map from point to tets for speeding point location.
  makepoint2tetmap();

  // Mark acutes vertices.
  markacutevertices();

  // Put all segments into the list.
  if (b->order == 4) {  // '-o4' option (for debug)
    // The sequential order.
    subsegpool->traversalinit();
    for (i = 0; i < subsegpool->items; i++) {
      sseg.sh = shellfacetraverse(subsegpool);
      sinfect(sseg);  // Only save it once.
      subsegstack->newindex((void **) &psseg);
      *psseg = sseg;
    }
  } else {
    // Randomly order the segments.
    subsegpool->traversalinit();
    for (i = 0; i < subsegpool->items; i++) {
      s = randomnation(i + 1);
      // Move the s-th seg to the i-th.
      subsegstack->newindex((void **) &psseg);
      *psseg = * (face *) fastlookup(subsegstack, s);
      // Put i-th seg to be the s-th.
      sseg.sh = shellfacetraverse(subsegpool);
      sinfect(sseg);  // Only save it once.
      psseg = (face *) fastlookup(subsegstack, s);
      *psseg = sseg;
    }
  }

  // Segments will be introduced.
  checksubsegs = 1;
  // Recover segments.
  delaunizesegments();

  // Put all subfaces into list (in sequential order).
  subfacepool->traversalinit();
  for (i = 0; i < subfacepool->items; i++) {
    ssub.sh = shellfacetraverse(subfacepool);
    sinfect(ssub);  // Only save it once.
    subfacstack->newindex((void **) &pssub);
    *pssub = ssub;
  }

  // Subfaces will be introduced.
  checksubfaces = 1;
  // Recover facets.
  constrainedfacets();

  // checksubsegs = 0;
  b->epsilon = bakeps;
}

#endif // #ifndef constrainCXX
