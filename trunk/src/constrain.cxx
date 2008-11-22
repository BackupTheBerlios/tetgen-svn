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
    //   and d lies above the horizon. We test the orientations of the
    //   search point with respect to three planes: abc (horizon plnae),
    //   bad (right plane), and acd (left plane). 
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
      tsspivot1(*searchtet, checkseg);  // SELF_CHECK
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

  enextfnextself(*searchtet); // Go the opposite face.
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
    tsspivot1(*searchtet, checkseg);
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
        pd = oppo(neightet);
        // dir = tri_edge_inter(pa, pb, pc, startpt, endpt, pd, &pos);
        if (dir != DISJOINT) break;
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
        pd = oppo(neightet);
        // dir = tri_edge_inter(pa, pb, pc, startpt, endpt, pd, &pos);
        if (dir != DISJOINT) break;
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
      tsspivot1(*searchtet, checkseg);
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

  if (pointtype(ei) == ACUTEVERTEX) {
    if (pointtype(ej) == ACUTEVERTEX) {
      // Both ei and ej are ACUTEVERTEX.
      stype = 0;
    } else {
      // ej is either a RIDGEVERTEX or a STEINERVERTEX.
      stype = 1;
    }
  } else {
    if (pointtype(ei) == RIDGEVERTEX) {
      if (pointtype(ej) == ACUTEVERTEX) {
        stype = 1; sign = -1;
      } else {
        if (pointtype(ej) == RIDGEVERTEX) {
          // Both ei and ej are non-acute.
          stype = 0;
        } else {
          // ej is a STEINERVETEX.
          ek = farsdest(*sseg);
          if (pointtype(ek) == ACUTEVERTEX) {
            stype = 1; sign = -1;
          } else {
            stype = 0;
          }
        }
      }
    } else {
      // ei is a STEINERVERTEX.
      if (pointtype(ej) == ACUTEVERTEX) {
        stype = 1; sign = -1;
      } else {
        ek = farsorg(*sseg);
        if (pointtype(ej) == RIDGEVERTEX) {
          if (pointtype(ek) == ACUTEVERTEX) {
            stype = 1;
          } else {
            stype = 0;
          }
        } else {
          // Both ei and ej are STEINERVETEXs. ei has priority.
          if (pointtype(ek) == ACUTEVERTEX) {
            stype = 1;
          } else {
            ek = farsdest(*sseg);
            if (pointtype(ek) == ACUTEVERTEX) {
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
      d3 = DIST(ei, ej);
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
      pointtype(pa) = acuteflag ? ACUTEVERTEX : RIDGEVERTEX;
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
  REAL bakeps;
  int s, i;

  shellface sptr;

  if (!b->quiet) {
    printf("Delaunizing segments.\n");
  }

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
  // Bakup the epsilon.
  bakeps = b->epsilon;
  b->epsilon = 0;

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
        pointtype(newpt) = STEINERVERTEX;
        sspivot(splitshs[0], sseg);
        sspivot(splitshs[1], nsseg);
        // Some subfaces may be non-Delaunay.
        lawsonflip();
        // Insert newpt into the DT.
        insertvertex(newpt, &searchtet, true);
      } else {
        // Split the segment by refpt.
        flipn2nf(refpt, splitshs, 1);
        if (pointtype(refpt) != ACUTEVERTEX) {
          pointtype(refpt) = RIDGEVERTEX;
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

  checksubsegs = 0;
  b->epsilon = bakeps;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutsubface()    Look for a given subface in the tetrahedralization T.   //
//                                                                           //
// 'ssub' is the subface, denoted as abc. If abc exists in T, it is 'locked' //
// at the place where the two tets sharing at it.                            //
//                                                                           //
// The returned value indicates onre of the following cases:                 //
//   - COPLANAR, abc exists and is inserted;                                 //
//   - ACROSSVERT, a vertex (the origin of 'searchtet') lies on the facet.   //
//   - ACROSSEDGE, an edge (in 'searchtet') intersects ab;                   //
//   - ACROSSFACE, a face (in 'searchtet') intersects ab;                    //
//   - ACROSSTET, a tet (in 'searchtet') crosses the facet containg abc.     //
//                                                                           //
// If the returned value is ACROSSEDGE or ACROSSFACE, it means none of the   //
// three edges of abc exists in T.                                           //
//                                                                           //
// If the retunred value is ACROSSTET, let 'searchtet' be abde. The edge de  //
// intersects the facet containing abc. The vertex d lies exactly below the  //
// facet, while the vertex e may lie exactly on the facet, i.e., a, b, c,    //
// and e are coplanar.                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::intersection tetgenmesh::scoutsubface(face* ssub,
  triface* searchtet)
{
  triface spintet;
  face checksh;
  point pa, pb, pc, pd, pe;
  enum intersection dir;
  REAL ori;
  int i;

  tetrahedron ptr;
  int *iptr, tver;

  // Search an edge of 'ssub' in tetrahedralization.
  for (i = 0; i < 3; i++) {
    pa = (point) ssub->sh[3 + i];
    pb = (point) ssub->sh[3 + (i + 1) % 3];
    
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
      pd = dest(*searchtet); 
      if (pd == pb) {
        // Found the edge. Break the loop.
        break;
      } else {
        // A vertex lies on the search edge. Return it.
        enextself(*searchtet);
        return ACROSSVERT;
      }
    }
  }

  if (i == 3) {
    // None of the three edges exists.
    return dir;  // ACROSSEDGE or ACROSSFACE.
  }

  ssub->shver = (i << 1);
  pc = sapex(*ssub);

  if (b->verbose > 1) {
    printf("    Scout subface (%d, %d, %d).\n", pointmark(pa), pointmark(pb),
      pointmark(pc));
  }

  // Searchtet holds edge pa->pb. Search a face with apex pc.
  spintet = *searchtet;
  while (1) {
    fnextself(spintet);
    pd = apex(spintet);
    if (pd == pc) {
      // Found! Insert the subface.
      tspivot(spintet, checksh);
      assert(checksh.sh == NULL); // SELF_CHECK
      tsbond(spintet, *ssub);
      symedgeself(spintet);
      tspivot(spintet, checksh);
      assert(checksh.sh == NULL); // SELF_CHECK
      tsbond(spintet, *ssub);
      return SHAREFACE;
    }
    if (pd == apex(*searchtet)) break;
  }

  // Search an edge crossing the facet containing abc.
  if (searchtet->ver & 01) {
    // Adjust to 0th edge ring.
    esymself(*searchtet);
    sesymself(*ssub);
    pd = pa; pa = pb; pb = pd;
  }

  // Get a face containing ab and its apex lies below abc.
  spintet = *searchtet;
  while (1) {
    pd = apex(spintet);
    ori = orient3d(pa, pb, pc, pd);
    if (ori > 0) break;
    fnextself(spintet);
    assert(pd != apex(*searchtet)); // SELF_CHECK
  }
  // Search a tet whose apex->oppo crosses the facet containig abc.
  while (1) {
    pe = oppo(spintet);
    ori = orient3d(pa, pb, pc, pe);
    if (ori <= 0) break;  // stop at pd->pe.    
    fnextself(spintet);
  }

  if (ori == 0) {
    if (pointtype(pe) == VOLVERTEX) {
      // A vertex lies on the facet.
      enext0fnext(spintet, *searchtet);
      esymself(*searchtet);  // b->a.
      enext2self(*searchtet);  // e->a.
      return ACROSSVERT;
    }
  }

  // Return a crossing tet.
  if (b->verbose > 1) {
    printf("    Found a crossing tet (%d, %d, %d, %d).\n", pointmark(pa),
      pointmark(pb), pointmark(apex(spintet)), pointmark(pe));
  }

  *searchtet = spintet;
  return ACROSSTET;
}

#endif // #ifndef constrainCXX
