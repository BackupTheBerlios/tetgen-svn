#ifndef constrainCXX
#define constrainCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// finddirection()    Find the tet on the path from one point to another.    //
//                                                                           //
// The path starts from 'searchtet''s origin and ends at 'endpt'. On finish, //
// 'searchtet' returns the tet on the path, its origin does not change.      //
//                                                                           //
// The return value notes whether the path is collinear with an edge, or cr- //
// osses an edge, or crosses an face of the tet, 'searchtet' represents the  //
// according edge or face, respectively.                                     //
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
    return COLLINEAR;
  }
  if (pc == endpt) {
    // pa->pc is the search edge.
    enext2self(*searchtet);
    esymself(*searchtet);
    return COLLINEAR;
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
      return COLLINEAR;
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
              return COLLINEAR;
            }
            if (lori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pc.
              enext2self(*searchtet);
              esymself(*searchtet);
              return COLLINEAR;
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
              return COLLINEAR;
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
      return COLLINEAR;
    }
    if (lori == 0) {
      // pa->'endpt' is COLLINEAR with pa->pc.
      enext2self(*searchtet);
      esymself(*searchtet);
      return COLLINEAR;
    }
    // pa->'endpt' crosses the edge pb->pc.
    enextself(*searchtet);
    return ACROSSEDGE;
  }
  if (rori == 0) {
    if (lori == 0) {
      // pa->'endpt' is COLLINEAR with pa->pd.
      enext0fnextself(*searchtet); // face abd.
      enext2self(*searchtet);
      esymself(*searchtet);
      return COLLINEAR;
    }
    // pa->'endpt' crosses the edge pb->pd.
    enext0fnextself(*searchtet); // face abd.
    enextself(*searchtet);
    return ACROSSEDGE;
  }
  if (lori == 0) {
    // pa->'endpt' crosses the edge pc->pd.
    enext2fnextself(*searchtet);  // face cad
    enext2self(*searchtet);
    return ACROSSEDGE;
  }
  // pa->'endpt' crosses the face bcd.
  enextfnextself(*searchtet);
  return ACROSSFACE;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutsegment()    Look for a given segment in the tetrahedralization.     //
//                                                                           //
// Search an edge in current tetrahedralization that matches the given segm- //
// ment. Return TRUE if such edge exists, and the segment is attached to the //
// surrounding tets. Otherwise return FALSE.                                 //
//                                                                           //
// If 'searchtet' != NULL, it's origin must be the origin of 'sseg'.  It is  //
// used as the starting tet for searching the edge.                          //
//                                                                           //
// If the segment does not exist (return FALSE), 'refpt' returns a reference //
// point for the segment. 'searchtet' returns the tet containing the 'refpt'.//
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
  REAL ori1, ori2;
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
            orgflag = true;
            // Found. Set pa as its origin.
            switch (i) {
              case 4: searchtet->loc = 0; searchtet->ver = 0; break;
              case 5: searchtet->loc = 0; searchtet->ver = 2; break;
              case 6: searchtet->loc = 0; searchtet->ver = 4; break;
              case 7: searchtet->loc = 1; searchtet->ver = 2; break;
            }
            sseg->shver = shver;
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

  if (dir == COLLINEAR) {
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
      return dir;
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

  // The segment is missing. Search a reference point for it.
  symedgeself(*searchtet);
  if (oppo(*searchtet) == dummypoint) {
    // Enter outside! Move back to inside.
    esymself(*searchtet);
    symedgeself(*searchtet); // This is the old one.
    fnextself(*searchtet);
    assert(oppo(*searchtet) != dummypoint);  // SELF_CHECK
  }

  pa = org(*searchtet);
  angmax = interiorangle(pa, startpt, endpt, NULL);
  *refpt = pa;
  pb = dest(*searchtet);
  ang = interiorangle(pb, startpt, endpt, NULL);
  if (ang > angmax) {
    angmax = ang;
    *refpt = pb;
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
        dir = tri_edge_inter(pa, pb, pc, startpt, endpt, pd, &pos);
        if (dir != DISJOINT) break;
      }
      assert(i < 3);  // SELF_CHECK
    } else { // if (dir == ACROSSEDGE)
      // Check whether two segments are intersecting.
      tsspivot1(*searchtet, checkseg);
      if (checkseg.sh != NULL) {
        printf("  Invalid PLC!  Two segments intersect each other.\n");
        printf("    1st: (%d, %d), 2nd: (%d, %d).\n", 
          pointmark(sorg(*sseg)), pointmark(sdest(*sseg)), 
          pointmark(org(*searchtet)), pointmark(dest(*searchtet)));
        terminatetetgen(1);
      }
      // Find the tet containing the face intersected by the segment.
      ori1 = orient3d(pa, pb, pc, endpt);
      assert(ori1 != 0);  // SELF_CHECK
      neightet = *searchtet;
      while (1) {
        fnextself(neightet);
        pc = apex(neightet);
        if (pc == endpt) break;  // Reach the end point.
        if (pc != dummypoint) {
          ang = interiorangle(pc, startpt, endpt, NULL);
          if (ang > angmax) {
            angmax = ang;
            *refpt = pc;
            reftet = neightet;
          }
        }
        // Check if this face intersects the segment.
        ori2 = orient3d(pa, pb, pc, endpt);
        if ((ori2 == 0) || (ori1 * ori2 < 0)) break;
        ori1 = ori2;
      }
      assert(neightet.tet != searchtet->tet); // SELF_CHECK
      // Stop if we reached the end point.
      if (pc == endpt) break;
      // Where we are?
      pd = oppo(neightet);
      if (pd == dummypoint) {
        // We entered outside! Move back.
        esymself(neightet);
        symedgeself(neightet);
      }
      // One of the two side faces in 'searchtet' intersects the segment.
      *searchtet = neightet;
      neightet.ver = 0;
      for (i = 0; i < 2; i++) {
        neightet.loc = locverpivot[searchtet->loc][searchtet->ver][i];
        pa = org(neightet);
        pb = dest(neightet);
        pc = apex(neightet);
        pd = oppo(neightet);
        dir = tri_edge_inter(pa, pb, pc, startpt, endpt, pd, &pos);
        if (dir != DISJOINT) break;
      }
      assert(i < 2);  // SELF_CHECK
    }

    if (dir == ACROSSVERT) {
      // This segment passing a vertex. Choose it and return.
      for (i = 0; i < pos; i++) {
        enextself(neightet);
      }
      pd = org(neightet);
      if (b->verbose > 2) {
        angmax =  interiorangle(pd, startpt, endpt, NULL);
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
    // Update the pa, pb, pc.
    pa = org(*searchtet);
    pb = dest(*searchtet);
    pc = apex(*searchtet);

  } // while (1)

  // dir is either ACROSSVERT, or ACROSSEDGE, or ACROSSFACE.
  if (b->verbose > 2) {
    printf("      Choose refpoint %d, ang(%g), dir(%d).\n", pointmark(*refpt),
      angmax / PI * 180.0, (int) dir);
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
  REAL split, L, d, d1, d2;
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
    d = DIST(ei, refpt);
    if (d < 0.5 * L) {
      // Choose ei as center.
      split = d / L;
      for (i = 0; i < 3; i++) {
        vt[i] = ei[i] + split * (ej[i] - ei[i]);
      }
    } else {
      // Choose ej as center.
      d = DIST(ej, refpt);
      split = d / L;
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
      if (d1 < 0.5 * d) {
        d -= d1;
      } else {
        d /= 2.0;
      }
      split = d / L;
      for (i = 0; i < 3; i++) {
        vt[i] = ek[i] + split * (ej[i] - ek[i]);
      }
    }
    d1 > d2 ? r3count++ : r2count++;
  }

  if (b->verbose > 1) {
    printf("    split = %g.\n", split);
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

void tetgenmesh::markacutevertices(int* idx2seglist, face* segperverlist)
{
  point pa, pb, pc;
  REAL anglimit, ang;
  bool acuteflag;
  int acutecount;
  int idx, i, j;

  if (b->verbose) {
    printf("  Marking acute vertices.\n");
  }

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
  int s, i;

  face *segperverlist;
  int *idx2seglist;

  shellface sptr;

  if (!b->quiet) {
    printf("Delaunizing segments.\n");
  }

  // Construct a map from point to tets for speeding point location.
  makepoint2tetmap();

  // Construct a map from points to segments.
  makepoint2submap(subsegpool, idx2seglist, segperverlist);

  // Mark acutes vertices.
  markacutevertices(idx2seglist, segperverlist);

  // Initialize the pool for tet-subseg connections.
  tet2subpool = new memorypool(6*sizeof(shellface), SUBPERBLOCK, POINTER, 0);

  // Calculate the log-2 of the number of segments.
  s = 1; i = 0;
  while (s < subsegpool->items) {
    s <<= 1; i++;
  }
  // Initialize the stack of segments to be recovered.
  subsegstack = new arraypool(sizeof(face), i + 1);
  // Put all segments into the list.
  subsegpool->traversalinit();
  for (i = 0; i < subsegpool->items; i++) {
    sseg.sh = shellfacetraverse(subsegpool);
    sinfect(sseg);  // Only save it once.
    subsegstack->newindex((void **) &psseg);
    *psseg = sseg;
  }

  delete [] idx2seglist;
  delete [] segperverlist;

  // Segments will be introduced.
  checksubsegs = 1;

  // Loop until 'subsegstack' is empty.
  while (subsegstack->objects > 0l) {
    // Default seglist is used as a stack.
    subsegstack->objects--;
    psseg = (face *) fastlookup(subsegstack, subsegstack->objects);
    sseg = *psseg;
    suninfect(sseg);

    // Insert the segment.
    searchtet.tet = NULL;
    dir = scoutsegment(&sseg, &searchtet, &refpt);

    if (dir != COLLINEAR) {
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
      sinfect(sseg);
      subsegstack->newindex((void **) &psseg);
      *psseg = sseg;
      sinfect(nsseg);
      subsegstack->newindex((void **) &psseg);
      *psseg = nsseg;
    }
  }

  if (b->verbose) {
    printf("  %d protecting points.\n", r1count + r2count + r3count);
  }

  // Clear the pointers from tets to subsegs.
  tetrahedronpool->traversalinit();
  searchtet.tet = tetrahedrontraverse(tetrahedronpool);
  while (searchtet.tet != NULL) {
    searchtet.tet[8] = NULL;
    searchtet.tet = tetrahedrontraverse(tetrahedronpool);  
  }

  hulltetrahedronpool->traversalinit();
  searchtet.tet = tetrahedrontraverse(hulltetrahedronpool);
  while (searchtet.tet != NULL) {
    searchtet.tet[8] = NULL;
    searchtet.tet = tetrahedrontraverse(hulltetrahedronpool);  
  }

  checksubsegs = 0;

  delete subsegstack;
  delete tet2subpool;
}

#endif // #ifndef constrainCXX
