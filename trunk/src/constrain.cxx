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
  int *iptr;

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
              // pa->'endpt' is collinear with pa->pb.
              return COLLINEAR;
            }
            if (lori == 0) {
              // pa->'endpt' is collinear with pa->pc.
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
              // pa->'endpt' is collinear with pa->pd.
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

    // Move to the next tet.
    if (nextmove = RMOVE) {
      enext0fnextself(*searchtet);
    } else if (nextmove = LMOVE) {
      enext2fnextself(*searchtet);
    }
    symself(*searchtet);
    pa = org(*searchtet);
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
      // pa->'endpt' is collinear with pa->pb.
      return COLLINEAR;
    }
    if (lori == 0) {
      // pa->'endpt' is collinear with pa->pc.
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
      // pa->'endpt' is collinear with pa->pd.
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
// If 'searchtet' != NULL, it's origin must be the origin of 'sseg'. It //
// can be used as the starting tet for searching the edge. If such edge does //
// not found (return FALSE), 'searchtet' contains the tet whose edge or face //
// is crossed by the segment.                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::scoutsegment(face* sseg, triface* searchtet, point* refpt)
{
  triface neightet;
  face checkseg;
  point startpt, endpt;
  point pa, pb, pc, pd;
  enum location loc;
  enum intersection dir;
  REAL angmax, ang;
  bool orgflag;
  int shver, i;

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
      if (searchtet->tet[4] != NULL) {
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
    if (pd != endpt) {
      // Split the segment.
      flipn2nf(pd, sseg, 1);
      lawsonflip();
    }
    // Found! Insert (the first part of) the segment.
    neightet = *searchtet;
    do {
      tssbond1(neightet, *sseg);
      fnextself(neightet);
    } while (neightet.tet != searchtet->tet);
    if (pd != endpt) {
      // Search the scond part of the segment.
      senext2self(*sseg);
      spivotself(*sseg);
      sseg->shver = 0;
      if (sorg(*sseg) != pd) sesymself(*sseg);
      enextself(*searchtet);
      return scoutsegment(sseg, searchtet, refpt);
    } else {
      // The job is done.
      return true;
    }
  }

  if (b->verbose > 1) {
    printf("    Scout ref point of seg (%d, %d).\n", pointmark(startpt), 
      pointmark(endpt));
  }

  // The segment is missing. Search a reference point for it.
  symedgeself(*searchtet);
  assert(oppo(*searchtet) != dummypoint); // SELF_CHECK

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

  if (dir == ACROSSEDGE) {
    // Check whether two segments are intersecting.
    tsspivot1(*searchtet, checkseg);
    if (checkseg.sh != NULL) {
      printf("  Invalid PLC!  Two segments intersect each other.\n");
      printf("    1st: (%d, %d), 2nd: (%d, %d).\n", 
        pointmark(sorg(*sseg)), pointmark(sdest(*sseg)), 
        pointmark(org(*searchtet)), pointmark(dest(*searchtet)));
      terminatetetgen(1);
    }
    // Here we should check all apexes. NOT done yet.
  }

  return false;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// delaunizesegments()    Recover segments in a Delaunay tetrahedralization. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
/*
void tetgenmesh::delaunizesegments()
{


  if (!b->quiet) {
    printf("Delaunizing segments.\n");
  }

  // Segments will be attached to tets. Initialize the connection pool.
  tet2subpool = new memorypool(6*sizeof(shellface), SUBPERBLOCK, POINTER, 0);



  delete tet2subpool;
}*/

#endif // #ifndef constrainCXX
