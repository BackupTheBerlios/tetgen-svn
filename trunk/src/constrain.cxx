#ifndef constrainCXX
#define constrainCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// finddirection()    Find the tet on the path from on point to another.     //
//                                                                           //
// The path starts from 'searchtet''s origin and ends at 'endpt'. On finish, //
// 'searchtet' returns the tet on the path, its origin does not change.      //
//                                                                           //
// The return value notes whether the path is collinear with an edge, or cr- //
// osses an edge, or crosses an face of the tet, 'searchtet' represents the  //
// according edge or face, respectively.                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::direction tetgenmesh::finddirection(triface* searchtet, 
  point endpt)
{
  triface neightet;
  point pa, pb, pc, pd, pn;
  REAL hori, rori, lori;
  enum {HMOVE, RMOVE, LMOVE} nextmove;
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

  if (b->verbose > 1) {
     printf("    Search direction from (%d, %d, %d) to %d.\n", 
       pointmark(pa), pointmark(pb), pointmark(pc), pointmark(endpt));
  }

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
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutsegment()    Look for a given segment in the tetrahedralization.     //
//                                                                           //
// Search an edge in current tetrahedralization that matches the given segm- //
// ment. Return TRUE if such edge exists, and the segment is attached to the //
// surrounding tets. Otherwise return FALSE.                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::scoutsegment(face* insseg)
{
  triface searchtet;
  point pa, pb;
  enum location loc;
  bool searchflag;
  int m, i;

  if (b->verbose > 1) {
    printf("    Scout segment (%d, %d).\n", pointmark(insseg->sh[3]), 
      pointmark(insseg->sh[4]));
  }

  searchflag = false;

  // Search a tet containing an endpoint of this segment as origin.
  for (m = 0; m < 2 && !searchflag; m++) {
    pa = (point) insseg->sh[m + 3];
    decode(point2tet(pa), searchtet);
    if (searchtet.tet[4] != NULL) {
      // Check if this tet contains pa.
      for (i = 4; i < 8 && !searchflag; i++) {
        if ((point) searchtet.tet[i] == pa) {
          searchflag = true;
          // Found. Set pa as its origin.
          switch (i) {
            case 4: searchtet.loc = 0; searchtet.ver = 0; break;
            case 5: searchtet.loc = 0; searchtet.ver = 2; break;
            case 6: searchtet.loc = 0; searchtet.ver = 4; break;
            case 7: searchtet.loc = 1; searchtet.ver = 2; break;
          }
          pb = (point) insseg->sh[(m + 1) % 2 + 3];
        }
      }
    }
  }

  if (!searchflag) {
    // Locate pa in tetrahedralization.
    pa = (point) insseg->sh[3];
    randomsample(pa, &searchtet);
    loc = locate(pa, &searchtet);
    assert(loc == ONVERTEX);  // SELF_CHECK
    pb = (point) insseg->sh[4];
    force_ptloc_count++;
  }


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
