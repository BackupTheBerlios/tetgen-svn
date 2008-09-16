#ifndef delaunayCXX
#define delaunayCXX

#include "tetgen.h"

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

  insphere_count++;

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
// randomnation()    Generate a random number between 0 and 'choices' - 1.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

unsigned long tetgenmesh::randomnation(unsigned long choices)
{
  unsigned long newrandom;

  if (choices >= 714025l) {
    newrandom = (randomseed * 1366l + 150889l) % 714025l;
    randomseed = (newrandom * 1366l + 150889l) % 714025l;
    newrandom = newrandom * (choices / 714025l) + randomseed;
    if (newrandom >= choices) {
      return newrandom - choices;
    } else {
      return newrandom;
    }
  } else {
    randomseed = (randomseed * 1366l + 150889l) % 714025l;
    return randomseed % choices;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// randomsample()    Randomly sample the tetrahedra for point loation.       //
//                                                                           //
// This routine implements Muecke's Jump-and-walk point location algorithm.  //
// It improves the simple walk-through by "jumping" to a good starting point //
// via random sampling.  Searching begins from one of handles:  the input    //
// 'searchtet', a recently encountered tetrahedron 'recenttet',  or from one //
// chosen from a random sample.  The choice is made by determining which one //
// 's origin is closest to the point we are searcing for.  Having chosen the //
// starting tetrahedron, the simple Walk-through algorithm is executed.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::randomsample(point searchpt, triface *searchtet)
{
  tetrahedron *firsttet, *tetptr;
  point torg;
  void **sampleblock;
  long sampleblocks, samplesperblock, samplenum;
  unsigned long alignptr;
  REAL searchdist, dist;
  int tetblocks, i, j;

  if ((searchtet->tet != NULL) && (searchtet->tet[4] != NULL)) {
    // Get the distance from the suggested starting tet to the search point.
    if ((point) searchtet->tet[7] != dummypoint) {
      torg = org(*searchtet);
    } else {
      torg = (point) searchtet->tet[4];
    }
    searchdist = NORM2(searchpt[0] - torg[0], searchpt[1] - torg[1], 
                       searchpt[2] - torg[2]);
  } else {
    searchdist = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
  }

  // If a recently encountered tetrahedron has been recorded and has not
  //   been deallocated, test it as a good starting point.
  if ((recenttet.tet != NULL) && (recenttet.tet[4] != NULL)) {
    if ((point) recenttet.tet[7] != dummypoint) {
      torg = org(recenttet);
    } else {
      torg = (point) recenttet.tet[4];
    }
    dist = NORM2(searchpt[0] - torg[0], searchpt[1] - torg[1],
                 searchpt[2] - torg[2]);
    if (dist < searchdist) {
      *searchtet = recenttet;
      searchdist = dist;
    }
  }

  // Select "good" candidate using k random samples, taking the closest one.
  //   The number of random samples taken is proportional to the fourth root
  //   of the number of tetrahedra in the mesh.
  while (samples * samples * samples * samples < tetrahedronpool->items) {
    samples++;
  }
  // Find how much blocks in current tet pool.
  tetblocks = (tetrahedronpool->maxitems + ELEPERBLOCK - 1) / ELEPERBLOCK;
  // Find the average samles per block. Each block at least have 1 sample.
  samplesperblock = (samples + tetblocks - 1) / tetblocks;
  sampleblocks = samples / samplesperblock;
  sampleblock = tetrahedronpool->firstblock;
  for (i = 0; i < sampleblocks; i++) {
    alignptr = (unsigned long) (sampleblock + 1);
    firsttet = (tetrahedron *)
               (alignptr + (unsigned long) tetrahedronpool->alignbytes
               - (alignptr % (unsigned long) tetrahedronpool->alignbytes));
    for (j = 0; j < samplesperblock; j++) {
      if (i == tetblocks - 1) {
        // This is the last block.
        samplenum = randomnation((int)
                      (tetrahedronpool->maxitems - (i * ELEPERBLOCK)));
      } else {
        samplenum = randomnation(ELEPERBLOCK);
      }
      tetptr = (tetrahedron *)
               (firsttet + (samplenum * tetrahedronpool->itemwords));
      if (tetptr[4] != (tetrahedron) NULL) {
        torg = (point) tetptr[4];
        dist = NORM2(searchpt[0] - torg[0], searchpt[1] - torg[1],
	             searchpt[2] - torg[2]);
        if (dist < searchdist) {
          searchtet->tet = tetptr;
          searchtet->loc = 0;
	  searchtet->ver = 0;
          searchdist = dist;
        }
      }
    }
    sampleblock = (void **) *sampleblock;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// locate()    Find a simplex containing a given point.                      //
//                                                                           //
// This routine implements the simple Walk-through point location algorithm. //
// Begins its search from 'searchtet', assume there is a line segment L from //
// the origin of 'searchtet' to the query point 'searchpt', and simply walk  //
// towards 'searchpt' by traversing all faces intersected by L.              //
//                                                                           //
// On completion, 'searchtet' is a tetrahedron that contains 'searchpt'. The //
// returned value indicates one of the following cases:                      //
//   - ONVERTEX, the search point lies on the origin of 'searchtet'.         //
//   - ONEDGE, the search point lies on an edge of 'searchtet'.              //
//   - ONFACE, the search point lies on a face of 'searchtet'.               //
//   - INTET, the search point lies in the interior of 'searchtet'.          //
//   - OUTSIDE, the search point lies outside the mesh. 'searchtet' is a     //
//     hull tetrahedron whose base face is visible by the search point.      //
//                                                                           //
// WARNING: This routine is designed for convex triangulations, and will not //
// generally work after the holes and concavities have been carved.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::locateresult tetgenmesh::locate(point searchpt,
  triface* searchtet)
{
  tetrahedron ptr;
  triface neightet;
  point torg, tdest, tapex, toppo, ntoppo;
  enum {ORGMOVE, DESTMOVE, APEXMOVE} nextmove;
  REAL ori, oriorg, oridest, oriapex;
  REAL searchdist, dist;
  int *iptr;

  if ((point) searchtet->tet[7] == dummypoint) {
    // A hull tet. Choose the neighbor of its base face.
    searchtet->loc = 0;
    symself(*searchtet);
  } else {
    // Keep the CCW edge ring.
    if (searchtet->ver &= 01) esymself(*searchtet);
  }
  // Let searchtet be the face such that 'searchpt' lies above to it.
  for (; ; searchtet->loc = (searchtet->loc + 1) % 4) { 
    torg = org(*searchtet);
    tdest = dest(*searchtet);
    tapex = apex(*searchtet);
    ori = orient3d(torg, tdest, tapex, searchpt);
    if (ori < 0) {
      // searchpt lies above searchtet's face.
      break;
    } else if (ori > 0) {
      // searchpt lies below searchtet's face.
      symself(*searchtet);
      torg = org(*searchtet);
      tdest = dest(*searchtet);
      tapex = apex(*searchtet);
      break;
    }
    // searchpt is coplanar with searchtet's face. Go to the next face.
  }

  ptloc_trav_tets_count = 1l;  // Algorithimic count.

  // Walk through tetrahedra to locate the point.
  while (true) {

    toppo = oppo(*searchtet);
    
    // Check if we have walked out of the domain.
    if (toppo == dummypoint) {
      return OUTSIDE;
    }
    
    // Check if the vertex is we seek.
    if (toppo == searchpt) {
      // Adjust the origin of searchtet to be searchpt.
      enext0fnextself(*searchtet);
      esymself(*searchtet);
      enext2self(*searchtet);
      return ONVERTEX;
    }

    // We enter from serarchtet's base face. There are three other faces in
    //   searchtet (all connecting to toppo), which one is the exit?
    oriorg = orient3d(tdest, tapex, toppo, searchpt); 
    oridest = orient3d(tapex, torg, toppo, searchpt);
    oriapex = orient3d(torg, tdest, toppo, searchpt);

    // Now decide which face to move. It is possible there are more than one
    //   faces are viable moves. Use the opposite points of thier neighbors
    //   to discriminate, i.e., we choose the face whose opposite point has
    //   the shortest distance to searchpt.
    if (oriorg < 0) {
      if (oridest < 0) {
        if (oriapex < 0) {
          // Any of the three faces is a viable move. 
          nextmove = ORGMOVE;
          enextfnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            searchdist = NORM2(searchpt[0] - ntoppo[0],
                               searchpt[1] - ntoppo[1],
                               searchpt[2] - ntoppo[2]);
          } else {
            searchdist = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
          }
          enext2fnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            dist = NORM2(searchpt[0] - ntoppo[0], searchpt[1] - ntoppo[1],
                         searchpt[2] - ntoppo[2]);
          } else {
            dist = searchdist;
          }
          if (dist < searchdist) {
            nextmove = DESTMOVE;
            searchdist = dist;
          }
          enext0fnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            dist = NORM2(searchpt[0] - ntoppo[0], searchpt[1] - ntoppo[1],
                         searchpt[2] - ntoppo[2]);
          } else {
            dist = searchdist;
          }
          if (dist < searchdist) {
            nextmove = APEXMOVE;
            searchdist = dist;
          }
        } else {
          // Two faces, opposite to origin and destination, are viable.
          nextmove = ORGMOVE;
          enextfnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            searchdist = NORM2(searchpt[0] - ntoppo[0],
                               searchpt[1] - ntoppo[1],
                               searchpt[2] - ntoppo[2]);
          } else {
            searchdist = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
          }
          enext2fnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            dist = NORM2(searchpt[0] - ntoppo[0], searchpt[1] - ntoppo[1],
                         searchpt[2] - ntoppo[2]);
          } else {
            dist = searchdist;
          }
          if (dist < searchdist) {
            nextmove = DESTMOVE;
            searchdist = dist;
          }
        }
      } else {
        if (oriapex < 0) {
          // Two faces, opposite to origin and apex, are viable.
          nextmove = ORGMOVE;
          enextfnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            searchdist = NORM2(searchpt[0] - ntoppo[0],
                               searchpt[1] - ntoppo[1],
                               searchpt[2] - ntoppo[2]);
          } else {
            searchdist = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
          }
          enext0fnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            dist = NORM2(searchpt[0] - ntoppo[0], searchpt[1] - ntoppo[1],
                         searchpt[2] - ntoppo[2]);
          } else {
            dist = searchdist;
          }
          if (dist < searchdist) {
            nextmove = APEXMOVE;
            searchdist = dist;
          }
        } else {
          // Only the face opposite to origin is viable.
          nextmove = ORGMOVE;
        }
      }
    } else {
      if (oridest < 0) {
        if (oriapex < 0) {
          // Two faces, opposite to destination and apex, are viable.
          nextmove = DESTMOVE;
          enext2fnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            searchdist = NORM2(searchpt[0] - ntoppo[0],
                               searchpt[1] - ntoppo[1],
                               searchpt[2] - ntoppo[2]);
          } else {
            searchdist = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
          }
          enext0fnext(*searchtet, neightet);
          symself(neightet);
          ntoppo = oppo(neightet);
          if (ntoppo != dummypoint) {
            dist = NORM2(searchpt[0] - ntoppo[0], searchpt[1] - ntoppo[1],
                         searchpt[2] - ntoppo[2]);
          } else {
            dist = searchdist;
          }
          if (dist < searchdist) {
            nextmove = APEXMOVE;
            searchdist = dist;
          }
        } else {
          // Only the face opposite to destination is viable.
          nextmove = DESTMOVE;
        }
      } else {
        if (oriapex < 0) {
          // Only the face opposite to apex is viable.
          nextmove = APEXMOVE;
        } else {
          // The point we seek must be on the boundary of or inside this
          //   tetrahedron. Check for boundary cases.
          if (oriorg == 0) {
            // Go to the face opposite to origin.
            enextfnextself(*searchtet);
            if (oridest == 0) {
              enextself(*searchtet); // edge apex->oppo
              return ONEDGE;
            }
            if (oriapex == 0) {
              enext2self(*searchtet);
              return ONEDGE;
            }
            return ONFACE;
          }
          if (oridest == 0) {
            // Go to the face opposite to destination.
            enext2fnextself(*searchtet);
            if (oriapex == 0) {
              enextself(*searchtet);
              return ONEDGE;
            }
            return ONFACE;
          }
          if (oriapex == 0) {
            // Go to the face opposite to apex
            enext0fnextself(*searchtet);
            return ONFACE;
          }
          return INTET;
        }
      }
    }
    
    // Move to the selected face.
    if (nextmove == ORGMOVE) {
      enextfnextself(*searchtet);
    } else if (nextmove == DESTMOVE) {
      enext2fnextself(*searchtet);
    } else {
      enext0fnextself(*searchtet);
    }
    // Move to the adjacent tetrahedron (maybe a hull tetrahedron).
    symself(*searchtet);
    // Retreat the three vertices of the base face.
    torg = org(*searchtet);
    tdest = dest(*searchtet);
    tapex = apex(*searchtet);

    ptloc_trav_tets_count++;  // Algorithimic count.

  } // while (true)
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initialDT()    Create an initial Delaunay tetrahedralization.             //
//                                                                           //
// The tetrahedralization contains only one tetrahedron abcd, and four hull  //
// tetrahedra.  The points pa, pb, pc, and pd must be linearly independent.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::initialDT(point pa, point pb, point pc, point pd)
{
  triface firsttet, tetopa, tetopb, tetopc, tetopd;
  triface worktet, worktet1;
  int *iptr;

  if (b->verbose > 1) {
    printf("    Create init tet (%d, %d, %d, %d)\n", pointmark(pa),
      pointmark(pb), pointmark(pc), pointmark(pd));
  }

  // Create the first tetrahedron.
  maketetrahedron(tetrahedronpool, &firsttet);
  setvertices(firsttet, pa, pb, pc, pd);
  // Create four hull tetrahedra.
  maketetrahedron(hulltetrahedronpool, &tetopa);
  setvertices(tetopa, pb, pc, pd, dummypoint);
  maketetrahedron(hulltetrahedronpool, &tetopb);
  setvertices(tetopb, pc, pa, pd, dummypoint);
  maketetrahedron(hulltetrahedronpool, &tetopc);
  setvertices(tetopc, pa, pb, pd, dummypoint);
  maketetrahedron(hulltetrahedronpool, &tetopd);
  setvertices(tetopd, pb, pa, pc, dummypoint);

  // Connect hull tetrahedra to firsttet (at four faces of firsttet).
  bond(firsttet, tetopd);
  enext0fnext(firsttet, worktet);
  bond(worktet, tetopc);
  enextfnext(firsttet, worktet);
  bond(worktet, tetopa);
  enext2fnext(firsttet, worktet);
  bond(worktet, tetopb);

  // Connect hull tetrahedra together (at six edges of firsttet).
  enext0fnext(tetopc, worktet); // edge ab
  enext0fnext(tetopd, worktet1);
  bond(worktet, worktet1);
  enext0fnext(tetopa, worktet); // edge bc
  enext2fnext(tetopd, worktet1);
  bond(worktet, worktet1);
  enext0fnext(tetopb, worktet); // edge ca
  enextfnext(tetopd, worktet1);
  bond(worktet, worktet1);
  enext2fnext(tetopc, worktet); // edge da
  enextfnext(tetopb, worktet1);
  bond(worktet, worktet1);
  enext2fnext(tetopa, worktet); // edge db
  enextfnext(tetopc, worktet1);
  bond(worktet, worktet1);
  enext2fnext(tetopb, worktet); // edge dc
  enextfnext(tetopa, worktet1);
  bond(worktet, worktet1);

  // Set the vertex type.
  if (pointtype(pa) == UNUSEDVERTEX) {
    pointtype(pa) = VOLVERTEX;
  }
  if (pointtype(pb) == UNUSEDVERTEX) {
    pointtype(pb) = VOLVERTEX;
  }
  if (pointtype(pc) == UNUSEDVERTEX) {
    pointtype(pc) = VOLVERTEX;
  }
  if (pointtype(pd) == UNUSEDVERTEX) {
    pointtype(pd) = VOLVERTEX;
  }

  // Remember the first tetrahedron.
  recenttet = firsttet;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertvertex()    Insert a point (p) into current tetrahedralization (T). //
//                                                                           //
// The point p will be first located in T. 'searchtet' is a suggested start- //
// tetrahedron, it can be NULL. Note that p may lies outside T. In such case,//
// the convex hull of T will be updated to include p as a vertex.            //
//                                                                           //
// If 'bowyerwatson' is TRUE, the Bowyer-Watson algorithm is used to recover //
// the Delaunayness of T. If 'bowyerwatson' is FALSE and 'incrflip' is TRUE, //
// then the incremental flip algorithm is used. If both of these two options //
// are FALSE, only p is inserted, do nothing with regard to the Delaunayness //
// of T (T may be non-Delaunay after this function).                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::insertvertex(point insertpt, triface *searchtet,
  bool bowyerwatson, bool increflip)
{
  list *cavetetlist, *cavebdrylist;
  tetrahedron ptr;
  triface *cavetet, neightet, newtet;
  point *pts;
  enum locateresult loc;
  REAL sign, ori;
  bool enqflag;
  int copcount;
  int *iptr, i, j;

  if (b->verbose > 1) {
    printf("    Insert point %d\n", pointmark(insertpt));
  }

  // Locate the point p.
  if (searchtet->tet == NULL) {
    randomsample(insertpt, searchtet);
  }
  loc = locate(insertpt, searchtet);

  // Update algorithmic counts.
  if (ptloc_max_tets_count < ptloc_trav_tets_count) {
    ptloc_max_tets_count = ptloc_trav_tets_count;
  }
  
  if (loc == ONVERTEX) {
    // The point already exists. Mark it and do nothing on it.
    point2ppt(insertpt) = (tetrahedron) org(*searchtet);
    pointtype(insertpt) = DUPLICATEDVERTEX;
    dupverts++;
    return;
  }

  if (b->verbose > 1) {
    printf("    Located (%d) tet (%d, %d, %d, %d).\n", (int) loc,
      pointmark(org(*searchtet)), pointmark(dest(*searchtet)), 
      pointmark(apex(*searchtet)), pointmark(oppo(*searchtet)));
  }

  // Initialize working lists.
  cavetetlist = new list(sizeof(triface));
  cavebdrylist = new list(sizeof(triface));

  // Add the searchtet into list.
  infect(*searchtet);
  cavetetlist->append(searchtet);
  if (loc == ONFACE) {
    // Add its neighbor tet into list.
    sym(*searchtet, neightet);
    infect(neightet);
    cavetetlist->append(&neightet);
  } else if (loc == ONEDGE) {
    // Add all tets sharing its edge into list.
    fnext(*searchtet, neightet);
    while (neightet.tet != searchtet->tet) {
      infect(neightet);
      cavetetlist->append(&neightet);
      fnextself(neightet);
    }
  }

  // There is no coplanar cavity boundary face yet.
  copcount = 0;
  
  // Form a star-shaped cavity with respect to the inserting point.
  for (i = 0; i < cavetetlist->len(); i++) {
    cavetet = (triface *) cavetetlist->get(i);
    for (cavetet->loc = 0; cavetet->loc < 4; cavetet->loc++) {
      sym(*cavetet, neightet);
      // Do check if it is not included in 'cavetetlist' yet.
      if (!infected(neightet)) {
        enqflag = false;
        pts = (point *) neightet.tet;
        if (pts[7] != dummypoint) {
          // A volume tet. Operate on it if it has not been tested yet.
          if (!marktested(neightet)) {
            if (bowyerwatson) {
              // Use Bowyer-Watson algorithm, do Delaunay check.
              sign = insphere_sos(pts[4], pts[5], pts[6], pts[7], insertpt);
              enqflag = (sign < 0.0);
            }
            marktest(neightet); // Only test it once.
          }
        } else {
          // It is a hull tet. Check if its base face is visible by p. 
          //   This happens when p lies lies outside the hull face.
          ori = orient3d(pts[4], pts[5], pts[6], insertpt);
          enqflag = (ori < 0.0);
          // Check if this face is coplanar with p. This case may create
          //   a degenerate tet (zero volume). 
          // Note: for convex domain, it only happen at hull face.
          if (ori == 0.0) copcount++;
        }
        if (enqflag) {
          // Found a tet in the cavity.
          infect(neightet);
          cavetetlist->append(&neightet);
        } else {
          // Found a boubdary face of the cavity. It may be a face of a hull
          //   tet which contains 'dummypoint'. Choose the edge in the face 
          //   such that its endpoints are not 'dummypoint', while its apex
          //   may be 'dummypoint' (see Fig. 1.4).
          cavetet->ver = 4;
          cavebdrylist->append(cavetet);
        }
      } // if (!infected(neightet)
    } // for (cavetet->loc = 0;
  } // for (i = 0;
  
  // Create new tetrahedra in the Bowyer-Watson cavity. Connect them to the
  //   tetrahedra at outside of the cavity.
  for (i = 0; i < cavebdrylist->len(); i++) {
    cavetet = (triface *) cavebdrylist->get(i);
    sym(*cavetet, neightet); // neightet may be a hull tet.
    if (apex(*cavetet) != dummypoint) {
      // Create a new tet in the cavity (see Fig. bowyerwatson 1 or 3).
      unmarktest(neightet); // Unmark it.
      neightet.ver = 0; // Choose the 0th edge ring.
      maketetrahedron(tetrahedronpool, &newtet);
      setorg(newtet, dest(neightet));
      setdest(newtet, org(neightet));
      setapex(newtet, apex(neightet));
      setoppo(newtet, insertpt);
    } else {
      // Create a new hull tet (see Fig. bowyerwatson 2).
      maketetrahedron(hulltetrahedronpool, &newtet);
      assert(cavetet->ver == 4); // SELF_CHECK (see Fig. 1.4).
      setorg(newtet, dest(*cavetet));
      setdest(newtet, org(*cavetet));
      setapex(newtet, insertpt);
      setoppo(newtet, dummypoint);
      // Note: the acvity boundary face is at the enext0fnext place.
      enext0fnextself(newtet);
    }
    // Connect newtet <==> neightet, this also disconnect the old bond at
    //   neightet (while cavetet still holds a pointer to neightet).
    bond(newtet, neightet);
    // Replace the old boundary face with the new tet in list.
    *cavetet = newtet;
  }
  
  // Connect the set of new tetrahedra together.
  for (i = 0; i < cavebdrylist->len(); i++) {
    cavetet = (triface *) cavebdrylist->get(i);
    cavetet->ver = 0;
    for (j = 0; j < 3; j++) {
      enext0fnext(*cavetet, newtet); // Go to the face.
      // Operate on it if it is open.
      if (newtet.tet[newtet.loc] == NULL) {
        // Find its adjacent face by rotating faces around the edge of
        //   cavetet. The rotating direction is opposite to newtet.
        //   Stop the rotate at a face which is open.
        esym(*cavetet, neightet); // Set the rotate dir.
        do {
          fnextself(neightet); // Go to the face in the adjacent tet. 
        } while (neightet.tet[neightet.loc] != NULL);
        bond(newtet, neightet); // Connect newtet <==> neightet.
      }
      enextself(*cavetet);
    }
  }

  // Delete the old tetrahedra of the cavity.
  for (i = 0; i < cavetetlist->len(); i++) {
    cavetet = (triface *) cavetetlist->get(i);
    tetrahedrondealloc(cavetet->tet);
  }
  
  if (bowyerwatson && (copcount > 0)) {
    // There may exist zero volume tetrahedra. Check and remove them.
    bowyerwatsonpostproc(cavebdrylist);
  }

  // If the flip option is used.
  if (!bowyerwatson && increflip) {
    lawsonflip(cavebdrylist);
  }

  // Set the point type.
  if (pointtype(insertpt) == UNUSEDVERTEX) {
    pointtype(insertpt) = VOLVERTEX;
  }

  // Set a handle for the point location.
  if (bowyerwatson) {
    // Soome new tets may be dead, find a live one.
    for (i = 0; i < cavebdrylist->len(); i++) {
      cavetet = (triface *) cavebdrylist->get(i);
      if (cavetet->tet[4] != NULL) break;
    }
    recenttet = * cavetet;
  }
  
  delete cavetetlist;
  delete cavebdrylist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// bowyerwatsonpostproc()    Remove degenerate tets from the cavity.         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::bowyerwatsonpostproc(list *cavebdrylist)
{
  queue *removeque;
  tetrahedron ptr;
  triface oldtets[4], newtets[2];
  triface *cavetet, neightet;
  point *pts;
  REAL ori;
  int *iptr, i, j;

  removeque = new queue(sizeof(triface));
  
  // Queue all degenerate tets.
  for (i = 0; i < cavebdrylist->len(); i++) {
    cavetet = (triface *) cavebdrylist->get(i);
    pts = (point *) cavetet->tet;
    if (pts[7] != dummypoint) {
      // Check if the new tet is degenerate.
      ori = orient3d(pts[4], pts[5], pts[6], pts[7]);
      if (ori == 0) {
        removeque->push(cavetet);
      }
    }
  }

  if (b->verbose > 1) {
    printf("    Removing %ld degenerate tets.\n", removeque->len());
  }
  
  while (!removeque->empty()) {
    cavetet = (triface *) removeque->pop();
    if (b->verbose > 1) {
      pts = (point *) cavetet->tet;
      printf("    Remove tet (%d, %d, %d, %d).\n", pointmark(pts[4]),
        pointmark(pts[5]), pointmark(pts[6]), pointmark(pts[7]));
    }
    // Find the hull edge in cavetet.
    cavetet->ver = 0;
    for (j = 0; j < 3; j++) {
      enext0fnext(*cavetet, neightet);
      symself(neightet);
      if ((point) neightet.tet[7] == dummypoint) break;
      enextself(*cavetet);
    }
    // Because of existing multiple degenerate cases. It is possible
    //   that the other hull face is not pop yet.
    if (j < 3) {
      // Collect tets for flipping the edge.
      oldtets[0] = *cavetet;
      for (j = 0; j < 3; j++) {
        fnext(oldtets[j], oldtets[j + 1]);
      }
      if (oldtets[3].tet != cavetet->tet) {
        printf("Internal error in insertvertex(): Unknown flip case.\n");
        terminatetetgen(1);
      }
      // Do a 3-to-2 flip to remove the degenerate tet.
      flip32(oldtets, newtets, NULL);
      // Delete the old tets.
      tetrahedrondealloc(oldtets[0].tet);
      tetrahedrondealloc(oldtets[1].tet);
      tetrahedrondealloc(oldtets[2].tet);
    } else { 
      removeque->push(cavetet);  // Wait for the next round.
    } // if (j < 3)
  }

  delete removeque;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lawsonflip()    Flip non-locally Delaunay faces by primitive flips.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::lawsonflip(list *cavebdrylist)
{
  queue *flipque, *tmpque;
  triface oldtets[128], newtets[128];
  triface *fliptet, neightet, *tmptet;
  point *pts, pd, pe;
  REAL sign, ori;
  long flipcount;
  bool success;
  int n, i, j;

  flipque = new queue(sizeof(triface));
  tmpque = new queue(sizeof(triface)); // For flipnm().
  flipcount = flip23count + flip32count;

  // Add all boundary faces in queue.
  for (i = 0; i < cavebdrylist->len(); i++) {
    fliptet = (triface *) cavebdrylist->get(i);
    markface(*fliptet);
    flipque->push(fliptet);
  }
  // Save a tet for point location.
  recenttet = *fliptet;

  if (b->verbose > 1) {
    printf("    Lawson flip %ld faces.\n", flipque->len());
  }

  // Loop until the queue is empty.
  while (!flipque->empty()) {
  
    fliptet = (triface *) flipque->pop();
    
    if (fliptet->tet[4] == NULL) continue; // A dead tet.
    assert(facemarked(*fliptet)); // SELF_CHECK
    unmarkface(*fliptet);
    fliptet->ver = 4;
    if (apex(*fliptet) == dummypoint) continue; // A hull face.
    pts = (point *) fliptet->tet;
    if (pts[7] == dummypoint) continue; // A hull tet.
    sym(*fliptet, neightet);
    pe = oppo(neightet); 
    if (pe == dummypoint) continue; // A hull tet.
    
    // Delaunay test.
    sign = insphere_sos(pts[4], pts[5], pts[6], pts[7], pe);
    
    if (sign < 0.0) {
      // Try to flip the face.
      pd = oppo(*fliptet);
      for (i = 0; i < 3; i++) {
        ori = orient3d(org(*fliptet), dest(*fliptet), pd, pe);
        if (ori <= 0) break;
        enextself(*fliptet);
      }
      if (i == 3) {
        // A 2-to-3 flip is found.
        neightet.ver = 0;
        for (j = 0; j < 3; j++) { // Find the same edge.
          if (org(neightet) == dest(*fliptet)) break;
          enextself(neightet);
        }
        assert(j < 3); // SELF_CHECK
        oldtets[0] = *fliptet;
        oldtets[1] = neightet;
        flip23(oldtets, newtets, flipque);
        tetrahedrondealloc(oldtets[0].tet);
        tetrahedrondealloc(oldtets[1].tet);
        recenttet = newtets[0]; // for point location.
      } else {
        // Try to flip the edge org->dest of 'fliptet'.
        n = 0;
        oldtets[n] = *fliptet;
        do {
          fnext(oldtets[n], oldtets[n + 1]);
          n++;
        } while (oldtets[n].tet != fliptet->tet);
        // n is the total number of tets.
        if (n == 3) {
          flip32(oldtets, newtets, flipque);
          tetrahedrondealloc(oldtets[0].tet);
          tetrahedrondealloc(oldtets[1].tet);
          tetrahedrondealloc(oldtets[2].tet);
        } else {
          success = flipnm(n, oldtets, newtets, tmpque);
          if (success) {
            // Put all faces in tmpque into flipque.
            while (!tmpque->empty()) {
              tmptet = (triface *) tmpque->pop();
              if (tmptet->tet[4] != NULL) {
                flipque->push(tmptet);
              }
            }
            for (j = 0; j < n; j++) {
              tetrahedrondealloc(oldtets[j].tet);
            }
            recenttet = newtets[0]; // for point location.
          } else {
            // Clear all faces queued in tmpque.
            while (!tmpque->empty()) {
              tmptet = (triface *) tmpque->pop();
              if (tmptet->tet[4] != NULL) {
                unmarkface(*tmptet);
              }
            }
            markface(oldtets[0]);
            flipque->push(&(oldtets[0]));
            recenttet = oldtets[0]; // for point location.
          }
          tmpque->clear();
        } // if (n == 3)
      } // if (i == 3)
    } // if (sign < 0.0)
  } // while

  if (b->verbose > 1) {
    printf("    %ld flips.\n", flip23count + flip32count - flipcount);
  }

  delete flipque;
  delete tmpque;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// incrementaldelaunay()    Form a Delaunay tetrahedralization by increment- //
//                          ally inserting vertices.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::incrementaldelaunay()
{
  triface searchtet;
  point *permutarray, swapvertex;
  REAL v1[3], v2[3], n[3];
  REAL bboxsize, bboxsize2, bboxsize3, ori;
  int randindex, i, j;

  // Form a random permuation (uniformly at random) of the set of vertices.
  permutarray = new point[in->numberofpoints];
  pointpool->traversalinit();
  if (b->nojettison) { // '-J' option (only for debug)
    for (i = 0; i < in->numberofpoints; i++) {
      permutarray[i] = (point) pointpool->traverse();
    }
  } else {
    for (i = 0; i < in->numberofpoints; i++) {
      randindex = randomnation(i + 1);
      permutarray[i] = permutarray[randindex];
      permutarray[randindex] = (point) pointpool->traverse();
    }
  }

  // Calculate the diagonal size of its bounding box.
  bboxsize = sqrt(NORM2(xmax - xmin, ymax - ymin, zmax - zmin));
  bboxsize2 = bboxsize * bboxsize;
  bboxsize3 = bboxsize2 * bboxsize;
  
  // Make sure the second vertex is not identical with the first one.
  i = 1;
  while ((DIST(permutarray[0], permutarray[i]) / bboxsize) < b->epsilon) {
    i++;
    if (i == in->numberofpoints - 1) {
      printf("Exception:  All vertices are (nearly) identical (Tol = %g).\n",
             b->epsilon);
      terminatetetgen(1);
    }
  }
  if (i > 1) {
    // Swap to move the non-indetical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[1];
    permutarray[1] = swapvertex;
  }

  // Make sure the third vertex is not collinear with the first two.
  i = 2;
  for (j = 0; j < 3; j++) {
    v1[j] = permutarray[1][j] - permutarray[0][j];
    v2[j] = permutarray[i][j] - permutarray[0][j];
  }
  CROSS(v1, v2, n);
  while ((sqrt(NORM2(n[0], n[1], n[2])) / bboxsize2) < b->epsilon) {
    i++;
    if (i == in->numberofpoints - 1) {
      printf("Exception:  All vertices are (nearly) collinear (Tol = %g).\n",
             b->epsilon);
      terminatetetgen(1);
    }
    for (j = 0; j < 3; j++) {
      v2[j] = permutarray[i][j] - permutarray[0][j];
    }
    CROSS(v1, v2, n);
  }
  if (i > 2) {
    // Swap to move the non-indetical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[2];
    permutarray[2] = swapvertex;
  }

  // Make sure the fourth vertex is not coplanar with the first three.
  i = 3;
  ori = orient3d(permutarray[0], permutarray[1], permutarray[2], 
                 permutarray[i]);
  while ((fabs(ori) / bboxsize3) < b->epsilon) {
    i++;
    if (i == in->numberofpoints - 1) {
      printf("Exception:  All vertices are (nearly) coplanar (Tol = %g).\n",
             b->epsilon);
      terminatetetgen(1);
    }
    ori = orient3d(permutarray[0], permutarray[1], permutarray[2], 
                   permutarray[i]);
  }
  if (i > 3) {
    // Swap to move the non-indetical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[3];
    permutarray[3] = swapvertex;
  }

  // Orient the first four vertices in permutarray so that they follow the
  //   right-hand rule.
  if (ori > 0.0) {
    // Swap the first two vertices.
    swapvertex = permutarray[0];
    permutarray[0] = permutarray[1];
    permutarray[1] = swapvertex;
  }

  // Create the initial Delaunay tetrahedralization.
  initialDT(permutarray[0], permutarray[1], permutarray[2], permutarray[3]);

  if (b->verbose) {
    printf("  Incremental inserting vertices.\n");
  }

  for (i = 4; i < in->numberofpoints; i++) {
    searchtet.tet = NULL;  // Randomly sample tetrahedra.
    insertvertex(permutarray[i], &searchtet, (b->bowyerwatson > 0), true);
  }

  delete [] permutarray;
}

#endif // #ifndef delaunayCXX