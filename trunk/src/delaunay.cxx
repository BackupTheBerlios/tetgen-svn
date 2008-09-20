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
      } else {
        j--;
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
    ori = orient3d(torg, tdest, tapex, searchpt); orient3dcount++;
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

  // Walk through tetrahedra to locate the point.
  while (true) {

    ptloc_count++;  // Algorithimic count.

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
    orient3dcount+=3;

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
  bond(firsttet, tetopd); // ab
  enext0fnext(firsttet, worktet);
  bond(worktet, tetopc); // ab
  enextfnext(firsttet, worktet);
  bond(worktet, tetopa); // bc 
  enext2fnext(firsttet, worktet);
  bond(worktet, tetopb); // ca

  // Connect hull tetrahedra together (at six edges of firsttet).
  enext0fnext(tetopc, worktet); 
  enext0fnext(tetopd, worktet1);
  bond(worktet, worktet1); // ab
  enext0fnext(tetopa, worktet);
  enext2fnext(tetopd, worktet1);
  bond(worktet, worktet1); // bc
  enext0fnext(tetopb, worktet);
  enextfnext(tetopd, worktet1);
  bond(worktet, worktet1); // ca
  enext2fnext(tetopc, worktet);
  enextfnext(tetopb, worktet1);
  bond(worktet, worktet1); // da
  enext2fnext(tetopa, worktet);
  enextfnext(tetopc, worktet1);
  bond(worktet, worktet1); // db
  enext2fnext(tetopb, worktet);
  enextfnext(tetopa, worktet1);
  bond(worktet, worktet1); // dc

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
// the Delaunayness of T. Otherwise, do nothing with regard to the Delaunay- //
// ness of T (T may be non-Delaunay after this function).                    //
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
  long tetcount;
  bool enqflag;
  int copcount;
  int *iptr, i, j;

  clock_t loc_start, loc_end;

  if (b->verbose > 1) {
    printf("    Insert point %d\n", pointmark(insertpt));
  }

  // loc_start = clock();

  tetcount = ptloc_count;
  
  if (searchtet->tet == NULL) {
    randomsample(insertpt, searchtet);
  }
  loc = locate(insertpt, searchtet);

  // loc_end = clock();
  // tloctime += ((REAL) (loc_end - loc_start)) / CLOCKS_PER_SEC;

  if (b->verbose > 1) {
    printf("    Walk distance (# tets): %ld\n", ptloc_count - tetcount);
  }

  if (ptloc_max_count < (ptloc_count - tetcount)) {
    ptloc_max_count = (ptloc_count - tetcount);
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

  // The number of coplanar cavity boundary face.
  copcount = 0;

  // loc_start = clock();

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

  // Form a star-shaped cavity with respect to the inserting point.
  for (i = 0; i < cavetetlist->len(); i++) {
    cavetet = (triface *) cavetetlist->get(i);
    for (cavetet->loc = 0; cavetet->loc < 4; cavetet->loc++) {
      sym(*cavetet, neightet);
      // Do check if it is not deleted and not infected.
      if ((neightet.tet[4] != NULL) && !infected(neightet)) {
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
          //   This happens when p lies outside the hull face.
          ori = orient3d(pts[4], pts[5], pts[6], insertpt); orient3dcount++;
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
          // Found a boundary face of the cavity. It may be a face of a hull
          //   tet which contains 'dummypoint'. Choose the edge in the face 
          //   such that its endpoints are not 'dummypoint', while its apex
          //   may be 'dummypoint' (see Fig. 1.4).
          neightet.ver = 4;
          cavebdrylist->append(&neightet);
        }
      } // if (neightet.tet[4] != NULL)
    } // for (cavetet->loc = 0;
    // Delete the tet in the cavity.
    tetrahedrondealloc(cavetet->tet);
  } // for (i = 0;

  if (b->verbose > 1) {
    printf("    Size of the cavity: %d faces %d tets.\n", cavebdrylist->len(),
           cavetetlist->len());
  }
  
  // Create new tetrahedra in the Bowyer-Watson cavity. Connect them to the
  //   tetrahedra at outside of the cavity.
  for (i = 0; i < cavebdrylist->len(); i++) {
    cavetet = (triface *) cavebdrylist->get(i);
    neightet = *cavetet;
    if (apex(neightet) != dummypoint) {
      // Create a new tet in the cavity (see Fig. bowyerwatson 1 or 3).
      unmarktest(neightet); // Unmark it.
      maketetrahedron(tetrahedronpool, &newtet);
      setorg(newtet, dest(neightet));
      setdest(newtet, org(neightet));
      setapex(newtet, apex(neightet));
      setoppo(newtet, insertpt);
    } else {
      // Create a new hull tet (see Fig. bowyerwatson 2).
      maketetrahedron(hulltetrahedronpool, &newtet);
      setorg(newtet, org(neightet));
      setdest(newtet, dest(neightet));
      setapex(newtet, insertpt);
      setoppo(newtet, dummypoint);
      // Note: the cavity boundary face is at the enext0fnext place.
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

  // loc_end = clock();
  // tinserttime += ((REAL) (loc_end - loc_start)) / CLOCKS_PER_SEC;

  // Set a handle for the point location.
  recenttet = * (triface *) cavebdrylist->get(0);
  
  if (bowyerwatson && (copcount > 0)) {
    // There may exist degenerate tetrahedra. Check and remove them.
    bowyerwatsonpostproc(cavebdrylist);
  }

  // loc_start = clock();

  // If the flip option is used.
  if (!bowyerwatson && increflip) {
    lawsonflip(cavebdrylist);
  }
  
  // loc_end = clock();
  // tfliptime += ((REAL) (loc_end - loc_start)) / CLOCKS_PER_SEC;

  // Set the point type.
  if (pointtype(insertpt) == UNUSEDVERTEX) {
    pointtype(insertpt) = VOLVERTEX;
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
  triface fliptets[4];
  triface *cavetet, neightet;
  point *pts;
  REAL ori;
  int i, j;
  
  tetrahedron ptr;
  int *iptr;

  removeque = new queue(sizeof(triface));
  
  // Queue all degenerate tets.
  for (i = 0; i < cavebdrylist->len(); i++) {
    cavetet = (triface *) cavebdrylist->get(i);
    pts = (point *) cavetet->tet;
    if (pts[7] != dummypoint) {
      // Check if the new tet is degenerate.
      ori = orient3d(pts[4], pts[5], pts[6], pts[7]); orient3dcount++;
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
      fliptets[0] = *cavetet;
      for (j = 0; j < 3; j++) {
        fnext(fliptets[j], fliptets[j + 1]);
      }
      if (fliptets[3].tet != cavetet->tet) {
        printf("Internal error in insertvertex(): Unknown flip case.\n");
        terminatetetgen(1);
      }
      // Do a 3-to-2 flip to remove the degenerate tet.
      flip32(fliptets, 0);
      // Rememebr the new tet.
      recenttet = fliptets[0];
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
  triface fliptets[5], baktets[2];
  triface fliptet, neightet;
  point *pts, pd, pe;
  REAL sign, ori;
  long flipcount;
  int n, i;

  int *iptr;

  // Put all boundary faces (except hull faces) into a stack. 
  futureflip = (badface *) NULL;
  for (i = 0; i < cavebdrylist->len(); i++) {
    fliptet = * (triface *) cavebdrylist->get(i);
    sym(fliptet, neightet);
    if ((point) neightet.tet[7] != dummypoint) {
      futureflip = flippush(futureflip, &fliptet, oppo(fliptet));
    }
  }

  if (b->verbose > 1) {
    printf("    Lawson flip %ld faces.\n", flippool->items);
    flipcount = flip23count + flip32count;
  }

  while (futureflip != (badface *) NULL) {

    // Pop a face from the stack.
    fliptet = futureflip->tt;  // bace
    pd = futureflip->foppo;  // The new vertex.
    futureflip = futureflip->nextitem;

    pts = (point *) fliptet.tet;
    // Skip it if it is a dead or a hull tet.
    if ((pts[4] == NULL) || (pts[7] == dummypoint)) continue;
    // Skip it if it is not the same tet as we saved.
    if (oppo(fliptet) != pd) continue;
    // Skip it if the opposite tet does not exist.
    sym(fliptet, neightet);
    pe = oppo(neightet);
    if (pe == dummypoint) continue;

    sign = insphere_sos(pts[4], pts[5], pts[6], pts[7], pe);

    if (b->verbose > 2) {
      printf("  Insphere: (%d, %d, %d) %d, %d\n", pointmark(org(fliptet)),
        pointmark(dest(fliptet)), pointmark(apex(fliptet)),
        pointmark(oppo(fliptet)), pointmark(pe));
    }

    // Flip it if it is not locally Delaunay.
    if (sign < 0) {
      // Check the convexity of its three edges.
      fliptet.ver = 0;
      for (i = 0; i < 3; i++) {
        ori = orient3d(org(fliptet), dest(fliptet), pd, pe); orient3dcount++;
        if (ori <= 0) break;
        enextself(fliptet);
      }
      if (i == 3) {
        // A 2-to-3 flip is found.
        fliptets[0] = fliptet; // tet abcd, d is the new vertex.
        symedge(fliptets[0], fliptets[1]); // tet bace.
        flip23(fliptets, 1);
        recenttet = fliptets[0]; // for point location.
      } else {
        // A 3-to-2 or 4-to-4 may possible.
        enext0fnext(fliptet, fliptets[0]);
        esymself(fliptets[0]); // tet badc, d is the new vertex.
        n = 0;
        do {
          fnext(fliptets[n], fliptets[n + 1]);
          n++;
        } while ((fliptets[n].tet != fliptet.tet) && (n < 5));
        if (n == 3) {
          // Found a 3-to-2 flip.
          flip32(fliptets, 1);
          recenttet = fliptets[0]; // for point location.
        } else if ((n == 4) && (ori == 0)) {
          // Find a 4-to-4 flip.
          flipnmcount++;
          // First do a 2-to-3 flip.
          fliptets[0] = fliptet; // tet abcd, d is the new vertex.
          baktets[0] = fliptets[2];
          baktets[1] = fliptets[3];
          flip23(fliptets, 1);
          // Then do a 3-to-2 flip. 
          enextfnextself(fliptets[0]);  // fliptets[0] is edab.
          enextself(fliptets[0]);
          esymself(fliptets[0]);  // tet badc, d is the new vertex.
          fliptets[1] = baktets[0];
          fliptets[2] = baktets[1];
          flip32(fliptets, 1);
          recenttet = fliptets[0]; // for point location.
        } else {
          // An unflipable face. Ignore it. Will be flipped later.
        }
      } // if (i == 3)
    } // if (sign < 0)
  }

  if (b->verbose > 1) {
    printf("    %ld flips.\n", flip23count + flip32count - flipcount);
  }

  flippool->restart();
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