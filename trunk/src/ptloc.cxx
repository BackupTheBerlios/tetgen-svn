#ifndef ptlocCXX
#define ptlocCXX

#include "tetgen.h"

#define NORM2(x, y, z) (x) * (x) + (y) * (y) + (z) * (z)

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
// locate()    Find a simplex containing a given point.                      //
//                                                                           //
// This routine implements Muecke's Jump-and-walk point location algorithm.  //
// It improves the simple walk-through by "jumping" to a good starting point //
// via random sampling.  Searching begins from one of handles:  the input    //
// 'searchtet', a recently encountered tetrahedron 'recenttet',  or from one //
// chosen from a random sample.  The choice is made by determining which one //
// 's origin is closest to the point we are searcing for.  Having chosen the //
// starting tetrahedron, the simple Walk-through algorithm is executed.      //
//                                                                           //
// The return value indicates the location of the 'searchpt'. 'searchtet' is //
// adjusted to a tetrahedron corresponding to that value.                    //
//                                                                           //
// WARNING: This routine is designed for convex triangulations, and will not //
// generally work after the holes and concavities have been carved.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::locateresult tetgenmesh::locate(point searchpt, 
  triface *searchtet)
{
  tetrahedron *firsttet, *tetptr;
  point torg;
  void **sampleblock;
  long sampleblocks, samplesperblock, samplenum;
  unsigned long alignptr;
  REAL searchdist, dist;
  int tetblocks, i, j;

  if (searchtet->tet != NULL) {
    // Get the distance from the suggested starting tet to the search point.
    torg = org(*searchtet);
    searchdist = NORM2(searchpt[0] - torg[0], searchpt[1] - torg[1], 
                       searchpt[2] - torg[2]);
  } else {
    searchdist = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
  }

  // If a recently encountered tetrahedron has been recorded and has not
  //   been deallocated, test it as a good starting point.
  if ((recenttet.tet != NULL) && (recenttet.tet != searchtet->tet)) {
    torg = org(recenttet);
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
  samplesperblock = 1 + (samples / tetblocks);
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
  
  // Call simple walk-through to locate the point.
  return preciselocate(searchpt, searchtet); 
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// preciselocate()    Find a simplex containing a given point.               //
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

enum tetgenmesh::locateresult tetgenmesh::preciselocate(point searchpt,
  triface* searchtet)
{
  tetrahedron ptr;
  triface neightet;
  point torg, tdest, tapex, toppo, ntoppo;
  enum {ORGMOVE, DESTMOVE, APEXMOVE} nextmove;
  REAL ori, oriorg, oridest, oriapex;
  REAL searchdist, dist;
  int *iptr;

  // Let searchtet point to its face such that 'searchpt' lies above to it.
  searchtet->ver &= 01; // Keep the CCW edge ring.
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
    if (oriorg > 0) {
      if (oridest > 0) {
        if (oriapex > 0) {
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
        if (oriapex > 0) {
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
      if (oridest > 0) {
        if (oriapex > 0) {
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
        if (oriapex > 0) {
          // Only the face opposite to apex is viable.
          nextmove = APEXMOVE;
        } else {
          // The point we seek must be on the boundary of or inside this
          //   tetrahedron. Check for boundary cases.
          if (oriorg == 0) {
            // Go to the face opposite to origin.
            enextfnextself(*searchtet);
            if (oridest == 0) {
              enextself(*searchtet);
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

#endif // #ifndef ptlocCXX