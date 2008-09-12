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
  setorg(firsttet, pa);
  setdest(firsttet, pb);
  setapex(firsttet, pc);
  setoppo(firsttet, pd);
  // Create four hull tetrahedra.
  maketetrahedron(hulltetrahedronpool, &tetopa);
  setorg(tetopa, pb);
  setdest(tetopa, pc);
  setapex(tetopa, pd);
  setoppo(tetopa, dummypoint);
  maketetrahedron(hulltetrahedronpool, &tetopb);
  setorg(tetopb, pc);
  setdest(tetopb, pa);
  setapex(tetopb, pd);
  setoppo(tetopb, dummypoint);
  maketetrahedron(hulltetrahedronpool, &tetopc);
  setorg(tetopc, pa);
  setdest(tetopc, pb);
  setapex(tetopc, pd);
  setoppo(tetopc, dummypoint);
  maketetrahedron(hulltetrahedronpool, &tetopd);
  setorg(tetopd, pb);
  setdest(tetopd, pa);
  setapex(tetopd, pc);
  setoppo(tetopd, dummypoint);

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
  triface *cavetet, neightet, newtet;
  point *pts;
  enum locateresult loc;
  REAL sign, ori;
  bool enqflag;
  int *iptr;
  int i, j;

  if (b->verbose > 1) {
    printf("    Insert point %d\n", pointmark(insertpt));
  }

  // Locate the point p.
  loc = locate(insertpt, searchtet);
  
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
  
  // Connect the set of new tetrahedra of the cavity together.
  for (i = 0; i < cavebdrylist->len(); i++) {
    cavetet = (triface *) cavebdrylist->get(i);
    cavetet->ver = 0;
    for (j = 0; j < 3; j++) {
      // Go to the face needs to be connected.
      enext0fnext(*cavetet, newtet);
      // Operate on it if it has not yet been connected.
      if (newtet.tet[newtet.loc] == NULL) {
        // Find its adjacent face by rotating faces around the edge of
        //   cavetet. The rotating direction is opposite to newtet.
        //   Stop the rotate at a face which has no adjacent tet.
        esym(*cavetet, neightet); // Set the rotate dir.
        do {
          fnextself(neightet); // Go to the face in the adjacent tet.
        } while (neightet.tet[neightet.loc] != NULL);
        // Connect newtet <==> neightet.
        bond(newtet, neightet);
      }
      enextself(*cavetet);
    }
  }

  // Delete the old tetrahedra of the cavity.
  for (i = 0; i < cavetetlist->len(); i++) {
    cavetet = (triface *) cavetetlist->get(i);
    if ((point) cavetet->tet[7] != dummypoint) {
      tetrahedrondealloc(tetrahedronpool, cavetet->tet);
    } else {
      tetrahedrondealloc(hulltetrahedronpool, cavetet->tet);
    }
  }

  // If the flip option is used.
  if (!bowyerwatson && increflip) {
    // flip(cavebdrylist);
  }

  // Set the point type.
  if (pointtype(insertpt) == UNUSEDVERTEX) {
    pointtype(insertpt) = VOLVERTEX;
  }

  // Set a handle for the point location.
  recenttet = * (triface *) cavebdrylist->get(0);
  
  delete cavetetlist;
  delete cavebdrylist;
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
    // Skip the permutation.
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
    // Insert the point by Bowyer-Watson algorithm.
    searchtet.tet = NULL;
    insertvertex(permutarray[i], &searchtet, true, false);
  }

  delete [] permutarray;
}

#endif // #ifndef delaunayCXX