#ifndef delaunayCXX
#define delaunayCXX

#include "tetgen.h"

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

  // Remember the first tetrahedron.
  recenttet = firsttet;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// bowyerwatsoninsert()    Insert a point by Bowyer-Watson algorithm.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::bowyerwatsoninsert(point insertpt, triface* firsttet)
{
  list *cavetetlist, *cavebdrylist;
  triface *cavetet, neightet, newtet;
  point *pts;
  REAL sign, ori;
  bool enqflag;
  int *iptr;
  int i, j;

  // Initialize working lists.
  cavetetlist = new list(sizeof(triface));
  cavebdrylist = new list(sizeof(triface));

  // Collect all non-Delaunay tetrahedra. Form the Bowyer-Watson cavity.
  //   Starting from 'firsttet', do breath-first search around in its
  //   adjacent tetrahedra, and the adjacent of adjacent tetrahedra ...
  infect(*firsttet);
  cavetetlist->append(firsttet);
  for (i = 0; i < cavetetlist->len(); i++) {
    cavetet = (triface *) cavetetlist->get(i);
    for (cavetet->loc = 0; cavetet->loc < 4; cavetet->loc++) {
      sym(*cavetet, neightet);
      // Do check if it is not included in 'cavetetlist' yet.
      if (!infected(neightet)) {
        enqflag = false;
        pts = (point *) neightet.tet;
        if (pts[7] != dummypoint) {
          // A voolume tet. Do Delaunay check if it has not been tested yet. 
          if (!marktested(neightet)) {
            sign = insphere(pts[4], pts[5], pts[6], pts[7], insertpt);
            marktest(neightet); // Only test it once.
            enqflag = (sign < 0.0);
          }
        } else {
          // It is a hull tet. Check the following two cases:
          //   (1) its base face is visible by 'insertpt'. This happens
          //       when 'insertpt' lies outside the hull face;
          //   (2) its base face is coplanar with 'insertpt'.  This happens
          //       when 'insertpt' lies just on a hull face or a hull edge.
          // In both cases, the cuurent convex hull need to be updated.
          ori = orient3d(pts[4], pts[5], pts[6], insertpt);
          enqflag = (ori <= 0.0);
        } 
        if (enqflag) {
          // A non-Delaunay tet or a dead hull tet.
          infect(neightet);
          cavetetlist->append(&neightet);
        } else {
          // Found a boubdary face of the cavity.
          cavebdrylist->append(cavetet);
        }
      } // if (neightet.tet != NULL)
    } // for (cavetet->loc = 0;
  } // for (i = 0;
  
  // Create new tetrahedra in the Bowyer-Watson cavity. Connect them to the
  //   tetrahedra at outside of the cavity.
  for (i = 0; i < cavebdrylist->len(); i++) {
    cavetet = (triface *) cavebdrylist->get(i);
    sym(*cavetet, neightet); // neightet may be a hull tet.
    pts = (point *) cavetet->tet;
    if (pts[7] != dummypoint) {
      // Create a new tet in the cavity.
      unmarktest(neightet); // Unmark it.
      neightet.ver = 0; // Orient the face. 
      maketetrahedron(tetrahedronpool, &newtet);
      setorg(newtet, dest(neightet));
      setdest(newtet, org(neightet));
      setapex(newtet, apex(neightet));
      setoppo(newtet, insertpt);
      bond(newtet, neightet);
    } else {
      // Create a new hull tet (insertpt is on the hull).
      cavetet->ver = 4; // The edge's apex is 'dummypoint' (see Fig. 1.4).
      assert(apex(*cavetet) == dummypoint); // SELF_CHECK
      // Create a new hull tetrahedron.
      maketetrahedron(hulltetrahedronpool, &newtet);
      setorg(newtet, dest(*cavetet));  // Refer to Fig. 1.4.
      setdest(newtet, org(*cavetet));
      setapex(newtet, insertpt);
      setoppo(newtet, dummypoint);
      // Note: the hull face is at the enext0fnext place.
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
      // Go to the face needs to be connected.
      enext0fnext(*cavetet, newtet);
      // Do connection if it has not yet been connected.
      if (newtet.tet[newtet.loc] == NULL) {
        // Find its adjacent face by rotating faces around the edge of
        //   cavetet. The rotating direction is opposite to newtet.
        //   Stop the rotate at a face which has no adjacent tet.
        esym(*cavetet, neightet); // Set the rotate dir.
        do {
          // Go to the face in the adjacent tet.
          fnextself(neightet);
          // Continue if the adjacent tet exists.
        } while (neightet.tet[neightet.loc] != NULL);
        assert(apex(neightet) == dummypoint); // SELF_CHECK.
        // Connect newtet <==> neightet.
        bond(newtet, neightet);
      }
      enextself(*cavetet);
    }
  }

  // Delete the non-Delaunay tetrahedra.
  for (i = 0; i < cavetetlist->len(); i++) {
    cavetet = (triface *) cavetetlist->get(i);
    pts = (point *) cavetet->tet;
    if (pts[7] != dummypoint) {
      tetrahedrondealloc(tetrahedronpool, cavetet->tet);
    } else {
      tetrahedrondealloc(hulltetrahedronpool, cavetet->tet);
    }
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
  triface searchtet, neightet;
  point *permutarray, swapvertex;
  enum locateresult ptloc;
  REAL v1[3], v2[3], n[3];
  REAL bboxsize, bboxsize2, bboxsize3, ori;
  int randindex, i, j;

  // Form a random permuation (uniformly at random) of the set of vertices.
  permutarray = new point[in->numberofpoints];
  pointpool->traversalinit();
  for (i = 0; i < in->numberofpoints; i++) {
    randindex = randomnation(i + 1);
    permutarray[i] = permutarray[randindex];
    permutarray[randindex] = (point) pointpool->traverse();
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
      v1[j] = permutarray[1][j] - permutarray[0][j];
      v2[j] = permutarray[i][j] - permutarray[i][j];
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
  while ( (fabs(ori) / bboxsize3) < b->epsilon) {
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
    permutarray[i] = permutarray[2];
    permutarray[2] = swapvertex;
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
    // Locate the point in DT.
    searchtet.tet = NULL;
    ptloc = locate(permutarray[i], &searchtet);
    // Insert the point by Bowyer-Watson algorithm.
    bowyerwatsoninsert(permutarray[i], &searchtet);    
  } // for (i = 4; i < in->numberofpoints; i++)

  delete [] permutarray;
}

#endif // #ifndef delaunayCXX