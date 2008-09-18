#ifndef flipCXX
#define flipCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flippush()    Push a face (possibly will be flipped) into stack.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::badface* tetgenmesh::flippush(badface* flipstack, 
  triface* flipface, point pushpt)
{
  badface *newflipface;

  newflipface = (badface *) flippool->alloc();
  newflipface->tt = *flipface;
  newflipface->foppo = pushpt;
  newflipface->nextitem = flipstack;

  return newflipface;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip23()    Remove a face by tranforming 2-to-3 tetrahedra.               //
//                                                                           //
// This is the reverse operation of flip23(). 'oldtets' are two origin tets, //
// abcd and bace, 'newtets' returns three new tets, edab, edbc, and edca.    //
//                                                                           //
// In case there are hull tets involved in this flip.  There are two cases:  //
//   (1) If d is 'dummypoint', all three new tets are hull tets.  If e is    //
//       'dummypoint', we reconfigure e to d, i.e., turn it up-side down.    //
//   (2) If c is 'dummypoint' , two new tets edbc and edca are hull tets.    //
//       If a or b is 'dummypoint', we reconfigure it to c, i.e., rotate the //
//       three old tets counterclockwisely (right-hand rule) until a or b    //
//       is in c's position (see Fig.).                                      //
//                                                                           //
// If 'flipflag != 0', the convex hull faces will be checked and flipped if  //
// they are locally non-Delaunay. If 'flipflag == 1', we assume that d is a  //
// newly inserted vertex, so only the lower part of the convex hull will be  //
// cehcked (this avoids unnecessary insphere() testes).                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip23(triface* oldtets, triface* newtets, int flipflag)
{
  tetrahedron ptr;
  triface newface, casface;
  point pa, pb, pc, pd, pe;
  int dummyflag;  // range = {-1, 0, 1, 2}.
  int *iptr, i;

  // Check if e is dummypoint.
  if (oppo(oldtets[1]) == dummypoint) {
    // Swap the two old tets.
    newface = oldtets[0];
    oldtets[0] = oldtets[1];
    oldtets[1] = newface;
    dummyflag = -1;  // e is dummypoint.
  } else {
    // Check if a or b is dummypoint.
    if (org(oldtets[0]) == dummypoint) {
      dummyflag = 1;  // a is dummypoint.
    } else if (dest(oldtets[0]) == dummypoint) {
      dummyflag = 2;  // b is dummypoint.
    } else {
      dummyflag = 0;  // either c or d is dummypoint.
    }
    i = dummyflag;
    // Rotate i times.
    for (; i > 0; i--) {
      enextself(oldtets[0]);
      enext2self(oldtets[1]);
    }
  }

  pa = org(oldtets[0]);
  pb = dest(oldtets[0]);
  pc = apex(oldtets[0]);
  pd = oppo(oldtets[0]);
  pe = oppo(oldtets[1]);

  if (b->verbose > 1) {
    printf("    flip 2-to-3: (%d, %d, %d, %d, %d)\n", pointmark(pa),
      pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pe));
  }
  flip23count++;

  // Check if d is dummytet.
  if (pd != dummypoint) {
    // Create edab.
    maketetrahedron(tetrahedronpool, &(newtets[0]));
    setvertices(newtets[0], pe, pd, pa, pb);
    // Check if c is dummypoint.
    if (pc != dummypoint) {
      // Create edbc
      maketetrahedron(tetrahedronpool, &(newtets[1]));
      // Create edca
      maketetrahedron(tetrahedronpool, &(newtets[2]));
      setvertices(newtets[2], pe, pd, pc, pa);
    } else {
      // Create edbc
      maketetrahedron(hulltetrahedronpool, &(newtets[1]));
      // Create deac
      maketetrahedron(hulltetrahedronpool, &(newtets[2]));
      setvertices(newtets[2], pd, pe, pa, pc);
      // Adjust deac->edca
      enext0fnextself(newtets[2]);
      esymself(newtets[2]);
    }
    setvertices(newtets[1], pe, pd, pb, pc);
  } else {
    // d is dummypoint. Create three new hull tets.
    // Create abed
    maketetrahedron(hulltetrahedronpool, &(newtets[0]));
    setvertices(newtets[0], pa, pb, pe, pd);
    // Create bced
    maketetrahedron(hulltetrahedronpool, &(newtets[1]));
    setvertices(newtets[1], pb, pc, pe, pd);
    // Create caed
    maketetrahedron(hulltetrahedronpool, &(newtets[2]));
    setvertices(newtets[2], pc, pa, pe, pd);
    // Adjust abed->edab, bced->edbc, caed->edca
    for (i = 0; i < 3; i++) {
      enext2fnextself(newtets[i]);
      enext2self(newtets[i]);
      esymself(newtets[i]);
    }
  }

  // Bond three new tets together.
  for (i = 0; i < 3; i++) {
    enext0fnext(newtets[i], newface);
    bond(newface, newtets[(i + 1) % 3]);
  }
  // Bond top boundary faces (at abcd) to mesh.
  for (i = 0; i < 3; i++) {
    enextfnext(newtets[i], newface);
    enext0fnext(oldtets[0], casface);
    symself(casface);
    bond(newface, casface);
    enextself(oldtets[0]);
  }
  // Bond bottom boundary faces (at bace) to mesh.
  for (i = 0; i < 3; i++) {
    enext2fnext(newtets[i], newface);
    enext0fnext(oldtets[1], casface);
    symself(casface);
    bond(newface, casface);
    enext2self(oldtets[1]);
  }
  
  if (dummyflag != 0) {
    // Restore the original position of the points (for flipnm()).
    if (dummyflag == -1) { 
      // d and e were swapped.
      newface = oldtets[0];
      oldtets[0] = oldtets[1];
      oldtets[1] = newface;
      // Reverse the edge.
      for (i = 0; i < 3; i++) {
        enext0fnextself(newtets[i]);
        esymself(newtets[i]);
      }
      // Swap the last two new tets.
      newface = newtets[1];
      newtets[1] = newtets[2];
      newtets[2] = newface;
    } else {
      // either a or b were swapped.
      i = dummyflag;
      // Rotate i times.
      for (; i > 0; i--) {
        enext2self(oldtets[0]);
        enextself(oldtets[1]);
      }
      i = dummyflag;
      // Down-shift new tets i times.
      for (; i > 0; i--) {
        newface = newtets[0];
        newtets[0] = newtets[2];
        newtets[2] = newtets[1];
        newtets[1] = newface;
      }
    }
  }

  if (flipflag > 0) {
    // Queue faces which may be locally non-Delaunay.  
    pd = dest(newtets[0]);
    for (i = 0; i < 3; i++) {
      enext2fnext(newtets[i], newface);
      sym(newface, casface);
      if ((point) casface.tet[7] != dummypoint) {
        futureflip = flippush(futureflip, &casface, pd);
      }
    }
    if (flipflag > 1) {
      pe = org(newtets[0]);
      for (i = 0; i < 3; i++) {
        enextfnext(newtets[i], newface);
        sym(newface, casface);
        if ((point) casface.tet[7] != dummypoint) {
          futureflip = flippush(futureflip, &casface, pe);
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip32()    Remove an edge by transforming 3-to-2 tetrahedra.             //
//                                                                           //
// The three old tetrahedra, denoted edab, edbc, and edca, share at the edge //
// de, are given in 'oldtets'.  The flip32() operation replace them into two //
// new tetrahedra, abcd and bace, in 'newtets'.  At a result, the edge de is //
// replaced by the face abc within the convex hull of the points.            //
//                                                                           //
// In case there are hull tets involved in this flip.  There are two cases:  //
//   (1) If d is 'dummypoint', then abcd is hull tet, and bace is normal.    //
//       If e is 'dummypoint', we reconfigure e to d, i.e., turnover it.     //
//   (2) If c is 'dummypoint' then both abcd and bace are hull tets.         //
//       If a or b is 'dummypoint', we reconfigure it to c, i.e., rotate the //
//       three old tets counterclockwisely (right-hand rule) until a or b    //
//       is in c's position (see Fig.).                                      //
//                                                                           //
// If 'flipflag != 0', the convex hull faces will be checked and flipped if  //
// they are locally non-Delaunay. If 'flipflag == 1', we assume that a is a  //
// newly inserted vertex, so only the two opposite faces of the convex hull  //
// will be cehcked (this avoids unnecessary insphere() testes).              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip32(triface* oldtets, triface* newtets, int flipflag)
{
  tetrahedron ptr;
  triface newface, casface;
  point pa, pb, pc, pd, pe;
  int dummyflag;  // Rangle = {-1, 0, 1, 2}
  int *iptr, i;

  // Check if e is 'dummypoint'.
  if (org(oldtets[0]) == dummypoint) {
    // Reverse the edge.
    for (i = 0; i < 3; i++) {
      enext0fnextself(oldtets[i]);
      esymself(oldtets[i]);
    }
    // Swap the last two tets.
    newface = oldtets[1];
    oldtets[1] = oldtets[2];
    oldtets[2] = newface;
    dummyflag = -1; // e is dummypoint.
  } else {
    // Check if a or b is the 'dummypoint'.
    if (apex(oldtets[0]) == dummypoint) { 
      dummyflag = 1;  // a is dummypoint.
    } else if (apex(oldtets[1]) == dummypoint) {
      dummyflag = 2;  // b is dummypoint.
    } else {
      dummyflag = 0;  // either c or d is dummypoint.
    }
    i = dummyflag;
    // Down-shift i times.
    for (; i > 0; i--) {
      newface = oldtets[2];
      oldtets[2] = oldtets[1];
      oldtets[1] = oldtets[0];
      oldtets[0] = newface;
    }
  }

  pa = apex(oldtets[0]);
  pb = apex(oldtets[1]);
  pc = apex(oldtets[2]);
  pd = dest(oldtets[0]);
  pe = org(oldtets[0]);

  if (b->verbose > 1) {
    printf("    flip 3-to-2: (%d, %d, %d, %d, %d)\n", pointmark(pa),
      pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pe));
  }
  flip32count++;

  // Check if c is a dummypointc.
  if (pc != dummypoint) {
    // Create abcd. Check if d is a dummypoint.
    if (pd != dummypoint) {
      maketetrahedron(tetrahedronpool, &(newtets[0]));
    } else {
      maketetrahedron(hulltetrahedronpool, &(newtets[0])); // a hull tet.
    }
    setvertices(newtets[0], pa, pb, pc, pd);
    // Create bace.
    maketetrahedron(tetrahedronpool, &(newtets[1]));
    setvertices(newtets[1], pb, pa, pc, pe);
  } else {
    // c is dummypoint. The two new tets are hull tets.
    // Create badc.
    maketetrahedron(hulltetrahedronpool, &(newtets[0]));
    setvertices(newtets[0], pb, pa, pd, pc);
    // Create abec
    maketetrahedron(hulltetrahedronpool, &(newtets[1]));
    setvertices(newtets[1], pa, pb, pe, pc);
    // Adjust badc -> abcd.
    enext0fnextself(newtets[0]);
    esymself(newtets[0]);
    // Adjust abec -> bace.
    enext0fnextself(newtets[1]);
    esymself(newtets[1]);
  }
  
  // Bond abcd <==> bace.
  bond(newtets[0], newtets[1]);  
  // Bond other faces of abcd to mesh.
  for (i = 0; i < 3; i++) {
    enext0fnext(newtets[0], newface);
    enextfnext(oldtets[i], casface);
    symself(casface);
    bond(newface, casface);
    enextself(newtets[0]);
  }
  // Bond other faces of bace to mesh.
  for (i = 0; i < 3; i++) {
    enext0fnext(newtets[1], newface);
    enext2fnext(oldtets[i], casface);
    symself(casface);
    bond(newface, casface);
    enext2self(newtets[1]);
  }

  if (dummyflag != 0) {
    // Restore the original position of the points (for flipnm()).
    if (dummyflag == -1) {
      // e were dummypoint.
      for (i = 0; i < 3; i++) {
        enext0fnextself(oldtets[i]);
        esymself(oldtets[i]);
      }
      // Swap the last two tets.
      newface = oldtets[1];
      oldtets[1] = oldtets[2];
      oldtets[2] = newface;
      // Swap the two new tets.
      newface = newtets[0];
      newtets[0] = newtets[1];
      newtets[1] = newface;
    } else {
      // a or b was dummypoint.
      i = dummyflag;
      // Down-shift i times.
      for (; i > 0; i--) {
        newface = oldtets[0];
        oldtets[0] = oldtets[2];
        oldtets[2] = oldtets[1];
        oldtets[1] = newface;
      }
      i = dummyflag;
      // Rotate toward left i times.
      for (; i > 0; i--) {
        enextself(newtets[0]);
        enext2self(newtets[1]);
      }
    }
  }
  
  if (flipflag > 0) {
    // Queue faces which may be locally non-Delaunay.  
    pa = org(newtets[0]);
    enextfnext(newtets[0], newface);
    sym(newface, casface);
    if ((point) casface.tet[7] != dummypoint) {
      futureflip = flippush(futureflip, &casface, pa);
    }
    enext2fnext(newtets[1], newface);
    sym(newface, casface);
    if ((point) casface.tet[7] != dummypoint) {
      futureflip = flippush(futureflip, &casface, pa);
    }
    if (flipflag > 1) {
      pb = dest(newtets[0]);
      enext2fnext(newtets[0], newface);
      sym(newface, casface);
      if ((point) casface.tet[7] != dummypoint) {
        futureflip = flippush(futureflip, &casface, pb);
      }
      enextfnext(newtets[1], newface);
      sym(newface, casface);
      if ((point) casface.tet[7] != dummypoint) {
        futureflip = flippush(futureflip, &casface, pb);
      }
      pc = apex(newtets[0]);
      enext0fnext(newtets[0], newface);
      sym(newface, casface);
      if ((point) casface.tet[7] != dummypoint) {
        futureflip = flippush(futureflip, &casface, pc);
      }
      enext0fnext(newtets[1], newface);
      sym(newface, casface);
      if ((point) casface.tet[7] != dummypoint) {
        futureflip = flippush(futureflip, &casface, pc);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipnm()    Remove an edge by transforming n-to-m tetrahedra.             //
//                                                                           //
// This routine attemps to remove an edge, denoted as ab, by transforming a  //
// set K of n tetrahedra containing ab into a set K' of m tetrahedra, where  //
// K and K' have the same outer boundary, and ab is not in K', where n >= 3, //
// m = 2 * n - 4. It can be viwed as a n-to-m flip.                          //
//                                                                           //
// 'oldtets' contains n tets sharing at ab. Imaging that ab perpendicularly  //
// crosses your screen, b lies in front of a. Let the projections of the n   //
// apexes, denoted as p[0], p[1], ..., p[n-1], onto screen in counterclock-  //
// wise order (right-hand rule). The n tets are: abp[0]p[1], abp[1]p[2],..., //
// abp[n-1]p[0], respectively.  If one of p[i] is dummypoint, we reconfigure //
// it such that p[0] is dummypoint. a or b should not be dummypoint.         //
//                                                                           //
// The principle of the transformation is to recursively reduce the link of  //
// ab by perform 2-to-3 flips until the link size of ab is 3, hence a 3-to-2 //
// flip can be applied to remove the edge.                                   //
//                                                                           //
// Consider a face abp[i] (i in [0, n-1]), a flip23() can be applied if the  //
// edge p[(i-1) % n]p[(i+1) % n] crosses it in its interior. Relabel p[i] to //
// c, p[(i-1) % n] to e, and p[(i+1) % n] to d. This transforms the two tets,//
// abcd and bace, to three tets, edab, edbc, and edca. The tet edab remains  //
// in the star of ab, while edbc, and edca do not. As a result, p[i] (c) is  //
// not on the link of ab anymore. The link size is reduced by 1. A recursion //
// can be applied if n - 1 > 3.                                              //
//                                                                           //
// NOTE: In above, the flip23() can be applied even if the edge de crosses   //
// ab, i.e., a, b, e, and d are coplanar. Although this will creare a degen- //
// erate tet edab (has zero volume), it will be removed after ab is removed. //
//                                                                           //
// The number m can be counted as follows: there are n - 3 '2-to-3' flips, 1 //
// '3-to-2' flip. Each '2-to-3' flip produces two new tets, and the '3-to-2' //
// flip produces two new tets, hence m = (n - 3) * 2 + 2 = 2 * n - 4.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*
bool tetgenmesh::flipnm(int n, triface* oldtets, triface* newtets,
  bool delaunay, queue* flipque)
{
  triface *recuroldtets, tmpoldtets[2], baktet;
  point pa, pb, pc, pd, pe, pf;
  bool success, doflip;
  REAL sign, ori;
  int *iptr, i, j;

  // Check if any apex of ab is dummypoint.
  for (i = n - 1; i > 0; i--) {
    if (apex(oldtets[i]) == dummypoint) break;
  }
  // Now i is the shift distance to the first one.
  for (; i > 0; i--) {
    baktet = oldtets[0];
    for (j = 0; j < n - 1; j++) {
      oldtets[j] = oldtets[j + 1];
    }
    oldtets[n - 1] = baktet;
  }

  pa = org(oldtets[0]);
  pb = dest(oldtets[0]);

  if (b->verbose > 1) {
    printf("    flip %d-to-%d: (%d, %d)\n", n, 2 * n - 4, pointmark(pa),
           pointmark(pb));
  }
  flipnmcount++;

  success = false;

  for (i = 0; i < n; i++) {
    pc = apex(oldtets[i]);
    pd = apex(oldtets[(i + 1) % n]);
    pe = apex(oldtets[(i != 0) ? i - 1 : n - 1]);
    // Decide if the face abc is flipable.
    if (pc != dummypoint) {
      if ((pd != dummypoint) && (pc != dummypoint)) {
        ori = orient3d(pa, pb, pd, pe);
        if (ori >= 0) {  // Allow ori == 0, support 4-to-4 flip.
          ori = orient3d(pb, pc, pd, pe);
          if (ori > 0) {
            ori = orient3d(pc, pa, pd, pe);
          }
        }
        doflip = ori > 0;
        if (doflip && delaunay) {
          // Only do flip if abc is not a locally Delaunay face.
          sign = insphere_sos(pa, pb, pc, pd, pe);
          doflip = (sign < 0);
        }
      } else {
        doflip = false;
      }
    } else {
      // Two faces abd and abe are on the hull.
      doflip = true; // 2-to-2 flip is possible.
      if (delaunay) {
        // Only do flip if ab is not a locally Delaunay edge.
        pf = apex(oldtets[(i + 2) % n]);
        sign = insphere_sos(pa, pb, pd, pf, pe);
        doflip = (sign < 0);
      }
    }
    if (doflip) {
      // Get the two old tets.
      tmpoldtets[0] = oldtets[i];
      tmpoldtets[1] = oldtets[(i != 0) ? i - 1 : n - 1];
      // Adjust tmpoldtets[1] (abec) -> bace.
      enext0fnextself(tmpoldtets[1]);
      esymself(tmpoldtets[1]);
      // Adjust the tets so that newtets[2] will be edca (see Fig.).
      enextself(tmpoldtets[0]);
      enext2self(tmpoldtets[1]);
      // Do flip23() on abp[i] (abc).
      flip23(tmpoldtets, newtets, flipque);
      // Form the new star of ab which has n-1 tets.
      recuroldtets = new triface[n - 1];
      // Set the remaining n-2 tets around ab.
      for (j = 0; j < n - 2 ; j++) {
        recuroldtets[j] = oldtets[(i + 1 + j) % n];
      }
      // Put the last tet having ab. Leave a copy in baktet.
      recuroldtets[j] = baktet = newtets[2];
      // Adjust recuroldtets[j] to abp[0]p[1] (see Fig.).
      enext2fnextself(recuroldtets[j]);  // j == n - 2;
      enext2self(recuroldtets[j]);
      esymself(recuroldtets[j]);
      // Check the size of the link of ab.
      if (n > 4) {  // Actually, if (n - 1 > 3)
        // Recursively do flipnm().
        success = flipnm(n-1, recuroldtets, &(newtets[2]), delaunay, flipque);
        // Are we success?
        if (!success) {
          if (!delaunay) {
            // No! Reverse the flip23() operation.
            newtets[2] = baktet; // Do we really need this?
            flip32(newtets, tmpoldtets, flipque);
            // Adjust the tets back to original position (see Fig.).
            enext2self(tmpoldtets[0]);
            enextself(tmpoldtets[1]);
            // Adjust back to abp[(i-1) % n].
            enext0fnextself(tmpoldtets[1]);
            esymself(tmpoldtets[1]);
            // Delete the two original old tets first.
            tetrahedrondealloc(oldtets[i].tet);
            tetrahedrondealloc(oldtets[(i != 0) ? i - 1 : n - 1].tet);
            // Set the tets back to their original positions.
            oldtets[i] = tmpoldtets[0];
            oldtets[(i != 0) ? i - 1 : n - 1] = tmpoldtets[1];
            // Delete the three new tets.
            tetrahedrondealloc(newtets[0].tet);
            tetrahedrondealloc(newtets[1].tet);
            tetrahedrondealloc(newtets[2].tet); // = recuroldtets[j].tet
          } else {
            // Only delete the first two old tets. Their common face is
            //   not locally Delaunay and has been flipped.
            tetrahedrondealloc(oldtets[i].tet);
            tetrahedrondealloc(oldtets[(i != 0) ? i - 1 : n - 1].tet);
            // Here the tet recuroldtets[j] does not get deleted. It is
            //   a part of current mesh.
          }
        } else {
          // Delete the last tet in 'recuroldtets'. This tet is neither in
          //   'oldtets' nor in 'newtets'.
          tetrahedrondealloc(recuroldtets[j].tet);  // j == n - 2
        } 
      } else {
        // Remove ab by a flip32().
        flip32(recuroldtets, &(newtets[2]), flipque);
        // Delete the last tet in 'recuroldtets'. This one is just created
        //   by the flip 2-to-3 operation within this routine.
        tetrahedrondealloc(recuroldtets[j].tet);  // j == n - 2
        success = true;
      }      
      // The other tets in 'recuroldtets' are still in 'oldtets'.
      delete [] recuroldtets;
      break;
    } // if (doflip)
  } // for (i = 0; i < n; i++)
  
  return success;
}
*/

#endif // #ifndef flipCXX