#ifndef flipCXX
#define flipCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip32()    Perform a 3-to-2 tetrahedra transformation.                   //
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
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip32(triface* oldtets, triface* newtets, queue* flipque)
{
  tetrahedron ptr;
  triface newface, casface;
  point pa, pb, pc, pd, pe;
  int *iptr, i, j;

  // Check if e is 'dummypoint'.
  if (org(oldtets[0]) == dummypoint) {  // pe
    // Reconfigure the old tets such that d is 'dummypoint'.
    for (i = 0; i < 3; i++) {
      // edab->deba, edbc->decb, edca->deac.
      enext0fnextself(oldtets[i]);
      esymself(oldtets[i]);
    }
    // Swap the last two tets.
    newface = oldtets[1];
    oldtets[1] = oldtets[2];
    oldtets[2] = newface;
  } else {
    // Check if a or b is the 'dummypoint'.
    j = 0; // No shift.
    if (apex(oldtets[0]) == dummypoint) { 
      j = 2;  // Shift pa->pb->pc.
    } else if (apex(oldtets[1]) == dummypoint) {
      j = 1;  // Shift pb->pc.
    }
    for (; j > 0; j--) {
      newface = oldtets[2];
      for (i = 2; i > 0; i--) {
        oldtets[i] = oldtets[i - 1];
      }
      oldtets[0] = newface;
    }
  }

  pa = apex(oldtets[0]);
  pb = apex(oldtets[1]);
  pc = apex(oldtets[2]);
  pd = dest(oldtets[0]);
  pe = org(oldtets[0]);

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
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip23()    Perform a 2-to-3 tetrahedra transformation.                   //
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
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip23(triface* oldtets, triface* newtets, queue* flipque)
{
  tetrahedron ptr;
  triface newface, casface;
  point pa, pb, pc, pd, pe;
  int *iptr, i, j;

  // Check if e is dummypoint.
  if (oppo(oldtets[1]) == dummypoint) {
    // Swap the two old tets.
    newface = oldtets[0];
    oldtets[0] = oldtets[1];
    oldtets[1] = newface;
  } else {
    // Check if a or b is dummypoint.
    j = 0;
    if (org(oldtets[0]) == dummypoint) {
      j = 2;
    } else if (dest(oldtets[0]) == dummypoint) {
      j = 1;
    }
    for (; j > 0; j--) {
      enextself(oldtets[0]);
      enext2self(oldtets[1]);
    }
  }

  pa = org(oldtets[0]);
  pb = dest(oldtets[0]);
  pc = apex(oldtets[0]);
  pd = oppo(oldtets[0]);
  pe = oppo(oldtets[1]);

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
}

#endif // #ifndef flipCXX