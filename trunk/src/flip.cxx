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
// flip14()    Insert a vertex by transforming 1-to-4 tetrahedra.            //
//                                                                           //
// 'newpt' (p) lies in the interior of 'fliptets[0]' (abcd).  This routine   //
// removes abcd and replaces it by 4 new tets: abpd, bcpd, capd, and abcp in //
// 'fliptets[0], ..., fliptets[3]', resepctively.                            //
//                                                                           //
// If abcd is a hull tet, we adjust it and let d be the dummypoint. In this  //
// case, the three new tets: abpd, bcpd, and capd are hull tets.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip14(point newpt, triface* fliptets, int flipflag)
{
  triface castets[4]; // The outer boundary faces.
  triface newface, casface;
  point pa, pb, pc, pd;
  int i;

  int *iptr;

  // Check if the original tet is a hull tet.
  if ((point) fliptets[0].tet[7] == dummypoint) {
    fliptets[0].loc = 0;
  }
  fliptets[0].ver = 0;

  pa = org(fliptets[0]);
  pb = dest(fliptets[0]);
  pc = apex(fliptets[0]);
  pd = oppo(fliptets[0]);

  if (b->verbose > 1) {
    printf("    flip 1-to-4: %d (%d, %d, %d, %d)\n", pointmark(newpt), 
      pointmark(pa), pointmark(pb), pointmark(pc), pointmark(pd));
  }
  flip14count++;

  // Get the outer boundary faces.
  for(i = 0; i < 3; i++) {
    fnext(fliptets[0], castets[i]);
    enextself(fliptets[0]);
  }
  symedge(fliptets[0], castets[3]); // At face abc.

  // Delete the old tet.
  tetrahedrondealloc(fliptets[0].tet);

  // Check if d is dummytet.
  if (pd != dummypoint) {
    // Create three new tets for abpd, bcpd, and capd.
    maketetrahedron(tetrahedronpool, &(fliptets[0]));
    maketetrahedron(tetrahedronpool, &(fliptets[1]));
    maketetrahedron(tetrahedronpool, &(fliptets[2]));
  } else {
    // Create three hull tets for abpd, bcpd, and capd.
    maketetrahedron(hulltetrahedronpool, &(fliptets[0]));
    maketetrahedron(hulltetrahedronpool, &(fliptets[1]));
    maketetrahedron(hulltetrahedronpool, &(fliptets[2]));
  }
  // Create the new tet for abcp.
  maketetrahedron(tetrahedronpool, &(fliptets[3]));
  // Set the vertices.
  setvertices(fliptets[0], pa, pb, newpt, pd);
  setvertices(fliptets[1], pb, pc, newpt, pd);
  setvertices(fliptets[2], pc, pa, newpt, pd);
  setvertices(fliptets[3], pa, pb, pc, newpt);

  // Bond the new tets to outer boundary tets.
  for (i = 0; i < 3; i++) {
    enext0fnext(fliptets[i], newface);
    bond(newface, castets[i]);
  }
  bond(fliptets[3], castets[3]); // At face abc.
  // Bond the new tets together (at six interior faces).
  for (i = 0; i < 3; i++) {
    enext0fnext(fliptets[3], newface);
    bond(newface, fliptets[i]);
    enextself(fliptets[3]);
  }
  for (i = 0; i < 3; i++) {
    enextfnext(fliptets[i], newface);
    enext2fnext(fliptets[(i + 1) % 3], casface);
    bond(newface, casface);
  }

  if (flipflag > 0) {
    // Queue faces which may be locally non-Delaunay.
    for (i = 0; i < 3; i++) {
      enext0fnext(fliptets[i], newface);
      futureflip = flippush(futureflip, &newface, newpt);
    }
    futureflip = flippush(futureflip, &(fliptets[3]), newpt);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip26()    Insert a vertex by transforming 2-to-6 tetrahedra.            //
//                                                                           //
// The 'newpt' (p) lies exactly in the face (abc) of the tets 'fliptets[0]'  //
// (abcd) and 'fliptets[1]' (bace).  This routine removes the two tets and   //
// replaces them with 6 new tets: abpd, bcpd, capd (at top), bape, cbpe, and //
// acpe (at bottom), in 'fliptets[0], ..., fliptets[5]', respectively.       //
//                                                                           //
// On input, a, b, or c should not be dummypoint. If d is dummypoint, the 3  //
// top new tets are hull tets.  If e is dummypoint,  we reconfigure e to d,  //
// i.e., turn the two tets up-side down.                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip26(point newpt, triface* fliptets, int flipflag)
{
  triface topcastets[3], botcastets[3]; // The outer boundary faces.
  triface newface, casface;
  point pa, pb, pc, pd, pe;
  int i;

  int *iptr;

  // Check if e is dummypoint.
  if (oppo(fliptets[1]) == dummypoint) {
    // Swap the two old tets.
    newface = fliptets[0];
    fliptets[0] = fliptets[1];
    fliptets[1] = newface;
  }

  pa = org(fliptets[0]);
  pb = dest(fliptets[0]);
  pc = apex(fliptets[0]);
  pd = oppo(fliptets[0]);
  pe = oppo(fliptets[1]);

  if (b->verbose > 1) {
    printf("    flip 2-to-6: (%d, %d, %d, %d, %d)\n", pointmark(pa),
      pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pe));
  }
  flip26count++;

  // Get the outer boundary faces.
  for (i = 0; i < 3; i++) {
    fnext(fliptets[0], topcastets[i]);
    enextself(fliptets[0]);
  }
  for (i = 0; i < 3; i++) {
    fnext(fliptets[1], botcastets[i]);
    enext2self(fliptets[1]);
  }

  // Delete the old tets.
  tetrahedrondealloc(fliptets[0].tet);
  tetrahedrondealloc(fliptets[1].tet);

  // Check if d is dummytet.
  if (pd != dummypoint) {
    maketetrahedron(tetrahedronpool, &(fliptets[0])); // abpd
    maketetrahedron(tetrahedronpool, &(fliptets[1])); // bcpd
    maketetrahedron(tetrahedronpool, &(fliptets[2])); // capd
  } else {
    maketetrahedron(hulltetrahedronpool, &(fliptets[0])); // abpd
    maketetrahedron(hulltetrahedronpool, &(fliptets[1])); // bcpd
    maketetrahedron(hulltetrahedronpool, &(fliptets[2])); // capd
  }
  maketetrahedron(tetrahedronpool, &(fliptets[3])); // bape
  maketetrahedron(tetrahedronpool, &(fliptets[4])); // cbpe
  maketetrahedron(tetrahedronpool, &(fliptets[5])); // acpe
  // Set new vertices.
  setvertices(fliptets[0], pa, pb, newpt, pd);
  setvertices(fliptets[1], pb, pc, newpt, pd);
  setvertices(fliptets[2], pc, pa, newpt, pd);
  setvertices(fliptets[3], pb, pa, newpt, pe);
  setvertices(fliptets[4], pc, pb, newpt, pe);
  setvertices(fliptets[5], pa, pc, newpt, pe);

  // Bond the new tets to outer boundary faces.
  for (i = 0; i < 3; i++) {
    enext0fnext(fliptets[i], newface);
    bond(newface, topcastets[i]);
  }
  for (i = 0; i < 3; i++) {
    enext0fnext(fliptets[i + 3], newface);
    bond(newface, botcastets[i]);
  }
  // Bond the top and bottom new tets together.
  bond(fliptets[0], fliptets[3]);
  bond(fliptets[1], fliptets[4]);
  bond(fliptets[2], fliptets[5]);
  // Bond the new tets together.
  for (i = 0; i < 3; i++) {
    enextfnext(fliptets[i], newface);
    enext2fnext(fliptets[(i + 1) % 3], casface);
    bond(newface, casface);
  }
  for (i = 0; i < 3; i++) {
    enext2fnext(fliptets[3 + i], newface);
    enextfnext(fliptets[3 + (i + 1) % 3], casface);
    bond(newface, casface);
  }

  if (flipflag > 0) {
    for (i = 0; i < 6; i++) {
      enext0fnext(fliptets[i], newface);
      futureflip = flippush(futureflip, &newface, newpt);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipn2n()    Insert a vertex by transforming n-to-2n tetrahedra.          //
//                                                                           //
// The 'newpt'(p) lies on the edge (ab) of 'splitedge'(abcd). Let the n tets //
// containing ab be: abp[0]p[1], ..., abp[n-1]p[0], use an array 'fliptets': //
// 'fliptets[0]', ..., 'fliptets[n-1]', to store them.  This routine removes //
// the n tets and replaces them with 2n new tets: app[0]p[1], ..., app[n-1]- //
// p[0] (at top), pbp[0]p[1], ..., pbp[n-1]p[0] (at bottom) in 'fliptets[0]',//
// ..., 'fliptets[2n-1]', respectively.                                      //
//                                                                           //
// On input, a and b should not be dummypoint. If p[0] is dummypoint, the 2  //
// new tets connecting to p[0] are hull tets. If p[i], i != 0 is dummypoint, //
// we reconfigure p[i] to p[0], i.e., up-shift the tets i times.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flipn2n(point newpt, triface* splitedge, int flipflag)
{
  triface *fliptets, *bfliptets, *topcastets, *botcastets;
  triface newface, casface;
  point pa, pb, *pt;
  int dummyflag; // 0 or 1.
  int n, i;

  int *iptr;

  n = 0;
  dummyflag = 0;

  // First count the number n of tets having ab.
  splitedge->ver &= ~1;
  casface = *splitedge;
  do {
    if (apex(casface) == dummypoint) {
      // Remeber this face (it's apex is dummypoint).
      newface = casface;
      dummyflag = 1;
    }
    n++;
    fnextself(casface);
  } while (casface.tet != splitedge->tet);

  // Allocate spaces for fliptets.
  fliptets = new triface[n];
  bfliptets = new triface[n];
  topcastets = new triface[n];
  botcastets = new triface[n];
  pt = new point[n];

  // Get the n old tets.
  fliptets[0] = dummyflag == 0 ? *splitedge : newface;
  for (i = 0; i < n - 1; i++) {
    pt[i] = apex(fliptets[i]);
    fnext(fliptets[i], fliptets[i + 1]);
  }
  pt[i] = apex(fliptets[i]);
  pa = org(fliptets[0]);
  pb = dest(fliptets[0]);

  if (b->verbose > 1) {
    printf("    flip n-to-2n: (%d, %d) %d %d ..., (n = %d)\n", pointmark(pa),
      pointmark(pb), pointmark(pt[0]), pointmark(pt[1]), n);
  }
  flipn2ncount++;

  // Get the outer boundary faces.
  for (i = 0; i < n; i++) {
    enext2(fliptets[i], casface);  // at edge p[i]a.
    fnext(casface, topcastets[i]);
  }
  for (i = 0; i < n; i++) {
    enext(fliptets[i], casface);  // at edge bp[i].
    fnext(casface, botcastets[i]);
  }

  // Delete the n old tets.
  for (i = 0; i < n; i++) {
    tetrahedrondealloc(fliptets[i].tet);
  }

  // Create 2n new tets.
  if (pt[0] != dummypoint) {
    maketetrahedron(tetrahedronpool, &(fliptets[0])); // app[0]p[1]
    maketetrahedron(tetrahedronpool, &(fliptets[n - 1])); // app[n-1]p[0]
    maketetrahedron(tetrahedronpool, &(bfliptets[0])); // pbp[0]p[1]
    maketetrahedron(tetrahedronpool, &(bfliptets[n - 1])); // pbp[n-1]p[0]
  } else {
    maketetrahedron(hulltetrahedronpool, &(fliptets[0])); // app[0]p[1]
    maketetrahedron(hulltetrahedronpool, &(fliptets[n - 1])); // app[n-1]p[0]
    maketetrahedron(hulltetrahedronpool, &(bfliptets[0])); // pbp[0]p[1]
    maketetrahedron(hulltetrahedronpool, &(bfliptets[n - 1])); // pbp[n-1]p[0]
  }
  for (i = 1; i < n - 1; i++) {
    maketetrahedron(tetrahedronpool, &(fliptets[i])); // app[i]p[i+1]
    maketetrahedron(tetrahedronpool, &(bfliptets[i])); // pbp[i]p[i+1]
  }
  // Set new vertices.
  for (i = 0; i < n; i++) {
    setvertices(fliptets[i], pa, newpt, pt[i], pt[(i + 1) % n]);
    setvertices(bfliptets[i], pb, newpt, pt[i], pt[(i + 1) % n]);
  }

  // Bond new tets to outer boundary faces.
  for (i = 0; i < n; i++) {
    enext2fnext(fliptets[i], newface);
    bond(newface, topcastets[i]);
  }
  for (i = 0; i < n; i++) {
    enextfnext(bfliptets[i], newface);
    bond(newface, botcastets[i]);
  }
  // Bond top and bottom new tets togther
  for (i = 0; i < n; i++) {
    enextfnext(fliptets[i], newface);
    enext2fnext(bfliptets[i], casface);
    bond(newface, casface);
  }
  // Bond new tets together.
  for (i = 0; i < n; i++) {
    enext0fnext(fliptets[i], newface);
    bond(newface, fliptets[(i + 1) % n]);
  }
  for (i = 0; i < n; i++) {
    enext0fnext(bfliptets[i], newface);
    bond(newface, bfliptets[(i + 1) % n]);
  }

  if (flipflag > 0) {
    for (i = 0; i < n; i++) {
      enext2fnext(fliptets[i], newface);
      futureflip = flippush(futureflip, &newface, newpt);
    }
    for (i = 0; i < n; i++) {
      enextfnext(bfliptets[i], newface);
      futureflip = flippush(futureflip, &newface, newpt);
    }
  }

  delete [] fliptets;
  delete [] bfliptets;
  delete [] topcastets;
  delete [] botcastets;
  delete [] pt;
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

void tetgenmesh::flip23(triface* fliptets, int flipflag)
{
  triface topcastets[3], botcastets[3]; // The outer boundary faces.
  triface newface, casface;
  point pa, pb, pc, pd, pe;
  int dummyflag;  // range = {-1, 0, 1, 2}.
  int i;

  tetrahedron ptr;
  int *iptr;

  // Check if e is dummypoint.
  if (oppo(fliptets[1]) == dummypoint) {
    // Swap the two old tets.
    newface = fliptets[0];
    fliptets[0] = fliptets[1];
    fliptets[1] = newface;
    dummyflag = -1;  // e is dummypoint.
  } else {
    // Check if a or b is dummypoint.
    if (org(fliptets[0]) == dummypoint) {
      dummyflag = 1;  // a is dummypoint.
    } else if (dest(fliptets[0]) == dummypoint) {
      dummyflag = 2;  // b is dummypoint.
    } else {
      dummyflag = 0;  // either c or d is dummypoint.
    }
    i = dummyflag;
    // Rotate i times.
    for (; i > 0; i--) {
      enextself(fliptets[0]);
      enext2self(fliptets[1]);
    }
  }

  pa = org(fliptets[0]);
  pb = dest(fliptets[0]);
  pc = apex(fliptets[0]);
  pd = oppo(fliptets[0]);
  pe = oppo(fliptets[1]);

  if (b->verbose > 1) {
    printf("    flip 2-to-3: (%d, %d, %d, %d, %d)\n", pointmark(pa),
      pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pe));
  }
  flip23count++;

  // Get the outer boundary faces.
  for (i = 0; i < 3; i++) {
    fnext(fliptets[0], topcastets[i]);
    enextself(fliptets[0]);
  }
  for (i = 0; i < 3; i++) {
    fnext(fliptets[1], botcastets[i]);
    enext2self(fliptets[1]);
  }

  // Delete the old tets.
  tetrahedrondealloc(fliptets[0].tet);
  tetrahedrondealloc(fliptets[1].tet);

  // Check if d is dummytet.
  if (pd != dummypoint) {
    maketetrahedron(tetrahedronpool, &(fliptets[0])); // Create edab.
    setvertices(fliptets[0], pe, pd, pa, pb);
    // Check if c is dummypoint.
    if (pc != dummypoint) {
      maketetrahedron(tetrahedronpool, &(fliptets[1])); // Create edbc.
      maketetrahedron(tetrahedronpool, &(fliptets[2])); // Create edca.
      setvertices(fliptets[2], pe, pd, pc, pa);
    } else {
      maketetrahedron(hulltetrahedronpool, &(fliptets[1])); // Create edbc.
      maketetrahedron(hulltetrahedronpool, &(fliptets[2])); // Create deac.
      setvertices(fliptets[2], pd, pe, pa, pc);
      // Adjust deac->edca
      enext0fnextself(fliptets[2]);
      esymself(fliptets[2]);
    }
    setvertices(fliptets[1], pe, pd, pb, pc);
  } else {
    // d is dummypoint. Create three new hull tets.
    maketetrahedron(hulltetrahedronpool, &(fliptets[0])); // Create abed.
    maketetrahedron(hulltetrahedronpool, &(fliptets[1])); // Create bced.
    maketetrahedron(hulltetrahedronpool, &(fliptets[2])); // Create caed.
    setvertices(fliptets[0], pa, pb, pe, pd);
    setvertices(fliptets[1], pb, pc, pe, pd);
    setvertices(fliptets[2], pc, pa, pe, pd);
    // Adjust abed->edab, bced->edbc, caed->edca
    for (i = 0; i < 3; i++) {
      enext2fnextself(fliptets[i]);
      enext2self(fliptets[i]);
      esymself(fliptets[i]);
    }
  }

  // Bond three new tets together.
  for (i = 0; i < 3; i++) {
    enext0fnext(fliptets[i], newface);
    bond(newface, fliptets[(i + 1) % 3]);
  }
  // Bond to top outer boundary faces (at abcd).
  for (i = 0; i < 3; i++) {
    enextfnext(fliptets[i], newface);
    enextself(newface);
    bond(newface, topcastets[i]);
  }
  // Bond bottom outer boundary faces (at bace).
  for (i = 0; i < 3; i++) {
    enext2fnext(fliptets[i], newface);
    enext2self(newface);
    bond(newface, botcastets[i]);
  }
  
  if (dummyflag != 0) {
    // Restore the original position of the points (for flipnm()).
    if (dummyflag == -1) { 
      // Reverse the edge.
      for (i = 0; i < 3; i++) {
        enext0fnextself(fliptets[i]);
        esymself(fliptets[i]);
      }
      // Swap the last two new tets.
      newface = fliptets[1];
      fliptets[1] = fliptets[2];
      fliptets[2] = newface;
    } else {
      // either a or b were swapped.
      i = dummyflag;
      // Down-shift new tets i times.
      for (; i > 0; i--) {
        newface = fliptets[0];
        fliptets[0] = fliptets[2];
        fliptets[2] = fliptets[1];
        fliptets[1] = newface;
      }
    }
  }

  if (flipflag > 0) {
    // Queue faces which may be locally non-Delaunay.  
    pd = dest(fliptets[0]);
    for (i = 0; i < 3; i++) {
      if ((point) botcastets[i].tet[7] != dummypoint) {
        enext2fnext(fliptets[i], newface);
        futureflip = flippush(futureflip, &newface, pd);
      }
    }
    if (flipflag > 1) {
      pe = org(fliptets[0]);
      for (i = 0; i < 3; i++) {
        if ((point) topcastets[i].tet[7] != dummypoint) {
          enextfnext(fliptets[i], newface);
          futureflip = flippush(futureflip, &newface, pe);
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

void tetgenmesh::flip32(triface* fliptets, int flipflag)
{
  triface topcastets[3], botcastets[3]; // The outer boundary faces.
  triface newface, casface;
  point pa, pb, pc, pd, pe;
  int dummyflag;  // Rangle = {-1, 0, 1, 2}
  int i;

  tetrahedron ptr;
  int *iptr; 

  // Check if e is 'dummypoint'.
  if (org(fliptets[0]) == dummypoint) {
    // Reverse the edge.
    for (i = 0; i < 3; i++) {
      enext0fnextself(fliptets[i]);
      esymself(fliptets[i]);
    }
    // Swap the last two tets.
    newface = fliptets[1];
    fliptets[1] = fliptets[2];
    fliptets[2] = newface;
    dummyflag = -1; // e is dummypoint.
  } else {
    // Check if a or b is the 'dummypoint'.
    if (apex(fliptets[0]) == dummypoint) { 
      dummyflag = 1;  // a is dummypoint.
    } else if (apex(fliptets[1]) == dummypoint) {
      dummyflag = 2;  // b is dummypoint.
    } else {
      dummyflag = 0;  // either c or d is dummypoint.
    }
    i = dummyflag;
    // Down-shift i times.
    for (; i > 0; i--) {
      newface = fliptets[2];
      fliptets[2] = fliptets[1];
      fliptets[1] = fliptets[0];
      fliptets[0] = newface;
    }
  }

  pa = apex(fliptets[0]);
  pb = apex(fliptets[1]);
  pc = apex(fliptets[2]);
  pd = dest(fliptets[0]);
  pe = org(fliptets[0]);

  if (b->verbose > 1) {
    printf("    flip 3-to-2: (%d, %d, %d, %d, %d)\n", pointmark(pa),
      pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pe));
  }
  flip32count++;

  // Get the outer boundary faces.
  for (i = 0; i < 3; i++) {
    enextfnext(fliptets[i], casface);
    enextself(casface);
    symedge(casface, topcastets[i]);
  }
  for (i = 0; i < 3; i++) {
    enext2fnext(fliptets[i], casface);
    enext2self(casface);
    symedge(casface, botcastets[i]);
  }

  // Delete the old tets.
  tetrahedrondealloc(fliptets[0].tet);
  tetrahedrondealloc(fliptets[1].tet);
  tetrahedrondealloc(fliptets[2].tet);

  // Check if c is a dummypointc.
  if (pc != dummypoint) {
    // Create abcd. Check if d is a dummypoint.
    if (pd != dummypoint) {
      maketetrahedron(tetrahedronpool, &(fliptets[0]));
    } else {
      maketetrahedron(hulltetrahedronpool, &(fliptets[0])); // a hull tet.
    }
    setvertices(fliptets[0], pa, pb, pc, pd);
    // Create bace.
    maketetrahedron(tetrahedronpool, &(fliptets[1]));
    setvertices(fliptets[1], pb, pa, pc, pe);
  } else {
    // c is dummypoint. The two new tets are hull tets.
    // Create badc.
    maketetrahedron(hulltetrahedronpool, &(fliptets[0]));
    setvertices(fliptets[0], pb, pa, pd, pc);
    // Create abec
    maketetrahedron(hulltetrahedronpool, &(fliptets[1]));
    setvertices(fliptets[1], pa, pb, pe, pc);
    // Adjust badc -> abcd.
    enext0fnextself(fliptets[0]);
    esymself(fliptets[0]);
    // Adjust abec -> bace.
    enext0fnextself(fliptets[1]);
    esymself(fliptets[1]);
  }
  
  // Bond abcd <==> bace.
  bond(fliptets[0], fliptets[1]);
  // Bond new faces to top outer boundary faces (at abcd).
  for (i = 0; i < 3; i++) {
    enext0fnext(fliptets[0], newface);
    bond(newface, topcastets[i]);
    enextself(fliptets[0]);
  }
  // Bond new faces to bottom outer boundary faces (at bace).
  for (i = 0; i < 3; i++) {
    enext0fnext(fliptets[1], newface);
    bond(newface, botcastets[i]);
    enext2self(fliptets[1]);
  }

  if (dummyflag != 0) {
    // Restore the original position of the points (for flipnm()).
    if (dummyflag == -1) {
      // e were dummypoint. Swap the two new tets.
      newface = fliptets[0];
      fliptets[0] = fliptets[1];
      fliptets[1] = newface;
    } else {
      // a or b was dummypoint.
      i = dummyflag;
      // Rotate toward left i times.
      for (; i > 0; i--) {
        enextself(fliptets[0]);
        enext2self(fliptets[1]);
      }
    }
  }
  
  if (flipflag > 0) {
    // Queue faces which may be locally non-Delaunay.  
    pa = org(fliptets[0]);
    if ((point) topcastets[1].tet[7] != dummypoint) {
      enextfnext(fliptets[0], newface);
      futureflip = flippush(futureflip, &newface, pa);
    }
    if ((point) botcastets[1].tet[7] != dummypoint) {
      enext2fnext(fliptets[1], newface);
      futureflip = flippush(futureflip, &newface, pa);
    }
    if (flipflag > 1) {
      pb = dest(fliptets[0]);
      if ((point) topcastets[2].tet[7] != dummypoint) {
        enext2fnext(fliptets[0], newface);
        futureflip = flippush(futureflip, &newface, pb);
      }
      if ((point) botcastets[2].tet[7] != dummypoint) {
        enextfnext(fliptets[1], newface);
        futureflip = flippush(futureflip, &newface, pb);
      }
      pc = apex(fliptets[0]);
      if ((point) topcastets[0].tet[7] != dummypoint) {
        enext0fnext(fliptets[0], newface);
        futureflip = flippush(futureflip, &newface, pc);
      }
      if ((point) botcastets[0].tet[7] != dummypoint) {
        enext0fnext(fliptets[1], newface);
        futureflip = flippush(futureflip, &newface, pc);
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