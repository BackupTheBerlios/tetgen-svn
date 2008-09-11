#ifndef meshstatCXX
#define meshstatCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// printtet()    Print the detail contents of a tetrahedron.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::printtet(triface* t)
{
  triface tmpface, prtface;
  point tmppt;
  int facecount;

  printf("Tetra x%lx with loc(%i) and ver(%i):",
         (unsigned long)(t->tet), t->loc, t->ver);
  if ((point) t->tet[7] == dummypoint) {
    printf("  (hull tet)");
  }
  if (infected(*t)) {
    printf(" (infected)");
  }
  if (marktested(*t)) {
    printf(" (marktested)");
  }
  printf("\n");

  tmpface = *t;
  facecount = 0;
  while(facecount < 4) {
    tmpface.loc = facecount;
    sym(tmpface, prtface);
    if (prtface.tet == NULL) {
      printf("      [%i] Open !!\n", facecount);
    } else {
      printf("      [%i] x%lx  loc(%i).", facecount,
             (unsigned long)(prtface.tet), prtface.loc);
      if ((point) prtface.tet[7] == dummypoint) {
        printf("  (hull tet)");
      }
      if (infected(prtface)) {
        printf(" (infected)");
      }
      if (marktested(prtface)) {
        printf(" (marktested)");
      }
      printf("\n");
    }
    facecount++;
  }

  tmppt = org(*t);
  if(tmppt == (point) NULL) {
    printf("      Org [%i] NULL\n", locver2org[t->loc][t->ver]);
  } else {
    printf("      Org [%i] x%lx (%.12g,%.12g,%.12g) %d\n",
           locver2org[t->loc][t->ver], (unsigned long)(tmppt),
           tmppt[0], tmppt[1], tmppt[2], pointmark(tmppt));
  }
  tmppt = dest(*t);
  if(tmppt == (point) NULL) {
    printf("      Dest[%i] NULL\n", locver2dest[t->loc][t->ver]);
  } else {
    printf("      Dest[%i] x%lx (%.12g,%.12g,%.12g) %d\n",
           locver2dest[t->loc][t->ver], (unsigned long)(tmppt),
           tmppt[0], tmppt[1], tmppt[2], pointmark(tmppt));
  }
  tmppt = apex(*t);
  if(tmppt == (point) NULL) {
    printf("      Apex[%i] NULL\n", locver2apex[t->loc][t->ver]);
  } else {
    printf("      Apex[%i] x%lx (%.12g,%.12g,%.12g) %d\n",
           locver2apex[t->loc][t->ver], (unsigned long)(tmppt),
           tmppt[0], tmppt[1], tmppt[2], pointmark(tmppt));
  }
  tmppt = oppo(*t);
  if(tmppt == (point) NULL) {
    printf("      Oppo[%i] NULL\n", loc2oppo[t->loc]);
  } else {
    printf("      Oppo[%i] x%lx (%.12g,%.12g,%.12g) %d\n",
           loc2oppo[t->loc], (unsigned long)(tmppt),
           tmppt[0], tmppt[1], tmppt[2], pointmark(tmppt));
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkmesh()    Test mesh for geometrical and topological consistencies.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::checkmesh(memorypool* pool)
{
  triface tetloop;
  triface neightet, symtet;
  point pa, pb, pc, pd;
  REAL ori;
  int horrors;

  if (!b->quiet) {
    printf("  Checking consistency of mesh...\n");
  }

  horrors = 0;
  // Run through the list of tetrahedra, checking each one.
  pool->traversalinit();
  tetloop.tet = tetrahedrontraverse(pool);
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Check all four faces of the tetrahedron.
    for (tetloop.loc = 0; tetloop.loc < 4; tetloop.loc++) {
      pa = org(tetloop);
      pb = dest(tetloop);
      pc = apex(tetloop);
      pd = oppo(tetloop);
      if (tetloop.loc == 0) {  // Only test for inversion once.
        if (pd != dummypoint) {  // Only do test if it is not a hull tet.
          ori = orient3d(pa, pb, pc, pd);
          if (ori >= 0.0) {
            printf("  !! !! %s ", ori > 0.0 ? "Inverted" : "Degenerated");
            printf("  (%d, %d, %d, %d) (ori = %.17g)\n", pointmark(pa),
              pointmark(pb), pointmark(pc), pointmark(pd), ori);
            horrors++;
          }
        }
      }
      if (tetloop.tet[tetloop.loc] == NULL) {
        printf("  !! !! No neighbor at face (%d, %d, %d).\n", pointmark(pa),
            pointmark(pb), pointmark(pc));
        horrors++;
      } else {
        // Find the neighboring tetrahedron on this face.
        sym(tetloop, neightet);
        // Check that the tetrahedron's neighbor knows it's a neighbor.
        sym(neightet, symtet);
        if ((tetloop.tet != symtet.tet) || (tetloop.loc != symtet.loc)) {
          printf("  !! !! Asymmetric tetra-tetra bond:\n");
          if (tetloop.tet == symtet.tet) {
            printf("   (Right tetrahedron, wrong orientation)\n");
          }
          printf("    First:  (%d, %d, %d, %d)\n", pointmark(pa),
            pointmark(pb), pointmark(pc), pointmark(pd));
          printf("    Second: (%d, %d, %d, %d)\n", pointmark(org(neightet)),
            pointmark(dest(neightet)), pointmark(apex(neightet)),
            pointmark(oppo(neightet)));
          horrors++;
        }
        // Check bond() operation.
        tetloop.ver = 0; // Fix at the 0th edge.
        pa = org(tetloop);
        pb = dest(tetloop);
        // The edge of neightet should symmetric to this.
        if ((org(neightet) != pb) || (dest(neightet) != pa)) {
          printf("  !! !! Wrong edge-edge bond:\n");
          printf("    First:  (%d, %d, %d, %d)\n", pointmark(pa),
            pointmark(pb), pointmark(pc), pointmark(pd));
          printf("    Second: (%d, %d, %d, %d)\n", pointmark(org(neightet)),
            pointmark(dest(neightet)), pointmark(apex(neightet)),
            pointmark(oppo(neightet)));
          horrors++;
        }        
      }
    }
    tetloop.tet = tetrahedrontraverse(pool);
  }
  if (horrors == 0) {
    if (!b->quiet) {
      printf("  In my studied opinion, the mesh appears to be consistent.\n");
    }
  } else {
    printf("  !! !! !! !! %d %s witnessed.\n", horrors, 
      horrors > 1 ? "abnormity" : "abnormities");
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkdelaunay()    Ensure that the mesh is (constrained) Delaunay.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::checkdelaunay()
{
  triface tetloop;
  triface symtet;
  point pa, pb, pc, pd, pe;
  REAL sign;
  int horrors;

  if (!b->quiet) {
    printf("  Checking Delaunay property of the mesh...\n");
  }
  
  horrors = 0;
  // Run through the list of triangles, checking each one.
  tetrahedronpool->traversalinit();
  tetloop.tet = tetrahedrontraverse(tetrahedronpool);
  while (tetloop.tet != (tetrahedron *) NULL) {
    pa = (point) tetloop.tet[4];
    pb = (point) tetloop.tet[5];
    pc = (point) tetloop.tet[6];
    pd = (point) tetloop.tet[7];
    // Check all four faces of the tetrahedron.
    for (tetloop.loc = 0; tetloop.loc < 4; tetloop.loc++) {
      sym(tetloop, symtet);
      // Only do test if its adjoining tet is not a hull tet or its pointer
      //   is larger (to ensure that each pair isn't tested twice).
      if (((point) symtet.tet[7] != dummypoint)&&(tetloop.tet < symtet.tet)) {
        pe = oppo(symtet);
        sign = insphere_sos(pa, pb, pc, pd, pe);
        if (sign < 0.0) {
          printf("  !! Non-locally Delaunay (%d, %d, %d) - %d, %d\n",
            pointmark(pa), pointmark(pb), pointmark(pc), pointmark(pd),
            pointmark(pe));
        }
        horrors++;
      }
    }
    tetloop.tet = tetrahedrontraverse(tetrahedronpool);
  }

  if (horrors == 0) {
    if (!b->quiet) {
      printf("  The mesh is Delaunay.\n");
    }
  } else {
    printf("  !! !! !! !! Found %d non-Delaunay faces.\n", horrors);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// statistics()    Print all sorts of cool facts.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::statistics()
{
  printf("\nStatistics:\n\n");
  printf("  Input points: %d\n", in->numberofpoints);
  if (b->refine) {
    printf("  Input tetrahedra: %d\n", in->numberoftetrahedra);
  }
  if (b->plc) {
    printf("  Input facets: %d\n", in->numberoffacets);
    printf("  Input segments: %ld\n", insegments);
    printf("  Input holes: %d\n", in->numberofholes);
    printf("  Input regions: %d\n", in->numberofregions);
  }

  printf("\n  Mesh points: %ld\n", pointpool->items);
  printf("  Mesh tetrahedra: %ld\n", tetrahedronpool->items);
  printf("  Mesh faces: %ld\n", (tetrahedronpool->items * 4 + 
    hulltetrahedronpool->items) / 2l);
  printf("  Mesh edges: %ld\n", meshedges);

  if (b->plc || b->refine) {
    // printf("  Mesh boundary faces: %ld\n", subfaces->items);
    // printf("  Mesh boundary edges: %ld\n\n", subsegs->items);
  } else {
    printf("  Convex hull faces: %ld\n\n", hulltetrahedronpool->items);
  }
  if (b->verbose > 0) {
    // qualitystatistics();
    unsigned long totalmeshbytes;
    printf("Memory allocation statistics:\n\n");
    printf("  Maximum number of vertices: %ld\n", pointpool->maxitems);
    totalmeshbytes = pointpool->maxitems * pointpool->itembytes;
    printf("  Maximum number of tetrahedra: %ld\n", tetrahedronpool->maxitems);
    totalmeshbytes += tetrahedronpool->maxitems * tetrahedronpool->itembytes;
    //if (subfaces != (memorypool *) NULL) {
    //  printf("  Maximum number of subfaces: %ld\n", subfaces->maxitems);
    //  totalmeshbytes += subfaces->maxitems * subfaces->itembytes;
    //}
    //if (subsegs != (memorypool *) NULL) {
    //  printf("  Maximum number of segments: %ld\n", subsegs->maxitems);
    //  totalmeshbytes += subsegs->maxitems * subsegs->itembytes;
    //}
    printf("  Heap memory used by the mesh (K bytes): %g\n\n",
           (double) totalmeshbytes / 1024.0);
  }
}

#endif // #ifndef meshstatCXX