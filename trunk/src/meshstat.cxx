#ifndef meshstatCXX
#define meshstatCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Initialize fast look-up tables.                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::ve[6] = { 2, 5, 4, 1, 0, 3 };

int tetgenmesh::locver2org[4][6]  = {
  {0, 1, 1, 2, 2, 0},
  {0, 3, 3, 1, 1, 0},
  {1, 3, 3, 2, 2, 1},
  {2, 3, 3, 0, 0, 2} 
};

int tetgenmesh::locver2dest[4][6] = { 
  {1, 0, 2, 1, 0, 2},
  {3, 0, 1, 3, 0, 1},
  {3, 1, 2, 3, 1, 2},
  {3, 2, 0, 3, 2, 0}
};

int tetgenmesh::locver2apex[4][6] = { 
  {2, 2, 0, 0, 1, 1},
  {1, 1, 0, 0, 3, 3},
  {2, 2, 1, 1, 3, 3},
  {0, 0, 2, 2, 3, 3}
};

int tetgenmesh::loc2oppo[4] = { 3, 2, 0, 1 };

int tetgenmesh::locver2nextf[24] = {
  1, 5, 2, 5, 3, 5,
  3, 3, 2, 1, 0, 1,
  1, 3, 3, 1, 0, 3,
  2, 3, 1, 1, 0, 5
};

int tetgenmesh::locver2edge[4][6] = {
  {0, 0, 1, 1, 2, 2},
  {3, 3, 4, 4, 0, 0},
  {4, 4, 5, 5, 1, 1},
  {5, 5, 3, 3, 2, 2}
};

int tetgenmesh::edge2locver[6][2] = {
  {0, 0}, // 0  v0 -> v1
  {0, 2}, // 1  v1 -> v2
  {0, 4}, // 2  v2 -> v0
  {1, 0}, // 3  v0 -> v3
  {1, 2}, // 4  v1 -> v3
  {2, 2}  // 5  v2 -> v3
};

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
          horrors++;
        }
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
    
    // Report memory usages.
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

    // Report algorithmic performances.
    printf("Algorithmic statistics:\n\n");
    printf("  Maximal point location distance (# tets): %ld\n", 
           ptloc_max_tets_count);
    printf("  Number of orient3d tests: %ld\n", orient3dcount);
    printf("  Number of insphere tests: %ld\n", inspherecount);
    printf("  Number of symbolic insphere tests: %ld\n", insphere_sos_count);
    if (b->bowyerwatson == 0) {
      printf("  Number of 2-to-3 flips: %ld\n", flip23count);
      printf("  Number of 3-to-2 flips: %ld\n", flip32count);
      printf("  Number of n-to-m flips: %ld\n", flipnmcount);
      printf("  Number of total primitive flips: %ld\n",
             flip23count + flip32count);
    }
    printf("\n");
  }
}

#endif // #ifndef meshstatCXX