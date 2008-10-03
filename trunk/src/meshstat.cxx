#ifndef meshstatCXX
#define meshstatCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Initialize fast look-up tables.                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::ve[6] = { 2, 5, 4, 1, 0, 3 };
int tetgenmesh::vo[6] = { 0, 1, 1, 2, 2, 0 };
int tetgenmesh::vd[6] = { 1, 0, 2, 1, 0, 2 };
int tetgenmesh::va[6] = { 2, 2, 0, 0, 1, 1 };

int tetgenmesh::verver2zero[6][6] = {
  {0, 0, 2, 2, 4, 4},
  {0, 0, 2, 2, 4, 4},
  {2, 2, 4, 4, 0, 0},
  {2, 2, 4, 4, 0, 0},
  {4, 4, 0, 0, 2, 2},
  {4, 4, 0, 0, 2, 2}
};

int tetgenmesh::ver2zero[6] = {0, 0, 2, 2, 4, 4};

int tetgenmesh::zero2ver[6][6] = {
  {0, 0, 2, 2, 4, 4}, 
  {0, 0, 2, 2, 4, 4},
  {4, 4, 0, 0, 2, 2},
  {4, 4, 0, 0, 2, 2},
  {2, 2, 4, 4, 0, 0},
  {2, 2, 4, 4, 0, 0}
};

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

int tetgenmesh::locver2nextf[32] = {
  1, 5, 2, 5, 3, 5, 0, 0,
  3, 3, 2, 1, 0, 1, 0, 0,
  1, 3, 3, 1, 0, 3, 0, 0,
  2, 3, 1, 1, 0, 5, 0, 0
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

int tetgenmesh::locpivot[4][3] = {
  {1, 2, 3},
  {0, 2, 3},
  {0, 1, 3},
  {0, 1, 2}
};

int tetgenmesh::locverpivot[4][6][2] = {
  {{2, 3}, {2, 3}, {1, 3}, {1, 3}, {1, 2}, {1, 2}},
  {{0, 2}, {0, 2}, {0, 3}, {0, 3}, {2, 3}, {2, 3}},
  {{0, 3}, {0, 3}, {0, 1}, {0, 1}, {1, 3}, {1, 3}},
  {{0, 1}, {0, 1}, {0, 2}, {0, 2}, {1, 2}, {1, 2}}
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
  tetloop.ver = 0;
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
        symedge(tetloop, neightet);
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
        // Check if they have the same edge (the bond() operation).
        if ((org(neightet) != pb) || (dest(neightet) != pa)) {
          printf("  !! !! Wrong edge-edge bond:\n");
          printf("    First:  (%d, %d, %d, %d)\n", pointmark(pa),
            pointmark(pb), pointmark(pc), pointmark(pd));
          printf("    Second: (%d, %d, %d, %d)\n", pointmark(org(neightet)),
            pointmark(dest(neightet)), pointmark(apex(neightet)),
            pointmark(oppo(neightet)));
          horrors++;
        }
        // Check if they have the same opposite.
        if (oppo(neightet) == pd) {
          printf("  !! !! Two identical tetra:\n");
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
// checkshells()    Test the surface mesh for consistencies.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::checkshells()
{
  face shloop, spinsh, nextsh;
  point pa, pb;
  int horrors, i;

  if (!b->quiet) {
    printf("  Checking consistency of the mesh boundary...\n");
  }
  horrors = 0;

  subfacepool->traversalinit();
  shloop.sh = shellfacetraverse(subfacepool);
  while (shloop.sh != NULL) {
    shloop.shver = 0;
    for (i = 0; i < 3; i++) {
      // Check the face ring at this edge.
      pa = sorg(shloop);
      pb = sdest(shloop);
      spinsh = shloop;
      spivot(spinsh, nextsh);
      while (nextsh.sh != shloop.sh) {
        // check if they have the same edge.
        if (!((sorg(nextsh) == pa) && (sdest(nextsh) == pb) ||
             (sorg(nextsh) == pb) && (sdest(nextsh) == pa))) {
           printf("  !! !! Wrong subface-subface connection.\n");
           printf("    First: x%lx (%d, %d, %d).\n", (unsigned long) spinsh.sh,
             pmark(sorg(spinsh)), pmark(sdest(spinsh)), pmark(sapex(spinsh)));
           printf("    Scond: x%lx (%d, %d, %d).\n", (unsigned long) nextsh.sh,
             pmark(sorg(nextsh)), pmark(sdest(nextsh)), pmark(sapex(nextsh)));
           horrors++;
        }
        spinsh = nextsh;
        spivot(spinsh, nextsh);
      }
      senextself(shloop);
    }
    shloop.sh = shellfacetraverse(subfacepool);
  }

  if (horrors == 0) {
    if (!b->quiet) {
      printf("  Mesh boundaries connected correctly.\n");
    }
  } else {
    printf("  !! !! !! !! %d boundary connection viewed with horror.\n",
           horrors);
    return;
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
  tetloop.ver = 0;
  // Run through the list of triangles, checking each one.
  tetrahedronpool->traversalinit();
  tetloop.tet = tetrahedrontraverse(tetrahedronpool);
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Check all four faces of the tetrahedron.
    for (tetloop.loc = 0; tetloop.loc < 4; tetloop.loc++) {
      sym(tetloop, symtet);
      // Only do test if its adjoining tet is not a hull tet or its pointer
      //   is larger (to ensure that each pair isn't tested twice).
      if (((point) symtet.tet[7] != dummypoint)&&(tetloop.tet < symtet.tet)) {
        pa = org(tetloop);
        pb = dest(tetloop);
        pc = apex(tetloop);
        pd = oppo(tetloop);
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
    printf("  Mesh boundary faces: %ld\n", subfacepool->items);
    printf("  Mesh boundary edges: %ld\n", subsegpool->items);
  } else {
    printf("  Convex hull faces: %ld\n", hulltetrahedronpool->items);
  }
  printf("\n");

  if (b->verbose > 0) {
    // qualitystatistics();
    // Report memory usages.
    unsigned long totalmeshbytes;
    printf("Memory allocation statistics:\n\n");
    printf("  Maximum number of vertices: %ld\n", pointpool->maxitems);
    totalmeshbytes = pointpool->maxitems * pointpool->itembytes;
    printf("  Maximum number of tetrahedra: %ld\n", tetrahedronpool->maxitems);
    totalmeshbytes += tetrahedronpool->maxitems * tetrahedronpool->itembytes;
    if (subfacepool != (memorypool *) NULL) {
      printf("  Maximum number of subfaces: %ld\n", subfacepool->maxitems);
      totalmeshbytes += subfacepool->maxitems * subfacepool->itembytes;
    }
    if (subsegpool != (memorypool *) NULL) {
      printf("  Maximum number of segments: %ld\n", subsegpool->maxitems);
      totalmeshbytes += subsegpool->maxitems * subsegpool->itembytes;
    }
    printf("  Heap memory used by the mesh (K bytes): %g\n\n",
           ((double) totalmeshbytes) / 1024.0);

    // Report algorithmic performances.
    printf("Algorithmic statistics:\n\n");
    printf("  Number of orient3d tests: %ld\n", orient3dcount);
    printf("  Number of insphere tests: %ld\n", inspherecount);
    printf("  Number of symbolic insphere tests: %ld\n", insphere_sos_count);
    if (b->bowyerwatson) {
      printf("  Total size of cavities (# tets): %ld\n", totalbowatcavsize);
      printf("  Maximal size of a single cavity: %ld\n", maxbowatcavsize);
    } else {
      printf("  Number of 1-to-4 flips: %ld\n", flip14count);
      printf("  Number of 2-to-6 flips: %ld\n", flip26count);
      printf("  Number of n-t-2n flips: %ld\n", flipn2ncount);
      printf("  Number of n-to-m flips: %ld\n", flipnmcount);
      printf("  Number of 2-to-3 flips: %ld\n", flip23count);
      printf("  Number of 3-to-2 flips: %ld\n", flip32count);
      printf("  Number of total primitive flips: %ld\n",
             flip23count + flip32count);
    }
    if (b->plc) {
      printf("  Number of 2-to-2 flips: %ld\n", flip22count);
    }
    printf("  Total point location distance (# tets): %ld\n", ptloc_count);
    printf("  Maximal point location distance (# tets): %ld\n", 
           ptloc_max_count);
    // printf("  Total point location time (millisec):  %g\n", 
    //        tloctime * 1e+3);
    // printf("  Total point insertion time (millisec):  %g\n",
    //        tinserttime * 1e+3);
    // if (b->bowyerwatson == 0) {
    //   printf("  Total flip time (millisec):  %g\n", tfliptime * 1e+3);
    // }
    printf("\n");
  }
}

#endif // #ifndef meshstatCXX