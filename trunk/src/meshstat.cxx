#ifndef meshstatCXX
#define meshstatCXX

#include "tetgen.h"

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

#endif // #ifndef meshstatCXX