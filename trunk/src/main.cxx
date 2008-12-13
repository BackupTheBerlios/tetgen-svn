#ifndef mainCXX
#define mainCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// test_tri_tri()    Triangle-triangle intersection tests.                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void test_tri_tri(tetgenbehavior *b, tetgenio *in)
{
  tetgenmesh m;
  tetgenmesh::face s1, s2;
  tetgenmesh::point U[3], V[3], T1[3], T2[3];
  tetgenmesh::intersection dir;
  int types[2], pos[4];
  int i, j;

  // For checking results.
  tetgenmesh::point pa, pb;
  int horror = 0;

  m.b = b;
  m.in = in;
  exactinit();
  m.initializepools();
  m.transfernodes();
  m.meshsurface();

  m.subfacepool->traversalinit();
  s1.sh = m.shellfacetraverse(m.subfacepool);
  s2.sh = m.shellfacetraverse(m.subfacepool);

  V[0] = (tetgenmesh::point) s1.sh[3]; // P
  V[1] = (tetgenmesh::point) s1.sh[4]; // Q
  V[2] = (tetgenmesh::point) s1.sh[5]; // R

  U[0] = (tetgenmesh::point) s2.sh[3]; // A
  U[1] = (tetgenmesh::point) s2.sh[4]; // B
  U[2] = (tetgenmesh::point) s2.sh[5]; // C

  // The permuations of {0, 1, 2}
  int perm[6][3] = { 
    {0, 1, 2}, 
    {2, 0, 1}, 
    {1, 2, 0},
    {1, 0, 2}, 
    {2, 1, 0}, 
    {0, 2, 1}};

  b->epsilon = 0;
  b->verbose = 3;

  for (j = 0; j < 36; j++) {

    // Get the permuation of the vertices.
    T1[0] = U[perm[j/6][0]];
    T1[1] = U[perm[j/6][1]];
    T1[2] = U[perm[j/6][2]];
    
    T2[0] = V[perm[j%6][0]];
    T2[1] = V[perm[j%6][1]];
    T2[2] = V[perm[j%6][2]];

    i = m.tri_tri_test(T1[0], T1[1], T1[2], T2[0], T2[1], T2[2], NULL, 1, 
                       types, pos);
    if (i == 0) {
      printf("  [%d] DISJOINT.\n", j);
      continue;
    }

    // Report the intersection types and positions.
    for (i = 0; i < 2; i++) {
      dir = (enum tetgenmesh::intersection) types[i];
      switch (dir) {
      case tetgenmesh::DISJOINT: 
        printf("  [%d] DISJOINT\n", j); break;
      case tetgenmesh::SHAREVERT: 
        printf("  [%d] SHAREVERT %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      case tetgenmesh::SHAREEDGE: 
        printf("  [%d] SHAREEDGE %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      case tetgenmesh::SHAREFACE: 
        printf("  [%d] SHAREFACE\n", j); break;
      case tetgenmesh::TOUCHEDGE: 
        printf("  [%d] TOUCHEDGE %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      case tetgenmesh::TOUCHFACE: 
        printf("  [%d] TOUCHFACE %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      case tetgenmesh::ACROSSVERT: 
        printf("  [%d] ACROSSVERT %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      case tetgenmesh::ACROSSEDGE: 
        printf("  [%d] ACROSSEDGE %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      case tetgenmesh::ACROSSFACE: 
        printf("  [%d] ACROSSFACE %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      case tetgenmesh::ACROSSTET: 
        printf("  [%d] ACROSSTET\n", j); break;
      case tetgenmesh::TRIEDGEINT:
        printf("  [%d] TRIEDGEINT %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      case tetgenmesh::EDGETRIINT:
        printf("  [%d] EDGETRIINT %d %d\n", j, pos[i*2], pos[i*2+1]); break;
      }
    }

    horror = 0;

    // Check the correctness of this test. Insert the code here
    if (types[0] != (int) tetgenmesh::TOUCHEDGE) {
      printf("  !! Wrong case.\n");
      horror++;
    }
    //<<<<<<<<<<<<< The first intersection 
    /*if (pos[0] != 3) {
      printf("  !! Wrong report\n");
      horror++;
    }*/
    pa = T1[pos[0]];
    pb = T1[(pos[0] + 1) % 3];
    if (!(//m.pmark(pa) == 15
        ((m.pmark(pa) == 14) && (m.pmark(pb) == 16)) ||
        ((m.pmark(pa) == 16) && (m.pmark(pb) == 14))
       )) {
      printf("  !! Wrong report\n");
      horror++;
    }
    // >>>>>>>>>>>>>>>>
    /*if (pos[1] != 3) {
      printf("  !! Wrong report\n");
      horror++;
    }*/
    pa = T2[pos[1]];
    pb = T2[(pos[1] + 1) % 3];
    if (!(m.pmark(pa) == 12
         //((m.pmark(pa) == 11) && (m.pmark(pb) == 12)) ||
         //((m.pmark(pa) == 12) && (m.pmark(pb) == 11))
       )) {
      printf("  !! Wrong report\n");
      horror++;
    }
    // <<<<<<<<<<<<< The second intersection
    if (types[1] != (int) tetgenmesh::DISJOINT) {
      printf("  !! Wrong case.\n");
      horror++;
    }
    /*if (pos[2] != 3) {
      printf("  !! Wrong report\n");
      horror++;
    }*/
    /*pa = T1[pos[2]];
    pb = T1[(pos[2] + 1) % 3];
    if (!(m.pmark(pa) == 15
        //((m.pmark(pa) == 15) && (m.pmark(pb) == 14)) ||
        //((m.pmark(pa) == 14) && (m.pmark(pb) == 15))
       )) {
      printf("  !! Wrong report\n");
      horror++;
    }*/
    // >>>>>>>>>>>>>
    /*if (pos[3] != 3) {
      printf("  !! Wrong report\n");
      horror++;
    }*/
    /*pa = T2[pos[3]];
    pb = T2[(pos[3] + 1) % 3];
    if (!(m.pmark(pa) == 12
        //((m.pmark(pa) == 12) && (m.pmark(pb) == 11)) ||
        //((m.pmark(pa) == 11) && (m.pmark(pb) == 12))
       )) {
      printf("  !! Wrong report\n");
      horror++;
    }*/
    // >>>>>>>>>>>>>>
    if (horror == 0) {
      printf("  test passed.\n");
    } else {
      //exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedralize()    The interface for users using TetGen library to       //
//                     generate tetrahedral meshes with all features.        //
//                                                                           //
// The sequence is roughly as follows.  Many of these steps can be skipped,  //
// depending on the command line switches.                                   //
//                                                                           //
// - Initialize constants and parse the command line.                        //
// - Read the vertices from a file and either                                //
//   - tetrahedralize them (no -r), or                                       //
//   - read an old mesh from files and reconstruct it (-r).                  //
// - Insert the PLC segments and facets (-p).                                //
// - Read the holes (-p), regional attributes (-pA), and regional volume     //
//   constraints (-pa).  Carve the holes and concavities, and spread the     //
//   regional attributes and volume constraints.                             //
// - Enforce the constraints on minimum quality bound (-q) and maximum       //
//   volume (-a). Also enforce the conforming Delaunay property (-q and -a). //
// - Promote the mesh's linear tetrahedra to higher order elements (-o).     //
// - Write the output files and print the statistics.                        //
// - Check the consistency and Delaunay property of the mesh (-C).           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetrahedralize(tetgenbehavior *b, tetgenio *in, tetgenio *out,
  tetgenio *addin, tetgenio *bgmin)
{
  tetgenmesh m;
  // Variables for timing the performance of TetGen (defined in time.h).
  clock_t tv[17];

  tv[0] = clock();

  m.b = b;
  m.in = in;
  exactinit();
  m.initializepools();
  m.transfernodes();

  tv[1] = clock();

  if (b->refine) {
    // m.reconstructmesh();
  } else {
    m.incrementaldelaunay();
  }

  tv[2] = clock();

  if (!b->quiet) {
    if (b->refine) {
      printf("Mesh reconstruction seconds:");
    } else {
      printf("Delaunay seconds:");
    }
    printf("  %g\n", (tv[2] - tv[1]) / (REAL) CLOCKS_PER_SEC);
  }

  if (b->plc) {
    m.meshsurface();
    m.delaunizesegments();
  }

  tv[3] = clock();

  if (!b->quiet) {
    if (b->plc) {
      if (b->diagnose != 1) {
        printf("Segment and facet ");
      } else {
        printf("Intersection "); 
      }
      printf("seconds:  %g\n", (tv[3] - tv[2]) / (REAL) CLOCKS_PER_SEC);
      if ((b->diagnose != 1) && (b->verbose > 0)) {
        printf("  Surface mesh seconds:  %g\n",
               (tv[3] - tv[2]) / (REAL) CLOCKS_PER_SEC);
        /*printf("  Segment recovery seconds:  %g\n",
               (tv[16] - tv[15]) / (REAL) CLOCKS_PER_SEC);
        printf("  facet recovery seconds:  %g\n",
               (tv[4] - tv[16]) / (REAL) CLOCKS_PER_SEC);*/
      }
    }
  }

  printf("\n");

  if (out != (tetgenio *) NULL) {
    out->firstnumber = in->firstnumber;
    out->mesh_dim = in->mesh_dim;
  }

  if (b->nonodewritten || b->noiterationnum) {
    if (!b->quiet) {
      printf("NOT writing a .node file.\n");
    }
  } else {
    m.outnodes(out);
  }

  if (b->noelewritten == 1) {
    if (!b->quiet) {
      printf("NOT writing an .ele file.\n");
    }
    // m.numberedges();
  } else {
    m.outelements(out);
  }

  if (b->nofacewritten) {
    if (!b->quiet) {
      printf("NOT writing an .face file.\n");
    }
  } else {
    if (b->facesout) {
      if (m.tetrahedronpool->items > 0l) {
        m.outfaces(out);  // Output all faces.
      }
    } else {
      if (b->plc || b->refine) {
        if (m.subfacepool->items > 0l) {
          m.outsubfaces(out); // Output boundary faces.
        }
      } else {
        if (m.tetrahedronpool->items > 0l) {
          m.outhullfaces(out); // Output convex hull faces.
        }
      }
    }
  }

  if (b->edgesout) {
    if (b->edgesout > 1) {
      m.outedges(out); // -ee, output all mesh edges. 
    } else {
      m.outsubsegments(out); // -e, only output subsegments.
    }
  }

  if (b->neighout) {
    m.outneighbors(out);
  }

  if (b->voroout) {
    // m.outvoronoi(out);
  }

  tv[4] = clock();

  if (!b->quiet) {
    printf("\nOutput seconds:  %g\n",
           (tv[4] - tv[3]) / (REAL) CLOCKS_PER_SEC);
    printf("Total running seconds:  %g\n",
           (tv[4] - tv[0]) / (REAL) CLOCKS_PER_SEC);
  }

  if (b->docheck) {
    m.checkdelaunay();
  }

  if (!b->quiet) {
    m.statistics();
  }
}

#ifndef TETLIBRARY

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// main()    The entrance for running TetGen from command line.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])

#else // with TETLIBRARY

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedralize()    The entrance for calling TetGen from another program. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetrahedralize(char *switches, tetgenio *in, tetgenio *out,
  tetgenio *addin, tetgenio *bgmin)

#endif // not TETLIBRARY

{
  tetgenbehavior b;

#ifndef TETLIBRARY

  tetgenio in, addin, bgmin;
  
  if (!b.parse_commandline(argc, argv)) {
    terminatetetgen(1);
  }
  if (b.refine) {
    if (!in.load_tetmesh(b.infilename)) {
      terminatetetgen(1);
    }
  } else {
    if (!in.load_plc(b.infilename, (int) b.object)) {
      terminatetetgen(1);
    }
  }
  if (b.insertaddpoints) {
    if (!addin.load_node(b.addinfilename)) {
      addin.numberofpoints = 0l;
    }
  }
  if (b.metric) {
    if (!bgmin.load_tetmesh(b.bgmeshfilename)) {
      bgmin.numberoftetrahedra = 0l;
    }
  }

  // FOR DEBUG -S1 option.
  if (b.steiner > 0) {
    test_tri_tri(&b, &in);
    terminatetetgen(0);
  }

  if (bgmin.numberoftetrahedra > 0l) {
    tetrahedralize(&b, &in, NULL, &addin, &bgmin);
  } else {
    tetrahedralize(&b, &in, NULL, &addin, NULL);
  }

  return 0;

#else // with TETLIBRARY

  if (!b.parse_commandline(switches)) {
    terminatetetgen(1);
  }
  tetrahedralize(&b, in, out, addin, bgmin);

#endif // not TETLIBRARY
}

#endif // #ifndef mainCXX