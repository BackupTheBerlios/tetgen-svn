#ifndef mainCXX
#define mainCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Debug functions.                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::pmark(point p)
{
  return pointmark(p);
}

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
// terminatetetgen()    Terminate TetGen with a given exit code.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void terminatetetgen(int x)
{
#ifdef TETLIBRARY
  throw x;
#else
  exit(x);
#endif // #ifdef TETLIBRARY
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
        //if (m.subfaces->items > 0l) {
        //  m.outsubfaces(out); // Output boundary faces.
        //}
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
      // m.outsubsegments(out); // -e, only output subsegments.
    }
  }

  if (b->neighout) {
    m.outneighbors(out);
  }

  if (b->voroout) {
    // m.outvoronoi(out);
  }

  tv[3] = clock();

  if (!b->quiet) {
    printf("\nOutput seconds:  %g\n",
           (tv[3] - tv[2]) / (REAL) CLOCKS_PER_SEC);
    printf("Total running seconds:  %g\n",
           (tv[3] - tv[0]) / (REAL) CLOCKS_PER_SEC);
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