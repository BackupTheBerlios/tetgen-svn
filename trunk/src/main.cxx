#ifndef mainCXX
#define mainCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Debug functions.                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::ptet(triface* t)
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
      printf("      [%i] x%lx  loc(%i)", facecount,
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

int tetgenmesh::pmark(point p)
{
  return pointmark(p);
}

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