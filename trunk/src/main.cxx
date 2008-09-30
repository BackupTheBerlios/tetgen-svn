#ifndef mainCXX
#define mainCXX

#include "tetgen.h"

static REAL PI = 3.14159265358979323846264338327950288419716939937510582;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetgenbehavior()    Initialize veriables of 'tetgenbehavior'.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenbehavior::tetgenbehavior()
{
  // Initialize command line switches.
  plc = 0;
  quality = 0;
  refine = 0;
  coarse = 0;
  metric = 0;
  minratio = 2.0;
  goodratio = 0.0;
  minangle = 20.0;
  goodangle = 0.0;
  maxdihedral = 165.0;
  mindihedral = 5.0;
  varvolume = 0;
  fixedvolume = 0;
  maxvolume = -1.0;
  regionattrib = 0;
  bowyerwatson = 0;
  insertaddpoints = 0;
  diagnose = 0;
  conformdel = 0;
  zeroindex = 0;
  facesout = 0;
  edgesout = 0;
  neighout = 0;
  voroout = 0;
  meditview = 0;
  gidview = 0;
  geomview = 0;
  order = 1;
  nojettison = 0;
  nobound = 0;
  nonodewritten = 0;
  noelewritten = 0;
  nofacewritten = 0;
  noiterationnum = 0;
  nobisect = 0;
  steiner = -1;
  nomerge = 0;
  docheck = 0;
  quiet = 0;
  verbose = 0;
  useshelles = 0;
  epsilon = 1.0e-8;
  object = NONE;
  // Initialize strings
  commandline[0] = '\0';
  infilename[0] = '\0';
  outfilename[0] = '\0';
  addinfilename[0] = '\0';
  bgmeshfilename[0] = '\0';
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// versioninfo()    Print the version information of TetGen.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenbehavior::versioninfo()
{
  printf("Develop Version (Started on August 9, 2008).\n");
  printf("\n");
  printf("Copyright (C) 2002 - 2008\n");
  printf("Hang Si\n");
  printf("Mohrenstr. 39, 10117 Berlin, Germany\n");
  printf("si@wias-berlin.de\n");
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// syntax()    Print list of command line switches and exit the program.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenbehavior::syntax()
{
  printf("  tetgen [-prq_Ra_AiMYS_T_dzo_fenvgGOJBNEFICQVh] input_file\n");
  printf("    -p  Tetrahedralizes a piecewise linear complex (PLC).\n");
  printf("    -r  Reconstructs a previously generated mesh.\n");
  printf("    -q  Quality mesh generation (adding new mesh points to ");
  printf("improve mesh quality).\n");
  printf("    -R  Mesh coarsening (deleting redundant mesh points).\n");
  printf("    -a  Applies a maximum tetrahedron volume constraint.\n");
  printf("    -A  Assigns attributes to identify tetrahedra in different ");
  printf("regions.\n");
  printf("    -i  Inserts a list of additional points into mesh.\n");
  printf("    -M  Does not merge coplanar facets.\n");
  printf("    -Y  Suppresses boundary facets/segments splitting.\n");
  printf("    -S  Specifies maximum number of added points.\n");
  printf("    -T  Sets a tolerance for coplanar test (default 1e-8).\n");
  printf("    -d  Detects self-intersections of facets of the PLC.\n");
  printf("    -z  Numbers all output items starting from zero.\n");
  printf("    -o2 Generates second-order subparametric elements.\n");
  printf("    -f  Outputs all faces to .face file.");
  printf("file.\n");
  printf("    -e  Outputs all edges to .edge file.\n");
  printf("    -n  Outputs tetrahedra neighbors to .neigh file.\n");
  printf("    -v  Outputs Voronoi diagram to files.\n");
  printf("    -g  Outputs mesh to .mesh file for viewing by Medit.\n");
  printf("    -G  Outputs mesh to .msh file for viewing by Gid.\n");
  printf("    -O  Outputs mesh to .off file for viewing by Geomview.\n");
  printf("    -J  No jettison of unused vertices from output .node file.\n");
  printf("    -B  Suppresses output of boundary information.\n");
  printf("    -N  Suppresses output of .node file.\n");
  printf("    -E  Suppresses output of .ele file.\n");
  printf("    -F  Suppresses output of .face file.\n");
  printf("    -I  Suppresses mesh iteration numbers.\n");
  printf("    -C  Checks the consistency of the final mesh.\n");
  printf("    -Q  Quiet:  No terminal output except errors.\n");
  printf("    -V  Verbose:  Detailed information, more terminal output.\n");
  printf("    -h  Help:  A brief instruction for using TetGen.\n");
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// usage()    Print a brief instruction for using TetGen.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenbehavior::usage()
{
  printf("TetGen\n");
  printf("A Quality Tetrahedral Mesh Generator and 3D Delaunay ");
  printf("Triangulator\n");
  versioninfo();
  printf("\n");
  printf("What Can TetGen Do?\n");
  printf("\n");
  printf("  TetGen generates exact Delaunay tetrahedralizations, exact\n");
  printf("  constrained Delaunay tetrahedralizations, and quality ");
  printf("tetrahedral\n  meshes. The latter are nicely graded and whose ");
  printf("tetrahedra have\n  radius-edge ratio bounded, thus are suitable ");
  printf("for finite element and\n  finite volume analysis.\n"); 
  printf("\n");
  printf("Command Line Syntax:\n");
  printf("\n");
  printf("  Below is the command line syntax of TetGen with a list of ");
  printf("short\n");
  printf("  descriptions. Underscores indicate that numbers may optionally\n");
  printf("  follow certain switches.  Do not leave any space between a ");
  printf("switch\n");
  printf("  and its numeric parameter.  \'input_file\' contains input data\n");
  printf("  depending on the switches you supplied which may be a ");
  printf("  piecewise\n");
  printf("  linear complex or a list of nodes.  File formats and detailed\n");
  printf("  description of command line switches are found in user's ");
  printf("manual.\n");
  printf("\n");
  syntax();
  printf("\n");
  printf("Examples of How to Use TetGen:\n");
  printf("\n");
  printf("  \'tetgen object\' reads vertices from object.node, and writes ");
  printf("their\n  Delaunay tetrahedralization to object.1.node and ");
  printf("object.1.ele.\n");
  printf("\n");
  printf("  \'tetgen -p object\' reads a PLC from object.poly or object.");
  printf("smesh (and\n  possibly object.node) and writes its constrained ");
  printf("Delaunay\n  tetrahedralization to object.1.node, object.1.ele and ");
  printf("object.1.face.\n");
  printf("\n");
  printf("  \'tetgen -pq1.414a.1 object\' reads a PLC from object.poly or\n");
  printf("  object.smesh (and possibly object.node), generates a mesh ");
  printf("whose\n  tetrahedra have radius-edge ratio smaller than 1.414 and ");
  printf("have volume\n  of 0.1 or less, and writes the mesh to ");
  printf("object.1.node, object.1.ele\n  and object.1.face.\n");
  printf("\n");
  printf("Please send bugs/comments to Hang Si <si@wias-berlin.de>\n");
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// parse_commandline()    Read the command line, identify switches, and set  //
//                        up options and file names.                         //
//                                                                           //
// 'argc' and 'argv' are the same parameters passed to the function main()   //
// of a C/C++ program. They together represent the command line user invoked //
// from an environment in which TetGen is running.                           //
//                                                                           //
// When TetGen is invoked from an environment. 'argc' is nonzero, switches   //
// and input filename should be supplied as zero-terminated strings in       //
// argv[0] through argv[argc - 1] and argv[0] shall be the name used to      //
// invoke TetGen, i.e. "tetgen".  Switches are previously started with a     //
// dash '-' to identify them from the input filename.                        //
//                                                                           //
// When TetGen is called from within another program. 'argc' is set to zero. //
// switches are given in one zero-terminated string (no previous dash is     //
// required.), and 'argv' is a pointer points to this string.  No input      //
// filename is required (usually the input data has been directly created by //
// user in the 'tetgenio' structure).  A default filename 'tetgen-tmpfile'   //
// will be created for debugging output purpose.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenbehavior::parse_commandline(int argc, char **argv)
{
  int startindex;
  int increment;
  int meshnumber;
  int scount;
  int i, j, k;
  char workstring[1024];

  // First determine the input style of the switches.
  if (argc == 0) {
    startindex = 0;                    // Switches are given without a dash.
    argc = 1;                    // For running the following for-loop once.
    commandline[0] = '\0';
  } else {
    startindex = 1;
    strcpy(commandline, argv[0]);
    strcat(commandline, " ");
  }
  
  // Rcount used to count the number of '-R' be used.
  scount = 0;

  for (i = startindex; i < argc; i++) {
    // Remember the command line switches.
    strcat(commandline, argv[i]);
    strcat(commandline, " ");
    if (startindex == 1) {
      // Is this string a filename?
      if (argv[i][0] != '-') {
        strncpy(infilename, argv[i], 1024 - 1);
        infilename[1024 - 1] = '\0';
        // Go to the next string directly.
        continue;                     
      }
    }
    // Parse the individual switch from the string.
    for (j = startindex; argv[i][j] != '\0'; j++) {
      if (argv[i][j] == 'p') {
        plc = 1;
      } else if (argv[i][j] == 'r') {
        refine = 1;
      } else if (argv[i][j] == 'R') {
        coarse = 1;
      } else if (argv[i][j] == 'q') {
        quality++;
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          if (quality == 1) {
            minratio = (REAL) strtod(workstring, (char **) NULL);
          } else if (quality == 2) {
            mindihedral = (REAL) strtod(workstring, (char **) NULL);
          } else if (quality == 3) {
            maxdihedral = (REAL) strtod(workstring, (char **) NULL);
          }
        }
      } else if (argv[i][j] == 'm') {
        metric++;
      } else if (argv[i][j] == 'a') {
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          fixedvolume = 1;
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                 (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          maxvolume = (REAL) strtod(workstring, (char **) NULL);
        } else {
          varvolume = 1;
        }
      } else if (argv[i][j] == 'A') {
        regionattrib++;
      } else if (argv[i][j] == 'b') {
        bowyerwatson = 1;
      } else if (argv[i][j] == 'i') {
        insertaddpoints = 1;
      } else if (argv[i][j] == 'd') {
        diagnose = 1;
      } else if (argv[i][j] == 'z') {
        zeroindex = 1;
      } else if (argv[i][j] == 'f') {
        facesout = 1;
      } else if (argv[i][j] == 'e') {
        edgesout++;
      } else if (argv[i][j] == 'n') {
        neighout++;
      } else if (argv[i][j] == 'v') {
        voroout = 1;
      } else if (argv[i][j] == 'g') {
        meditview = 1;
      } else if (argv[i][j] == 'G') {
        gidview = 1;
      } else if (argv[i][j] == 'O') {
        geomview = 1;
      } else if (argv[i][j] == 'M') {
        nomerge = 1;
      } else if (argv[i][j] == 'Y') {
        nobisect++;
      } else if (argv[i][j] == 'J') {
        nojettison = 1;
      } else if (argv[i][j] == 'B') {
        nobound = 1;
      } else if (argv[i][j] == 'N') {
        nonodewritten = 1;
      } else if (argv[i][j] == 'E') {
        noelewritten = 1;
        if (argv[i][j + 1] == '2') {
          j++;
          noelewritten = 2;
        }
      } else if (argv[i][j] == 'F') {
        nofacewritten = 1;
      } else if (argv[i][j] == 'I') {
        noiterationnum = 1;
      } else if (argv[i][j] == 'o') {
        if (argv[i][j + 1] == '2') {
          j++;
          order = 2;
        }
      } else if (argv[i][j] == 'S') {
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                 (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          steiner = (int) strtol(workstring, (char **) NULL, 0);
        } 
      } else if (argv[i][j] == 'D') {
        conformdel++;
      } else if (argv[i][j] == 'T') {
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                 (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          epsilon = (REAL) strtod(workstring, (char **) NULL);
        } 
      } else if (argv[i][j] == 'C') {
        docheck++;
      } else if (argv[i][j] == 'Q') {
        quiet = 1;
      } else if (argv[i][j] == 'V') {
        verbose++;
      // } else if (argv[i][j] == 'v') {
        // versioninfo();
        // terminatetetgen(0);
      } else if ((argv[i][j] == 'h') || (argv[i][j] == 'H') ||
                 (argv[i][j] == '?')) {
        usage();
        terminatetetgen(0);
      } else {
        printf("Warning:  Unknown switch -%c.\n", argv[i][j]);
      }
    }
  }

  if (startindex == 0) {
    // Set a temporary filename for debugging output.
    strcpy(infilename, "tetgen-tmpfile");
  } else {
    if (infilename[0] == '\0') {
      // No input file name. Print the syntax and exit.
      syntax();
      terminatetetgen(0);
    }
    // Recognize the object from file extension if it is available.
    if (!strcmp(&infilename[strlen(infilename) - 5], ".node")) {
      infilename[strlen(infilename) - 5] = '\0';
      object = NODES;
    } else if (!strcmp(&infilename[strlen(infilename) - 5], ".poly")) {
      infilename[strlen(infilename) - 5] = '\0';
      object = POLY;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 6], ".smesh")) {
      infilename[strlen(infilename) - 6] = '\0';
      object = POLY;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".off")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = OFF;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".ply")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = PLY;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".stl")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = STL;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 5], ".mesh")) {
      infilename[strlen(infilename) - 5] = '\0';
      object = MEDIT;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".vtk")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = VTK;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".ele")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = MESH;
      refine = 1;
    }
  }
  plc = plc || diagnose;
  useshelles = plc || refine || coarse || quality;
  goodratio = minratio;
  goodratio *= goodratio;

  // Detect improper combinations of switches.
  if (plc && refine) {
    printf("Error:  Switch -r cannot use together with -p.\n");
    return false;
  }
  if (refine && (plc || noiterationnum)) {
    printf("Error:  Switches %s cannot use together with -r.\n",
           "-p, -d, and -I");
    return false;
  }
  if (diagnose && (quality || insertaddpoints || (order == 2) || neighout
      || docheck)) {
    printf("Error:  Switches %s cannot use together with -d.\n",
           "-q, -i, -o2, -n, and -C");
    return false;
  }

  // Be careful not to allocate space for element area constraints that 
  //   will never be assigned any value (other than the default -1.0).
  if (!refine && !plc) {
    varvolume = 0;
  }
  // Be careful not to add an extra attribute to each element unless the
  //   input supports it (PLC in, but not refining a preexisting mesh).
  if (refine || !plc) {
    regionattrib = 0;
  }
  // If '-a' or '-aa' is in use, enable '-q' option too.
  if (fixedvolume || varvolume) {
    if (quality == 0) {
      quality = 1;
    }
  }
  // Calculate the goodangle for testing bad subfaces.
  goodangle = cos(minangle * PI / 180.0);
  goodangle *= goodangle;

  increment = 0;
  strcpy(workstring, infilename);
  j = 1;
  while (workstring[j] != '\0') {
    if ((workstring[j] == '.') && (workstring[j + 1] != '\0')) {
      increment = j + 1;
    }
    j++;
  }
  meshnumber = 0;
  if (increment > 0) {
    j = increment;
    do {
      if ((workstring[j] >= '0') && (workstring[j] <= '9')) {
        meshnumber = meshnumber * 10 + (int) (workstring[j] - '0');
      } else {
        increment = 0;
      }
      j++;
    } while (workstring[j] != '\0');
  }
  if (noiterationnum) {
    strcpy(outfilename, infilename);
  } else if (increment == 0) {
    strcpy(outfilename, infilename);
    strcat(outfilename, ".1");
  } else {
    workstring[increment] = '%';
    workstring[increment + 1] = 'd';
    workstring[increment + 2] = '\0';
    sprintf(outfilename, workstring, meshnumber + 1);
  }
  // Additional input file name has the end ".a".
  strcpy(addinfilename, infilename);
  strcat(addinfilename, ".a");
  // Background filename has the form "*.b.ele", "*.b.node", ...
  strcpy(bgmeshfilename, infilename);
  strcat(bgmeshfilename, ".b");

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Debug functions.                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Print the detail informations of a tetrahedron.

void tetgenmesh::ptet(triface* t)
{
  triface tmpface, prtface;
  point *pts, tmppt;
  REAL ori;
  int facecount;

  printf("Tetra x%lx with loc(%i) ver(%i):",
         (unsigned long)(t->tet), t->loc, t->ver);
  if (t->tet == NULL) {
    printf("  !! NOT A VALID HANDLE\n");
    return;
  }
  pts = (point *) t->tet;
  if (pts[7] != dummypoint) {
    ori = orient3d(pts[4], pts[5], pts[6], pts[7]);
    printf("  ori = %g.\n", ori);
  } else {
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
      printf("      [%i] x%lx  loc(%i) ver(%i)", facecount,
             (unsigned long)(prtface.tet), prtface.loc, prtface.ver);
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
// Print the detail informations of a shellface.

void tetgenmesh::psh(face *s)
{
  face prtsh;
  triface prttet;
  point *pt, printpoint;
  REAL n[3];

  if (s->sh == NULL) {
    printf("Not a handle.\n");
    return;
  }

  if (s->sh[3] == NULL) {
    printf("A dead subface x%lx.\n", (unsigned long)(s->sh));
    return;
  }

  pt = (point *) s->sh;
  if (s->sh[5] != NULL) {
    printf("subface x%lx, ver %d, mark %d,", (unsigned long)(s->sh), s->shver,
      shellmark(*s));
    facenormal(pt[3], pt[4], pt[5], n, 1);
    printf(" area %g, edge lengths %g %g %g", sqrt(DOT(n, n)),
      DIST(pt[3], pt[4]), DIST(pt[4], pt[5]), DIST(pt[5], pt[3]));
  } else {
    printf("Subsegment x%lx, ver %d, mark %d,", (unsigned long)(s->sh),
      s->shver, shellmark(*s));
    printf(" length %g:", DIST(pt[3], pt[4]));
  }
  // if (sinfected(*sface)) {
  //   printf(" (infected)");
  // }
  // if (shell2badface(*sface)) {
  //   printf(" (queued)");
  // }
  // if (checkpbcs) {
  //   if (shellpbcgroup(*sface) >= 0) {
  //     printf(" (pbc %d)", shellpbcgroup(*sface));
  //   }
  // }
  printf("\n");

  sdecode(s->sh[0], prtsh);
  if (prtsh.sh == NULL) {
    printf("      [0] = No shell\n");
  } else {
    printf("      [0] = x%lx  %d\n", (unsigned long)(prtsh.sh), prtsh.shver);
  }
  sdecode(s->sh[1], prtsh);
  if (prtsh.sh == NULL) {
    printf("      [1] = No shell\n");
  } else {
    printf("      [1] = x%lx  %d\n", (unsigned long)(prtsh.sh), prtsh.shver);
  }
  sdecode(s->sh[2], prtsh);
  if (prtsh.sh == NULL) {
    printf("      [2] = No shell\n");
  } else {
    printf("      [2] = x%lx  %d\n", (unsigned long)(prtsh.sh), prtsh.shver);
  }

  printpoint = sorg(*s);
  if (printpoint == (point) NULL)
    printf("      Org [%d] = NULL\n", vo[s->shver]);
  else
    printf("      Org [%d] = x%lx  (%.12g,%.12g,%.12g) %d\n",
           vo[s->shver], (unsigned long)(printpoint), printpoint[0],
           printpoint[1], printpoint[2], pointmark(printpoint));
  printpoint = sdest(*s);
  if (printpoint == (point) NULL)
    printf("      Dest[%d] = NULL\n", vd[s->shver]);
  else
    printf("      Dest[%d] = x%lx  (%.12g,%.12g,%.12g) %d\n",
            vd[s->shver], (unsigned long)(printpoint), printpoint[0],
            printpoint[1], printpoint[2], pointmark(printpoint));

  if (sapex(*s) != NULL) {
    printpoint = sapex(*s);
    if (printpoint == (point) NULL)
      printf("      Apex[%d] = NULL\n", va[s->shver]);
    else
      printf("      Apex[%d] = x%lx  (%.12g,%.12g,%.12g) %d\n",
             va[s->shver], (unsigned long)(printpoint), printpoint[0],
             printpoint[1], printpoint[2], pointmark(printpoint));

    sdecode(s->sh[6], prtsh);
    if (prtsh.sh == NULL) {
      printf("      [6] = No subsegment\n");
    } else {
      printf("      [6] = x%lx  %d\n",
             (unsigned long)(prtsh.sh), prtsh.shver);
    }
    sdecode(s->sh[7], prtsh);
    if (prtsh.sh == NULL) {
      printf("      [7] = No subsegment\n");
    } else {
      printf("      [7] = x%lx  %d\n",
             (unsigned long)(prtsh.sh), prtsh.shver);
    }
    sdecode(s->sh[8], prtsh);
    if (prtsh.sh == NULL) {
      printf("      [8] = No subsegment\n");
    } else {
      printf("      [8] = x%lx  %d\n",
             (unsigned long)(prtsh.sh), prtsh.shver);
    }

    decode(s->sh[9], prttet);
    if (prttet.tet == NULL) {
      printf("      [9] = Outer space\n");
    } else {
      printf("      [9] = x%lx  %d\n",
             (unsigned long)(prttet.tet), prttet.loc);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Find and print the tetrahedron (or face or edge) with the given indices.
// Do not handle 'dummypoint' (-1).

void tetgenmesh::ptet(int i, int j, int k, int l)
{
  triface t;
  point *pts;
  int *marklist;
  int ii;

  marklist = new int[pointpool->items + 1];
  for (ii = 0; ii < pointpool->items + 1; ii++) marklist[ii] = 0;
  // Marke the given indices. 
  marklist[i] = marklist[j] = marklist[k] = marklist[l] = 1;

  t.loc = t.ver = 0;
  tetrahedronpool->traversalinit();
  t.tet = tetrahedrontraverse(tetrahedronpool);
  while (t.tet != NULL) {
    pts = (point *) t.tet;
    if (pts[7] != dummypoint) {
      if ((marklist[pointmark(pts[4])] + marklist[pointmark(pts[5])] +
           marklist[pointmark(pts[6])] + marklist[pointmark(pts[7])]) == 4) {
        ptet(&t);  // Find!
        break;
      }
    }
    t.tet = tetrahedrontraverse(tetrahedronpool);
  }
  
  if (t.tet == NULL) {
    printf("  !! Not exist.\n");
  }
  delete [] marklist;
}

void tetgenmesh::pface(int i, int j, int k)
{
  triface t, t1;
  point *pts;
  REAL sign;
  int *marklist;
  int ii;

  marklist = new int[pointpool->items + 1];
  for (ii = 0; ii < pointpool->items + 1; ii++) marklist[ii] = 0;
  // Marke the given indices. 
  marklist[i] = marklist[j] = marklist[k] = 1;

  t.ver = t1.ver = 0;
  tetrahedronpool->traversalinit();
  t.tet = tetrahedrontraverse(tetrahedronpool);
  while (t.tet != NULL) {
    pts = (point *) t.tet;
    if (pts[7] != dummypoint) {
      if ((marklist[pointmark(pts[4])] + marklist[pointmark(pts[5])] +
           marklist[pointmark(pts[6])] + marklist[pointmark(pts[7])]) == 3) {
        // Find a tet containing the search face.
        for (t.loc = 0; t.loc < 4; t.loc++) {
          sym(t, t1);
          pts = (point *) t1.tet;
          if ((marklist[pointmark(pts[4])] + marklist[pointmark(pts[5])] +
              marklist[pointmark(pts[6])] + 
              (pts[7] != dummypoint ? marklist[pointmark(pts[7])] : 0)) == 3)
            break;
        }
        assert(t.loc < 4);
        // Now t and t1 share the face.
        printf("  tet x%lx (%d, %d, %d, %d) %d\n", (unsigned long) t.tet,
          pointmark(org(t)), pointmark(dest(t)), pointmark(apex(t)), 
          pointmark(oppo(t)), t.loc);
        printf("  tet x%lx (%d, %d, %d, %d) %d\n", (unsigned long) t1.tet,
          pointmark(org(t1)), pointmark(dest(t1)), pointmark(apex(t1)),
          pointmark(oppo(t1)), t1.loc);
        if ((point) t1.tet[7] != dummypoint) {
          pts = (point *) t.tet;
          sign = insphere(pts[4], pts[5], pts[6], pts[7], oppo(t1));
          printf("  %s (sign = %.g).\n", sign > 0 ? "Delaunay" :
            (sign < 0 ? "Non-Delaunay" : "Cosphere"), sign);
          if (sign == 0) {
            sign = insphere_sos(pts[4], pts[5], pts[6], pts[7], oppo(t1));
            printf("  %s (symbolic).\n", sign > 0 ? "Delaunay":"Non-Delaunay");
          }
        }
        break;
      }
    }
    t.tet = tetrahedrontraverse(tetrahedronpool);
  }
  
  if (t.tet == NULL) {
    printf("  !! Not exist.\n");
  }
  delete [] marklist;
}

///////////////////////////////////////////////////////////////////////////////
// Print the index of the point.

int tetgenmesh::pmark(point p)
{
  return pointmark(p);
}

///////////////////////////////////////////////////////////////////////////////
// Geometrical tests.

REAL tetgenmesh::test_orient3d(int i, int j, int k, int l)
{
  point *idx2ptmap;
  REAL ori;
  int idx;

  idx = (int) pointpool->items;
  if ((i > idx) || (j > idx) || (k > idx) || (l > idx)) {
    printf("Input indices are invalid.\n");
    return 0;
  }

  makeindex2pointmap(idx2ptmap);
  ori = orient3d(idx2ptmap[i], idx2ptmap[j], idx2ptmap[k], idx2ptmap[l]);
  delete [] idx2ptmap;
  
  return ori;
}

REAL tetgenmesh::test_insphere(int i, int j, int k, int l, int m)
{
  point *idx2ptmap;
  REAL sign;
  int idx;

  idx = (int) pointpool->items;
  if ((i > idx) || (j > idx) || (k > idx) || (l > idx) || (m > idx)) {
    printf("Input indices are invalid.\n");
    return 0;
  }

  makeindex2pointmap(idx2ptmap);
  sign = insphere(idx2ptmap[i], idx2ptmap[j], idx2ptmap[k], idx2ptmap[l],
                  idx2ptmap[m]);
  if (sign == 0) {
    printf("  sign == 0.0! (symbolic perturbed) \n");
    sign = insphere_sos(idx2ptmap[i], idx2ptmap[j], idx2ptmap[k], idx2ptmap[l],
                        idx2ptmap[m]);
  }
  delete [] idx2ptmap;
  
  return sign;
}

///////////////////////////////////////////////////////////////////////////////
// Print an array of tetrahedra (in draw command)

void tetgenmesh::print_tetarray(int n, triface *tetarray)
{
  int i;

  for (i = 0; i < n; i++) {
    printf("p:draw_tet(%d, %d, %d, %d) -- %d\n", 
      pointmark(org(tetarray[i])), pointmark(dest(tetarray[i])),
      pointmark(apex(tetarray[i])), pointmark(oppo(tetarray[i])), i + 1);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Print a list of faces (in draw command)

void tetgenmesh::print_tetfacelist(list *tetfacelist)
{
  triface *f;
  int i;

  for (i = 0; i < tetfacelist->len(); i++) {
    f = (triface *) tetfacelist->get(i);
    if (f->tet[4] != NULL) {
      printf("p:draw_subface(%d, %d, %d) -- op (%d), x%lx (%d)\n", 
        pointmark(org(*f)), pointmark(dest(*f)), pointmark(apex(*f)), 
        pointmark(oppo(*f)), (unsigned long) f->tet, i+1);
    } else {
      printf("!! dead face x%lx (i)\n", (unsigned long) f->tet, i+1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Print current faces in flipstack (in draw command)

void tetgenmesh::print_flipstack()
{
  badface *traveface;
  int i;

  traveface = futureflip; 
  i = 0;
  while (traveface != NULL) {
    if (traveface->tt.tet[4] != NULL) {
      printf("%2d  (%d, %d, %d, %d) - %d\n",i+1,pointmark(org(traveface->tt)),
        pointmark(dest(traveface->tt)), pointmark(apex(traveface->tt)),
        pointmark(oppo(traveface->tt)), pointmark(traveface->foppo));
    }
    traveface = traveface->nextitem;
    i++;
  }
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

  if (b->plc) {
    m.meshsurface();
    tv[3] = clock();
  }

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
    m.numberedges();
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