#ifndef constrainCXX
#define constrainCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertsegment()    Insert a given segment into tetrahedralization.        //
//                                                                           //
// Search an edge in current tetrahedralization that matches the given segm- //
// ment. Return TRUE if such edge exists, and the segment is attached to the //
// surrounding tets. Otherwise return FALSE.                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::insertsegment(face* insseg)
{
  triface searchtet;
  point pa, pb;
  enum locateresult loc;
  bool searchflag;
  int m, i;

  if (b->verbose > 1) {
    printf("    Insert segment (%d, %d).\n", pointmark(insseg->sh[3]), 
      pointmark(insseg->sh[4]));
  }

  searchflag = false;

  // Search a tet containing an endpoint of this segment as origin.
  for (m = 0; m < 2 && !searchflag; m++) {
    pa = (point) insseg.sh[m + 3];
    decode(point2tet(pa), searchtet);
    if (searchtet.tet[4] != NULL) {
      // Check if this tet contains pa.
      for (i = 4; i < 8 && !searchflag; i++) {
        if (searchtet.tet[i] == pa) {
          searchflag = true;
          // Found. Set pa as its origin.
          switch (i) {
          case 4: searchtet.loc = 0; searchtet.ver = 0; break;
          case 5: searchtet.loc = 0; searchtet.ver = 2; break;
          case 6: searchtet.loc = 0; searchtet.ver = 4; break;
          case 7: searchtet.loc = 1; searchtet.ver = 2; break;
          }
          pb = (point) insseg.sh[(m + 1) % 2 + 3];
        }
      }
    }
  }

  if (!searchflag) {
    // Locate pa in tetrahedralization.
    pa = (point) insseg.sh[3];
    randomsample(pa, &searchtet);
    loc = locate(pa, &searchtet);
    assert(loc == ONVERTEX);  // SELF_CHECK
    pb = (point) insseg.sh[4];
    force_ptloc_count++;
  }

  if (b->verbose > 1) {
    printf("    Search from (%d, %d, %d, %d) to %d.\n", pointmark(pa), 
      pointmark(dest(searchtet)), pointmark(apex(searchtet)),
      pointmark(oppo(searchtet)), pointmark(pb));
  }

}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// delaunizesegments()    Recover segments in a Delaunay tetrahedralization. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::delaunizesegments()
{


  if (!b->quiet) {
    printf("Delaunizing segments.\n");
  }

  // Segments will be attached to tets. Initialize the connection pool.
  tet2subpool = new memorypool(6*sizeof(shellface), SUBPERBLOCK, POINTER, 0);



  delete tet2subpool;
}

#endif // #ifndef constrainCXX
