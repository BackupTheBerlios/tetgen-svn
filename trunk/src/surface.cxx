#ifndef surfaceCXX
#define surfaceCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// meshsurface()    Create a surface mesh of the input PLC.                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::meshsurface()
{
  list *ptlist, *conlist;
  memorypool *viri;
  point *idx2verlist;
  point tstart, tend, *cons;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  int *worklist;
  int end1, end2;
  int shmark, i, j;

  if (!b->quiet) {
    printf("Creating surface mesh.\n");
  }

  // Initialize working lists.
  ptlist = new list(sizeof(point *));
  conlist = new list(sizeof(point *) * 2);
  viri = new memorypool(sizeof(shellface *), 1024, POINTER, 0);
  worklist = new int[pointpool->items + 1];
  for (i = 0; i < pointpool->items + 1; i++) worklist[i] = 0;

  // Compute a mapping from indices to points.
  makeindex2pointmap(&idx2verlist);

  // Loop the facet list, triangulate each facet.
  for (shmark = 1; shmark <= in->numberoffacets; shmark++) {

    // Get a facet F.
    f = &in->facetlist[shmark - 1];

    // Process the duplicated points first, they are marked with type
    //   DUPLICATEDVERTEX.  If p and q are duplicated, and p'index > q's,
    //   then p is substituted by q.
    if (dupverts > 0l) {
      // Loop all polygons of this facet.
      for (i = 0; i < f->numberofpolygons; i++) {
        p = &(f->polygonlist[i]);
        // Loop other vertices of this polygon.
        for (j = 0; j < p->numberofvertices; j++) {
          end1 = p->vertexlist[j];
          tstart = idx2verlist[end1];
          if (pointtype(tstart) == DUPLICATEDVERTEX) {
            // Reset the index of vertex-j.
            tend = point2ppt(tstart);
            end2 = pointmark(tend);
            p->vertexlist[j] = end2;
          }
        }
      }
    }

    // Loop polygons of F, get the set of vertices and segments.
    for (i = 0; i < f->numberofpolygons; i++) {
      // Get a polygon.
      p = &(f->polygonlist[i]);
      // Get the first vertex.
      end1 = p->vertexlist[0];
      if ((end1 < in->firstnumber) || 
          (end1 >= in->firstnumber + in->numberofpoints)) {
        if (!b->quiet) {
          printf("Warning:  Invalid the 1st vertex %d of polygon", end1);
          printf(" %d in facet %d.\n", i + 1, shmark);
        }
        continue; // Skip this polygon.
      }
      tstart = idx2verlist[end1];
      // Add tstart to V if it haven't been added yet.
      if (worklist[end1] == 0) {
        ptlist->append(&tstart);
        worklist[end1] = 1;
      }
      // Loop other vertices of this polygon.
      for (j = 1; j <= p->numberofvertices; j++) {
        // get a vertex.
        if (j < p->numberofvertices) {
          end2 = p->vertexlist[j];
        } else {
          end2 = p->vertexlist[0];  // Form a loop from last to first.
        }
        if ((end2 < in->firstnumber) ||
            (end2 >= in->firstnumber + in->numberofpoints)) {
          if (!b->quiet) {
            printf("Warning:  Invalid vertex %d in polygon %d", end2, i + 1);
            printf(" in facet %d.\n", shmark);
          }
        } else {
          if (end1 != end2) {
            // 'end1' and 'end2' form a segment.
            tend = idx2verlist[end2];
            // Add tstart to V if it haven't been added yet.
            if (worklist[end2] == 0) {
              ptlist->append(&tend);
              worklist[end2] = 1;
            }
            // Save the segment in S (conlist).
            cons = (point *) conlist->append(NULL);
            cons[0] = tstart;
            cons[1] = tend;
            // Set the start for next continuous segment.
            end1 = end2;
            tstart = tend;
          } else {
            // Two identical vertices mean an isolated vertex of F.
            if (p->numberofvertices > 2) {
              // This may be an error in the input, anyway, we can continue
              //   by simply skipping this segment.
              if (!b->quiet) {
                printf("Warning:  Polygon %d has two identical verts", i + 1);
                printf(" in facet %d.\n", shmark);
              }
            } 
            // Ignore this vertex.
          }
        }
        // Is the polygon degenerate (a segment or a vertex)?
        if (p->numberofvertices == 2) break;
      }
    }
    // Unmark vertices.
    for (i = 0; i < ptlist->len(); i++) {
      tstart = * (point *) ptlist->get(i);
      end1 = pointmark(tstart);
      worklist[end1] = 0;
    }

    // Create a CDT of F.
// triangulate(shmark, ptlist, conlist, f->numberofholes, f->holelist, viri);

    // Clear working lists.
    ptlist->clear();
    conlist->clear();
    viri->restart();
  }

  // There are redundant segments in 'subsegpool', unify them, and build the
  //   face links of segments.
  // unifysegments();

  // Remember the number of input segments (for output).
  insegments = subsegpool->items;

  if (b->object == tetgenbehavior::STL) {
    // Remove redundant vertices (for .stl input mesh).
    jettisonnodes();
  }

  if (!b->nomerge && !b->nobisect) {
    // No '-M' switch - merge adjacent facets if they are coplanar.
    // mergefacets(flipqueue);
  }

  delete ptlist;
  delete conlist;
  delete viri;
  delete [] worklist;
  delete [] idx2verlist;
}

#endif // #ifndef surfaceCXX