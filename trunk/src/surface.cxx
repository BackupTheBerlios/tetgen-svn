#ifndef surfaceCXX
#define surfaceCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// triangulate()    Create a CDT for the facet.                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::triangulate(int shmark, arraypool* ptlist, arraypool* conlist,
  int holes, REAL* holelist)
{
  face newsh, newseg;
  point *ppa, *ppb, *ppc;
  int i;

  if (b->verbose > 1) {
    printf("    %ld vertices, %ld segments",ptlist->objects,conlist->objects);
    if (holes > 0) {
      printf(", %d holes", holes);
    }
    printf(", shmark: %d.\n", shmark);
  }

  if ((ptlist->objects == 3) && (conlist->objects == 3)) {
    // This CDT contains only one triangle.
    ppa = (point *) fastlookup(ptlist, 0);
    ppb = (point *) fastlookup(ptlist, 1);
    ppc = (point *) fastlookup(ptlist, 2);
    makeshellface(subfacepool, &newsh);
    setshvertices(newsh, *ppa, *ppb, *ppc);
    shellmark(newsh) = shmark;
    // Create three new segments.
    for (i = 0; i < 3; i++) {
      makeshellface(subsegpool, &newseg);
      setshvertices(newseg, sorg(newsh), sdest(newsh), NULL);
      ssbond(newsh, newseg);
      senextself(newsh);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// unifysegments()    Remove redundant segments and create face links.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::unifysegments()
{
  badface *facelink, *newlinkitem, *f1, *f2;
  face *facperverlist, sface;
  face subsegloop, testseg;
  point torg, tdest;
  REAL ori1, ori2, ori3;
  REAL n1[3], n2[3];
  int *idx2faclist;
  int segmarker;
  int idx, k, m;

  if (b->verbose) {
    printf("  Unifying segments.\n");
  }

  // Create a mapping from vertices to subfaces incident at them.
  makepoint2submap(subfacepool, idx2faclist, facperverlist);

  segmarker = 1;
  subsegloop.shver = 0;
  subsegpool->traversalinit();
  subsegloop.sh = shellfacetraverse(subsegpool);
  while (subsegloop.sh != (shellface *) NULL) {
    torg = sorg(subsegloop);
    tdest = sdest(subsegloop);

    idx = pointmark(torg) - in->firstnumber;
    // Loop through the set of subfaces containing 'torg'.  Get all the
    //   subfaces containing the edge (torg, tdest). Save and order them
    //   in 'sfacelist', the ordering is defined by the right-hand rule
    //   with thumb points from torg to tdest.
    for (k = idx2faclist[idx]; k < idx2faclist[idx + 1]; k++) {
      sface = facperverlist[k];
      // The face may be deleted if it is a duplicated face.
      if (sface.sh[3] == NULL) continue;
      // Search the edge torg->tdest.
      assert(sorg(sface) == torg); // SELF_CHECK
      if (sdest(sface) != tdest) {
        senext2self(sface);
        sesymself(sface);
      }
      if (sdest(sface) != tdest) continue;

      // Save the face f in facelink.
      if (flippool->items >= 2) {
        f1 = facelink;
        for (m = 0; m < flippool->items - 1; m++) {
          f2 = f1->nextitem;
          ori1 = orient3d(torg, tdest, sapex(f1->ss), sapex(f2->ss));
          ori2 = orient3d(torg, tdest, sapex(f1->ss), sapex(sface));
          if (ori1 > 0) {
            // apex(f2) is below f1.
            if (ori2 > 0) {
              // apex(f) is below f1 (see Fig.1). 
              ori3 = orient3d(torg, tdest, sapex(f2->ss), sapex(sface));
              if (ori3 > 0) {
                // apex(f) is below f2, insert it.
                break; 
              } else if (ori3 < 0) {
                // apex(f) is above f2, continue.
              } else { // ori3 == 0; 
                // f is coplanar and codirection with f2. 
                assert(0);
              }
            } else if (ori2 < 0) {
              // apex(f) is above f1 below f2, inset it (see Fig. 2).
              break;
            } else { // ori2 == 0;
              // apex(f) is coplanar with f1 (see Fig. 5).
              ori3 = orient3d(torg, tdest, sapex(f2->ss), sapex(sface));
              if (ori3 > 0) {
                // apex(f) is below f2, insert it.
                break; 
              } else {
                // f is coplanar and codirection with f1.
                assert(0);
              }
            }
          } else if (ori1 < 0) {
            // apex(f2) is above f1.
            if (ori2 > 0) {
              // apex(f) is below f1, continue (see Fig. 3).
            } else if (ori2 < 0) {
              // apex(f) is above f1 (see Fig.4).
              ori3 = orient3d(torg, tdest, sapex(f2->ss), sapex(sface));
              if (ori3 > 0) {
                // apex(f) is below f2, insert it.
                break;
              } else if (ori3 < 0) {
                // apex(f) is above f2, continue.
              } else { // ori3 == 0;
                // f is coplanar and codirection with f2.
                assert(0);
              }
            } else { // ori2 == 0;
              // f is coplanar and codirection with f1.
              assert(0);
            }
          } else { // ori1 == 0;
            // apex(f2) is coplanar with f1. By assumption, f1 is not
            //   coplanar and codirection with f2.
            if (ori2 > 0) {
              // apex(f) is below f1, continue (see Fig. 7).
            } else if (ori2 < 0) {
              // apex(f) is above f1, insert it (see Fig. 7).
              break;
            } else { // ori2 == 0.
              // apex(f) is coplanar with f1 (see Fig. 8).
              assert(0);  // Two more cases need to work out here.
            }
          }
          // Go to the next item;
          f1 = f2;
        } // for (m = 0; ...)
        // Insert sface between f1 and f2.
        newlinkitem = (badface *) flippool->alloc();
        newlinkitem->ss = sface;
        newlinkitem->nextitem = f1->nextitem;
        f1->nextitem = newlinkitem;
      } else if (flippool->items == 1) {
        f1 = facelink;
        // Make sure that f is not coplanar and codirection with f1.
        ori1 = orient3d(torg, tdest, sapex(f1->ss), sapex(sface));
        if (ori1 == 0) {
          // f is coplanar with f1 (see Fig. 8).
          facenormal(torg, tdest, sapex(f1->ss), n1, 1);
          facenormal(torg, tdest, sapex(sface), n2, 1);
          if (DOT(n1, n2) > 0) {
            // The two faces are codirectional as well.
            assert(0);
          }
        }
        // Add this face into link.
        newlinkitem = (badface *) flippool->alloc();
        newlinkitem->ss = sface;
        newlinkitem->nextitem = NULL;
        f1->nextitem = newlinkitem;
      } else {
        // The first face.
        newlinkitem = (badface *) flippool->alloc();
        newlinkitem->ss = sface;
        newlinkitem->nextitem = NULL;
        facelink = newlinkitem;
      }
    } // for (k = idx2faclist[idx]; ...)

    if (b->verbose > 1) {
      printf("    Found %ld segments at (%d  %d).\n", flippool->items,
             pointmark(torg), pointmark(tdest));
    }

    // Set the connection between this segment and faces containing it,
    //   at the same time, remove redundant segments.
    f1 = facelink;
    for (k = 0; k < flippool->items; k++) {
      sspivot(f1->ss, testseg);
      // If 'testseg' is not 'subsegloop' and is not dead, it is redundant.
      if ((testseg.sh != subsegloop.sh) && (testseg.sh[3] != NULL)) {
        shellfacedealloc(subsegpool, testseg.sh);
      }
      // Bonds the subface and the segment together.
      ssbond(f1->ss, subsegloop);
      f1 = f1->nextitem;
    }

    // Create the face ring at the segment.
    f1 = facelink;
    for (k = 1; k <= flippool->items; k++) {
      k < flippool->items ? f2 = f1->nextitem : f2 = facelink;
      if (b->verbose > 2) {
        printf("    Bond subfaces (%d, %d, %d) and (%d, %d, %d).\n",
               pointmark(torg), pointmark(tdest), pointmark(sapex(f1->ss)),
               pointmark(torg), pointmark(tdest), pointmark(sapex(f2->ss)));
      }
      sbond1(f1->ss, f2->ss);
      f1 = f2;
    }

    // Set the unique segment marker into the unified segment.
    shellmark(subsegloop) = segmarker;
    segmarker++;
    // sfacelist->clear();
    flippool->restart();

    subsegloop.sh = shellfacetraverse(subsegpool);
  }

  delete [] idx2faclist;
  delete [] facperverlist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// mergefacets()    Merge adjacent coplanar facets.                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::mergefacets()
{
  face parentsh, neighsh, neineighsh;
  face segloop;
  point eorg, edest;
  REAL ori;
  bool mergeflag;
  int* segspernodelist;
  int fidx1, fidx2;
  int i, j;

  if (b->verbose) {
    printf("  Merging adjacent coplanar facets.\n");
  }

  // Initialize 'segspernodelist'.
  segspernodelist = new int[pointpool->items + 1];
  for (i = 0; i < pointpool->items + 1; i++) segspernodelist[i] = 0;

  // Loop all segments, counter the number of segments sharing each vertex.
  subsegpool->traversalinit();
  segloop.sh = shellfacetraverse(subsegpool);
  while (segloop.sh != (shellface *) NULL) {
    // Increment the number of sharing segments for each endpoint.
    for (i = 0; i < 2; i++) {
      j = pointmark((point) segloop.sh[3 + i]);
      segspernodelist[j]++;
    }
    segloop.sh = shellfacetraverse(subsegpool);
  }

  // Loop all segments, merge adjacent coplanar facets.
  subsegpool->traversalinit();
  segloop.sh = shellfacetraverse(subsegpool);
  while (segloop.sh != (shellface *) NULL) {
    eorg = sorg(segloop);
    edest = sdest(segloop);
    spivot(segloop, parentsh);
    spivot(parentsh, neighsh);
    spivot(neighsh, neineighsh);
    if ((parentsh.sh != neighsh.sh) && (parentsh.sh == neineighsh.sh)) {
      // Exactly two subfaces at this segment.
      fidx1 = shellmark(parentsh) - 1;
      fidx2 = shellmark(neighsh) - 1;
      // Possible to merge them if they are not in the same facet.
      if (fidx1 != fidx2) {
        // Test if they are coplanar wrt the tolerance.
        ori = orient3d(eorg, edest, sapex(parentsh), sapex(neighsh));
        if ((ori == 0) || iscoplanar(eorg, edest, sapex(parentsh),
          sapex(neighsh), ori)) {
          // Found two adjacent coplanar facets.
          mergeflag = ((in->facetmarkerlist == NULL) || 
            (in->facetmarkerlist[fidx1] == in->facetmarkerlist[fidx2]));
          if (mergeflag) {
            if (b->verbose > 1) {
              printf("  Removing segment (%d, %d).\n", pointmark(eorg),
                     pointmark(edest));
            }
            ssdissolve(parentsh);
            ssdissolve(neighsh);
            shellfacedealloc(subsegpool, segloop.sh);
            j = pointmark(eorg);
            segspernodelist[j]--;
            if (segspernodelist[j] == 0) {
              pointtype(eorg) = FACETVERTEX;
            }
            j = pointmark(edest);
            segspernodelist[j]--;
            if (segspernodelist[j] == 0) {
              pointtype(edest) = FACETVERTEX;
            }
            // Add the edge to flip stack.
            futureflip = flipshpush(futureflip, &parentsh);
          }
        }
      }
    }
    segloop.sh = shellfacetraverse(subsegpool);
  }

  if (futureflip != NULL) {
    // Do Delaunay flip.
    lawsonflip();
  }

  delete [] segspernodelist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// meshsurface()    Create a surface mesh of the input PLC.                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::meshsurface()
{
  arraypool *ptlist, *conlist;
  point *idx2verlist;
  point tstart, tend, *pnewpt, *cons;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  int *worklist;
  int end1, end2;
  int shmark, i, j;

  if (!b->quiet) {
    printf("Creating surface mesh.\n");
  }

  // Initialize dynamic arrays (length: 2^8 = 256).
  ptlist = new arraypool(sizeof(point *), 8);
  conlist = new arraypool(2 * sizeof(point *), 8);
  worklist = new int[pointpool->items + 1];
  for (i = 0; i < pointpool->items + 1; i++) worklist[i] = 0;

  // Compute a mapping from indices to points.
  makeindex2pointmap(idx2verlist);

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
        ptlist->newindex((void **) &pnewpt);
        *pnewpt = tstart;
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
              ptlist->newindex((void **) &pnewpt);
              *pnewpt = tend;
              worklist[end2] = 1;
            }
            // Save the segment in S (conlist).
            conlist->newindex((void **) &cons);
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
    for (i = 0; i < ptlist->objects; i++) {
      pnewpt = (point *) fastlookup(ptlist, i);
      end1 = pointmark(*pnewpt);
      worklist[end1] = 0;
    }

    // Triangulate F into a CDT.
    triangulate(shmark, ptlist, conlist, f->numberofholes, f->holelist);

    // Clear working lists.
    ptlist->restart();
    conlist->restart();
  }

  delete ptlist;
  delete conlist;
  delete [] worklist;
  delete [] idx2verlist;

  // Remove redundant segments and build the face links.
  unifysegments();

  if (b->object == tetgenbehavior::STL) {
    // Remove redundant vertices (for .stl input mesh).
    jettisonnodes();
  }

  if (!b->nomerge && !b->nobisect) {
    // Merge adjacent coplanar facets.
    mergefacets();
  }

  // The total number of iunput segments.
  insegments = subsegpool->items;
}

#endif // #ifndef surfaceCXX