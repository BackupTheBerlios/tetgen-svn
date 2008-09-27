#ifndef surfaceCXX
#define surfaceCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// iscoplanar()    Check if four points are approximately coplanar.          //
//                                                                           //
// 'eps' is the relative error tolerance.  The coplanarity is determined by  //
// the equation q < tol, where q = fabs(6 * vol) / L^3, vol is the volume of //
// the tet klmn, and L is the average edge length of the tet.                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::iscoplanar(point k, point l, point m, point n, REAL tol)
{
  REAL ori, L, q;

  ori = orient3d(k, l, m, n);
  if (ori == 0.0) return true;

  L = DIST(k, l);
  L += DIST(l, m);
  L += DIST(m, k);
  L += DIST(k, n);
  L += DIST(l, n);
  L += DIST(m, n);
  assert(L > 0.0);  // SELF_CHECK
  L /= 6.0;
  q = fabs(ori) / (L * L * L);
  return q <= tol;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipshpush()    Push a subface edge into flip stack.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::badface* tetgenmesh::flipshpush(badface* flipstack, face* flipedge)
{
  badface *newflipface;

  newflipface = (badface *) flippool->alloc();
  newflipface->ss = *flipedge;
  newflipface->forg = sorg(*flipedge);
  newflipface->fdest = sdest(*flipedge);
  newflipface->nextitem = flipstack;

  return newflipface;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip22()    Remove an edge by transforming 2-to-2 subfaces.               //
//                                                                           //
// 'flipfaces' contains two faces: abc and bad. This routine removes these 2 //
// faces and replaces them by two new faces: cdb and dca.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip22(face* flipfaces, int flipflag)
{
  face bdedges[4], outfaces[4], infaces[4];
  face checkface, checkseg;
  point pa, pb, pc, pd;
  int i;

  pa = sorg(flipfaces[0]);
  pb = sdest(flipfaces[0]);
  pc = sapex(flipfaces[0]);
  pd = sapex(flipfaces[1]);

  if (b->verbose > 1) {
    printf("    flip 2-to-2: (%d, %d, %d, %d)\n", pointmark(pa),
      pointmark(pb), pointmark(pc), pointmark(pd));
  }

  // Collect the four boundary edges.
  senext(flipfaces[0], bdedges[0]);
  senext2(flipfaces[0], bdedges[1]);
  senext(flipfaces[1], bdedges[2]);
  senext2(flipfaces[1], bdedges[3]);

  // Collect outer boundary faces.
  for (i = 0; i < 4; i++) {
    spivot(bdedges[i], outfaces[i]);
    infaces[i] = outfaces[i];
    sspivot(bdedges[i], checkseg);
    if (checkseg.sh != NULL) {
      spivot(infaces[i], checkface);
      while (checkface.sh != bdedges[i].sh) {
        infaces[i] = checkface;
        spivot(infaces[i], checkface);
      }
    }
  }

  // Transform abc -> cdb.
  setshvertices(flipfaces[0], pc, pd, pb);
  // Transform bad -> dca.
  setshvertices(flipfaces[1], pd, pc, pa);

  // Reconnect boundary edges to outer boundary faces.
  for (i = 0; i < 4; i++) {
    sbond1(bdedges[i], outfaces[(3 + i) % 4]);
    sbond1(infaces[(3 + i) % 4], bdedges[i]);
    sspivot(outfaces[(3 + i) % 4], checkseg); // checkseg may be NULL.
    ssbond(bdedges[i], checkseg); // Clear the old bond as well.
  }

  if (flipflag) {
    // Put the boundary edges into flip stack.
    for (i = 0; i < 4; i++) {
      futureflip = flipshpush(futureflip, &bdedges[i]);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lawsonflip()    Flip non-locally Delaunay edges.                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::lawsonflip()
{
  face flipfaces[2];
  face checkseg;
  point pa, pb, pc, pd;
  REAL sign;
  long flipcount;

  if (b->verbose > 1) {
    printf("    Lawson flip %ld edges.\n", flippool->items);
  }
  flipcount = 0l;

  while (futureflip != (badface *) NULL) {
    // Pop an edge from the stack.
    flipfaces[0] = futureflip->ss;
    pa = futureflip->forg;
    pb = futureflip->fdest;
    futureflip = futureflip->nextitem;

    // Skip it if it is dead.
    if (flipfaces[0].sh[3] == NULL) continue;
    // Skip it if it is not the same edge as we saved.
    if ((sorg(flipfaces[0]) != pa) || (sdest(flipfaces[0]) != pb)) continue;
    // Skip it if it is a subsegment.
    sspivot(flipfaces[0], checkseg);
    if (checkseg.sh != NULL) continue;

    // Get the adjacent face.
    spivot(flipfaces[0], flipfaces[1]);
    sesymself(flipfaces[1]);
    pc = sapex(flipfaces[0]);
    pd = sapex(flipfaces[1]);

    // sign = incircle3d(pa, pb, pc, pd);

    if (sign < 0) {
      // It is non-locally Delaunay. Flip it.
      flip22(flipfaces, 1);
      flipcount++;
    }
  }

  if (b->verbose > 1) {
    printf("    %ld flips.\n", flipcount);
  }

  flippool->restart();
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// triangulate()    Create a CDT for the facet.                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::triangulate(int shmark, list* ptlist, list* conlist,
  int holes, REAL* holelist)
{
  face newsh, newseg;
  point *pts;
  int i;

  if (b->verbose > 1) {
    printf("    %d vertices, %d segments", ptlist->len(), conlist->len());
    if (holes > 0) {
      printf(", %d holes", holes);
    }
    printf(", shmark: %d.\n", shmark);
  }

  if ((ptlist->len() == 3) && (conlist->len() == 3)) {
    // This CDT contains only one triangle.
    pts = (point *) ptlist->base;
    makeshellface(subfacepool, &newsh);
    setshvertices(newsh, pts[0], pts[1], pts[2]);
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
  list *sfacelist;
  face *facperverlist;
  face subsegloop, testseg;
  face sface, sface1, sface2;
  point torg, tdest;
  REAL ori1, ori2, ori3;
  int *idx2faclist;
  int segmarker;
  int idx, k, m;

  if (b->verbose) {
    printf("  Unifying segments.\n");
  }

  // Initialize a list for storing the face link at a segment.
  sfacelist = new list(sizeof(face), NULL); 

  // Create a mapping from vertices to subfaces incident at them.
  makesubfacemap(idx2faclist, facperverlist);

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
      assert(sdest(sface) == tdest); // SELF_CHECK
      // Save the face f in 'sfacelist'.
      if (sfacelist->len() > 2) {
        for (m = 0; m < sfacelist->len() - 1; m++) {
          sface1 = * (face *)(* sfacelist)[m];
          sface2 = * (face *)(* sfacelist)[m + 1];
          ori1 = orient3d(torg, tdest, sapex(sface1), sapex(sface2));
          ori2 = orient3d(torg, tdest, sapex(sface1), sapex(sface));
          if (ori1 > 0) {
            if (ori2 > 0) {
              // Both apex(f), apex(f2) are below f1 (see Fig.1). 
              ori3 = orient3d(torg, tdest, sapex(sface2), sapex(sface));
              if (ori3 > 0) {
                // f is after both f1 and f2, continue. 
              } else if (ori3 < 0) {
                // f is between f1 and f2.
                break; 
              } else {
                // f is duplicated with f2. Not handled yet.
                assert(0);  // ori3 == 0; 
              }
            } else if (ori2 < 0) {
              // apex(f) is above f1, continue (see Fig. 2).
            } else { // ori2 == 0;
              // f is duplicated with f1. Not handled yet.
              assert(0); 
            }
          } else if (ori1 < 0) {
            if (ori2 > 0) {
              // apex(f) is between f1 and f2 (see Fig. 3).
              break;
            } else if (ori2 < 0) {
              // Both apex(f), apex(f2) are below f1 (see Fig.4).
              ori3 = orient3d(torg, tdest, sapex(sface2), sapex(sface));
              if (ori3 > 0) {
                // f is between f1 and f2 (see Fig.4).
                break;
              } else if (ori3 < 0) {
                // f is after f1 and f2, continue.
              } else { // ori3 == 0;
                // f is duplicated with f2. Not handled yet.
                assert(0);
              }
            } else { // ori2 == 0;
              // f is duplicated with f2. Not handled yet.
              assert(0);
            }
          } else { // ori1 == 0;
            // f is duplicated with f1. Not handled yet.
            assert(0);
          }
        } // for (m = 0; ...)
        sfacelist->insert(m + 1, &sface);
      } else {
        sfacelist->append(&sface);
      }
    } // for (k = idx2faclist[idx]; ...)

    if (b->verbose > 1) {
      printf("    Found %d segments at (%d  %d).\n", sfacelist->len(),
             pointmark(torg), pointmark(tdest));
    }

    // Set the connection between this segment and faces containing it,
    //   at the same time, remove redundant segments.
    for (k = 0; k < sfacelist->len(); k++) {
      sface = *(face *)(* sfacelist)[k];
      sspivot(sface, testseg);
      // If 'testseg' is not 'subsegloop' and is not dead, it is redundant.
      if ((testseg.sh != subsegloop.sh) && (testseg.sh[3] != NULL)) {
        shellfacedealloc(subsegpool, testseg.sh);
      }
      // Bonds the subface and the segment together.
      ssbond(sface, subsegloop);
    }
    // Set connection between these faces.
    sface = *(face *)(* sfacelist)[0];
    for (k = 1; k <= sfacelist->len(); k++) {
      if (k < sfacelist->len()) {
        sface1 = *(face *)(* sfacelist)[k];
      } else {
        sface1 = *(face *)(* sfacelist)[0];    // Form a face loop.
      }
      if (b->verbose > 2) {
        printf("    Bond subfaces (%d, %d, %d) and (%d, %d, %d).\n",
               pointmark(torg), pointmark(tdest), pointmark(sapex(sface)),
               pointmark(torg), pointmark(tdest), pointmark(sapex(sface1)));
      }
      sbond1(sface, sface1);
      sface = sface1;
    }

    // Set the unique segment marker into the unified segment.
    shellmark(subsegloop) = segmarker;
    segmarker++;
    sfacelist->clear();

    subsegloop.sh = shellfacetraverse(subsegpool);
  }

  // The total number of iunput segments.
  insegments = subsegpool->items;

  delete [] idx2faclist;
  delete [] facperverlist;
  delete sfacelist;
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
        if (iscoplanar(eorg, edest, sapex(parentsh), sapex(neighsh),
                       b->epsilon)) {
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
              pointtype(eorg) = FREESUBVERTEX;
            }
            j = pointmark(edest);
            segspernodelist[j]--;
            if (segspernodelist[j] == 0) {
              pointtype(edest) = FREESUBVERTEX;
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
  list *ptlist, *conlist;
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

    // Triangulate F into a CDT.
    triangulate(shmark, ptlist, conlist, f->numberofholes, f->holelist);

    // Clear working lists.
    ptlist->clear();
    conlist->clear();
  }

  // Remove redundant segments and build the face links.
  unifysegments();

  if (b->object == tetgenbehavior::STL) {
    // Remove redundant vertices (for .stl input mesh).
    jettisonnodes();
  }

  if (!b->nomerge && !b->nobisect) {
    // Merge adjacent coplanar facets.
    // mergefacets();
  }

  delete ptlist;
  delete conlist;
  delete [] worklist;
  delete [] idx2verlist;
}

#endif // #ifndef surfaceCXX