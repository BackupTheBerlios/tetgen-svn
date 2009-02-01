#ifndef surfaceCXX
#define surfaceCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// sinsertvertex()    Insert a vertex into a triangulation of a facet.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::sinsertvertex(point insertpt,face *searchsh, face *splitseg, 
  bool bwflag, bool hullflag)
{
  face *abfaces, *parysh, *pssub;
  face neighsh, newsh, casout, casin;
  face aseg, bseg, aoutseg, boutseg;
  face checkseg;
  triface neightet;
  point pa, pb, pc;
  enum location loc;
  REAL sign, area;
  int n, s, i, j;

  tetrahedron ptr;
  shellface sptr;

  if (splitseg == NULL) {
    assert(searchsh != NULL); // SELF_CHECK
    // loc = slocate(insertpt, searchsh);
  } else {
    spivot(*splitseg, *searchsh);
    loc = ONEDGE;
  }

  if (loc == OUTSIDE) {
    if (hullflag) {
      // return inserthullvertex();
    }
    assert(0); 
  }
  // The insert point should not lie on an unknown segment.
  if ((loc == ONEDGE) && (splitseg->sh == NULL)) {
    sspivot(*searchsh, checkseg);
    assert(checkseg.sh == NULL); // SELF_CHECK
  }
  if (loc == ONVERTEX) {
    assert(0);
  }

  if (b->verbose > 1) {
    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    pc = sapex(*searchsh);
    printf("    Insert point %d %s (%d, %d, %d)\n", pointmark(insertpt),
      loc == ONEDGE ? "on edge" : "in face", pointmark(pa), pointmark(pb), 
      pointmark(pc));
  }

  // Does 'insertpt' lie on a segment?
  if (splitseg != NULL) {
    splitseg->shver = 0;
    pa = sorg(*splitseg);
    // Count the number of faces at segment [a, b].
    n = 0;
    neighsh = *searchsh;
    do {
      spivotself(neighsh);
      n++;
    } while ((neighsh.sh != NULL) && (neighsh.sh != searchsh->sh));
    // n is at least 1.
    abfaces = new face[n];
    // Collect faces at seg [a, b].
    abfaces[0] = *searchsh;
    if (sorg(abfaces[0]) != pa) sesymself(abfaces[0]);
    for (i = 1; i < n; i++) {
      spivot(abfaces[i - 1], abfaces[i]);
      if (sorg(abfaces[i]) != pa) sesymself(abfaces[i]);
    }
  }

  // Initialize the cavity.
  smarktest(*searchsh);
  caveshlist->newindex((void **) &parysh);
  *parysh = *searchsh;
  if (loc == ONEDGE) {
    if (splitseg != NULL) {
      for (i = 1; i < n; i++) {
        smarktest(abfaces[i]);
        caveshlist->newindex((void **) &parysh);
        *parysh = abfaces[i];
      }
    } else {
      spivot(*searchsh, neighsh);
      if (neighsh.sh != NULL) {
        smarktest(neighsh);
        caveshlist->newindex((void **) &parysh);
        *parysh = neighsh;
      }
    }
  }

  // Form the Bowyer-Watson cavity.
  for (i = 0; i < caveshlist->objects; i++) {
    parysh = (face *) fastlookup(caveshlist, i);
    for (j = 0; j < 3; j++) {
      sspivot(*parysh, checkseg);
      if (checkseg.sh == NULL) {
        spivot(*parysh, neighsh);
        assert(neighsh.sh != NULL); // SELF_CHECK
        if (!smarktested(neighsh)) {
          if (bwflag) {
            pa = sorg(neighsh);
            pb = sdest(neighsh);
            pc = sapex(neighsh);
            sign = incircle3d(pa, pb, pc, insertpt);
            if (sign < 0) {
              smarktest(neighsh);
              caveshlist->newindex((void **) &pssub);
              *pssub = neighsh;
            }
          } else {
            sign = 1; // A boundary edge.
          }
        } else {
          sign = -1; // Not a boundary edge.
        }
      } else {
        sign = 1; // A segment!
      }
      if (sign >= 0) {
        // Add a boundary edge.
        caveshbdlist->newindex((void **) &pssub);
        *pssub = *parysh;
      }
      senextself(*parysh);
    }
  }

  // Creating new subfaces.
  for (i = 0; i < caveshbdlist->objects; i++) {
    parysh = (face *) fastlookup(caveshbdlist, i);
    sspivot(*parysh, checkseg);
    if ((parysh->shver & 01) != 0) sesymself(*parysh);
    pa = sorg(*parysh);
    pb = sdest(*parysh);
    // Create a new subface.
    makeshellface(subfacepool, &newsh); 
    setshvertices(newsh, pa, pb, insertpt);
    setshellmark(newsh, getshellmark(*parysh));
    if (checkconstraints) {
      area = areabound(*parysh);
      areabound(newsh) = area;
    }
    // Connect newsh to outer subfaces.
    spivot(*parysh, casout);
    if (casout.sh != NULL) {
      casin = casout;
      if (checkseg.sh != NULL) {
        spivot(casin, neighsh);
        while (neighsh.sh != parysh->sh) {
          casin = neighsh;
          spivot(casin, neighsh);
        }
      }
      sbond1(newsh, casout);
      sbond1(casin, newsh);
    }
    if (checkseg.sh != NULL) {
      ssbond(newsh, checkseg);
    }
    // Connect oldsh <== newsh (for connecting adjacent new subfaces).
    sbond1(*parysh, newsh);
  }

  // Connect adjacent new subfaces together.
  for (i = 0; i < caveshbdlist->objects; i++) {
    // Get an old subface at edge [a, b].
    parysh = (face *) fastlookup(caveshbdlist, i);
    sspivot(*parysh, checkseg);
    spivot(*parysh, newsh); // The new subface [a, b, p].
    senextself(newsh); // At edge [b, p].
    spivot(newsh, neighsh);
    if (neighsh.sh == NULL) {
      // Find the adjacent new subface at edge [b, p].
      pb = sdest(*parysh);
      neighsh = *parysh;
      while (1) {
        senextself(neighsh);
        spivotself(neighsh);
        if (!smarktested(neighsh)) break;
        if (sdest(neighsh) != pb) sesymself(neighsh);
      }
      // Now 'neighsh' is a new subface at edge [b, #].
      if (sorg(neighsh) != pb) sesymself(neighsh);
      assert(sorg(neighsh) == pb); // SELF_CHECK
      assert(sapex(neighsh) == insertpt); // SELF_CHECK
      senext2self(neighsh); // Go to the open edge [p, b].
      spivot(neighsh, casout); // SELF_CHECK
      assert(casout.sh == NULL); // SELF_CHECK
      sbond2(newsh, neighsh);
    }
    spivot(*parysh, newsh); // The new subface [a, b, p].
    senext2self(newsh); // At edge [p, a].
    spivot(newsh, neighsh);
    if (neighsh.sh == NULL) {
      // Find the adjacent new subface at edge [p, a].
      pa = sorg(*parysh);
      neighsh = *parysh;
      while (1) {
        senext2self(neighsh);
        spivotself(neighsh);
        if (!smarktested(neighsh)) break;
        if (sorg(neighsh) != pa) sesymself(neighsh);
      }
      // Now 'neighsh' is a new subface at edge [#, a].
      if (sdest(neighsh) != pa) sesymself(neighsh);
      assert(sdest(neighsh) == pa); // SELF_CHECK
      assert(sapex(neighsh) == insertpt); // SELF_CHECK
      senextself(neighsh); // Go to the open edge [a, p].
      spivot(neighsh, casout); // SELF_CHECK
      assert(casout.sh == NULL); // SELF_CHECK
      sbond2(newsh, neighsh);
    }
  }

  if (splitseg != NULL) {
    // Split the segment [a, b].
    aseg = *splitseg;
    pa = sorg(aseg);
    pb = sdest(aseg);
    if (b->verbose > 1) {
      printf("    Split seg (%d, %d) by %d.\n", pointmark(pa), pointmark(pb), 
        pointmark(insertpt));
    }
    // Insert the new point p.
    makeshellface(subsegpool, &bseg);
    setshvertices(bseg, insertpt, pb, NULL);
    setsdest(aseg, insertpt);
    setshellmark(bseg, getshellmark(aseg));
    if (checkconstraints) {
      areabound(bseg) = areabound(aseg);
    }
    // Connect [p, b]<->[b, #].
    senext(aseg, aoutseg);
    spivotself(aoutseg);
    if (aoutseg.sh != NULL) {
      senext(bseg, boutseg);
      sbond2(boutseg, aoutseg);
    }
    // Connect [a, p] <-> [p, b].
    senext(aseg, aoutseg);
    senext2(bseg, boutseg);
    sbond2(aoutseg, boutseg);
    // Connect subsegs [a, p] and [p, b] to the true new subfaces.
    for (i = 0; i < n; i++) {
      spivot(abfaces[i], newsh); // The faked new subface.
      if (sorg(newsh) != pa) sesymself(newsh);
      senext2(newsh, neighsh); // The edge [p, a] in newsh
      spivot(neighsh, casout);
      ssbond(casout, aseg);
      senext(newsh, neighsh); // The edge [b, p] in newsh
      spivot(neighsh, casout);
      ssbond(casout, bseg);
    }
    if (n > 1) {
      // Create the two face rings at [a, p] and [p, b].
      for (i = 0; i < n; i++) {
        spivot(abfaces[i], newsh); // The faked new subface.
        if (sorg(newsh) != pa) sesymself(newsh);
        spivot(abfaces[(i + 1) % n], neighsh); // The next faked new subface.
        if (sorg(neighsh) != pa) sesymself(neighsh);
        senext2(newsh, casout); // The edge [p, a] in newsh.
        senext2(neighsh, casin); // The edge [p, a] in neighsh.
        spivotself(casout);
        spivotself(casin);
        sbond1(casout, casin); // Let the i's face point to (i+1)'s face.
        senext(newsh, casout); // The edge [b, p] in newsh.
        senext(neighsh, casin); // The edge [b, p] in neighsh.
        spivotself(casout);
        spivotself(casin);
        sbond1(casout, casin);
      }
    } else { 
      // Only one subface contains this segment.
      // assert(n == 1);
      spivot(abfaces[0], newsh);  // The faked new subface.
      if (sorg(newsh) != pa) sesymself(newsh);
      senext2(newsh, casout); // The edge [p, a] in newsh.
      spivotself(casout);
      sdissolve(casout); // Disconnect to faked subface.
      senext(newsh, casout); // The edge [b, p] in newsh.
      spivotself(casout);
      sdissolve(casout); // Disconnect to faked subface.
    }
    // Delete the faked new subfaces.
    for (i = 0; i < n; i++) {
      spivot(abfaces[i], newsh); // The faked new subface.
      shellfacedealloc(subfacepool, newsh.sh);
    }
    if (checksubsegs) {
      // Add two subsegs into stack (for recovery).
      s = randomnation(subsegstack->objects + 1);
      subsegstack->newindex((void **) &parysh);
      *parysh = * (face *) fastlookup(subsegstack, s);
      sinfect(aseg); 
      parysh = (face *) fastlookup(subsegstack, s);
      *parysh = aseg;
      s = randomnation(subsegstack->objects + 1);
      subsegstack->newindex((void **) &parysh);
      *parysh = * (face *) fastlookup(subsegstack, s);
      sinfect(bseg);
      parysh = (face *) fastlookup(subsegstack, s);
      *parysh = bseg;
    }
    delete [] abfaces;
  }

  // Delete the old subfaces.
  for (i = 0; i < caveshlist->objects; i++) {
    parysh = (face *) fastlookup(caveshlist, i);
    if (checksubfaces) {
      // Disconnect in the neighbor tets.
      stpivot(*parysh, neightet);
      if (neightet.tet != NULL) {
        tsdissolve(neightet);
        symself(neightet);
        tsdissolve(neightet);
      }
    }
    shellfacedealloc(subfacepool, parysh->sh);
  }

  // Clean the working lists.
  caveshlist->restart();
  caveshbdlist->restart();
}

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
    setshellmark(newsh, shmark);
    // Create three new segments.
    for (i = 0; i < 3; i++) {
      makeshellface(subsegpool, &newseg);
      setshvertices(newseg, sorg(newsh), sdest(newsh), NULL);
      ssbond(newsh, newseg);
      senextself(newsh);
    }
  } else {
    printf("  This code does not do surface mesh yet.\n");
    // terminatetetgen(1);
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

  if (b->verbose > 1) {
    printf("  Unifying segments.\n");
  }

  // Create a mapping from vertices to subfaces.
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
              // f is coplanar and with f1 (see Fig. 7).
              ori3 = orient3d(torg, tdest, sapex(f2->ss), sapex(sface));
              if (ori3 > 0) {
                assert(0);  // f is also codirection with f1.
              } else {
                // f is above f2, continue.
              }
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
    if (flippool->items > 1) {
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
    }

    // Set the unique segment marker into the unified segment.
    setshellmark(subsegloop, segmarker);
    segmarker++;
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

  if (b->verbose > 1) {
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
    if (neighsh.sh != NULL) {
      spivot(neighsh, neineighsh);
      if (parentsh.sh == neineighsh.sh) {
        // Exactly two subfaces at this segment.
        fidx1 = getshellmark(parentsh) - 1;
        fidx2 = getshellmark(neighsh) - 1;
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
                setpointtype(eorg, FACETVERTEX);
              }
              j = pointmark(edest);
              segspernodelist[j]--;
              if (segspernodelist[j] == 0) {
                setpointtype(edest, FACETVERTEX);
              }
              // Add the edge to flip stack.
              futureflip = flipshpush(futureflip, &parentsh);
            }
          }
        }
      }
    } // if (neighsh.sh != NULL)
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
          if (getpointtype(tstart) == DUPLICATEDVERTEX) {
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