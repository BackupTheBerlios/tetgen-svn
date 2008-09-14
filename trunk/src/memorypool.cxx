#ifndef memorypoolCXX
#define memorypoolCXX

#include "tetgen.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// append()    Add a new item at the end of the list.                        //
//                                                                           //
// A new space at the end of this list will be allocated for storing the new //
// item. If the memory is not sufficient, reallocation will be performed. If //
// 'appitem' is not NULL, the contents of this pointer will be copied to the //
// new allocated space.  Returns the pointer to the new allocated space.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::list::append(void *appitem)
{
  // Do we have enough space?
  if (items == maxitems) {
    char* newbase = (char *) realloc(base, (maxitems + expandsize) * 
                                     itembytes);
    if (newbase == (char *) NULL) {
      printf("Error:  Out of memory.\n");
      terminatetetgen(1);
    }
    base = newbase;
    maxitems += expandsize;
  }
  if (appitem != (void *) NULL) {
    memcpy(base + items * itembytes, appitem, itembytes);
  }
  items++;
  return (void *) (base + (items - 1) * itembytes);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insert()    Insert an item before 'pos' (range from 0 to items - 1).      //
//                                                                           //
// A new space will be inserted at the position 'pos', that is, items lie    //
// after pos (including the item at pos) will be moved one space downwords.  //
// If 'insitem' is not NULL, its contents will be copied into the new        //
// inserted space. Return a pointer to the new inserted space.               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::list::insert(int pos, void* insitem)
{
  if (pos >= items) {
    return append(insitem);
  }
  // Do we have enough space.
  if (items == maxitems) {
    char* newbase = (char *) realloc(base, (maxitems + expandsize) *
                                     itembytes);
    if (newbase == (char *) NULL) {
      printf("Error:  Out of memory.\n");
      terminatetetgen(1);
    }
    base = newbase;
    maxitems += expandsize;
  }
  // Do block move.
  memmove(base + (pos + 1) * itembytes,   // dest
          base + pos * itembytes,         // src
          (items - pos) * itembytes);     // size in bytes
  // Insert the item.
  if (insitem != (void *) NULL) {
    memcpy(base + pos * itembytes, insitem, itembytes);
  }
  items++;
  return (void *) (base + pos * itembytes);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// del()    Delete an item at 'pos' (range from 0 to items - 1).             //
//                                                                           //
// The space at 'pos' will be overlapped by other item. If 'order' is 1, the //
// remaining items of the list have the same order as usual, i.e., items lie //
// after pos will be moved one space upwords. If 'order' is 0, the last item //
// of the list will be moved up to pos.                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::list::del(int pos, int order)
{
  // If 'pos' is the last item of the list, nothing need to do.
  if (pos >= 0 && pos < items - 1) {
    if (order == 1) {
      // Do block move. 
      memmove(base + pos * itembytes,       // dest
              base + (pos + 1) * itembytes, // src
              (items - pos - 1) * itembytes);
    } else {
      // Use the last item to overlap the del item.
      memcpy(base + pos * itembytes, // item at pos
             base + (items - 1) * itembytes, // item at last
             itembytes);
    }
  }
  if (items > 0) {
    items--;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// list()    Initialize the list.                                            //
//                                                                           //
// Set the size (in bytes) of each item, the maximum size allocated at once, //
// onece, the expand size when the list is full.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::list::list(int itbytes, int mitems, int exsize)
{
  // assert(itbytes > 0 && mitems > 0 && exsize > 0);
  itembytes = itbytes;
  maxitems = mitems;
  expandsize = exsize;
  base = (char *) malloc(maxitems * itembytes); 
  if (base == (char *) NULL) {
    printf("Error:  Out of memory.\n");
    terminatetetgen(1);
  }
  items = 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// restart()   Deallocate all items in this pool.                            //
//                                                                           //
// The pool is returned to its starting state, except that no memory is      //
// freed to the operating system.  Rather, the previously allocated blocks   //
// are ready to be reused.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::memorypool::restart()
{
  unsigned long alignptr;

  items = 0;
  maxitems = 0;

  // Set the currently active block.
  nowblock = firstblock;
  // Find the first item in the pool.  Increment by the size of (void *).
  alignptr = (unsigned long) (nowblock + 1);
  // Align the item on an `alignbytes'-byte boundary.
  nextitem = (void *)
    (alignptr + (unsigned long) alignbytes -
     (alignptr % (unsigned long) alignbytes));
  // There are lots of unallocated items left in this block.
  unallocateditems = itemsperblock;
  // The stack of deallocated items is empty.
  deaditemstack = (void *) NULL;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// alloc()   Allocate space for an item.                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::memorypool::alloc()
{
  void *newitem;
  void **newblock;
  unsigned long alignptr;

  // First check the linked list of dead items.  If the list is not 
  //   empty, allocate an item from the list rather than a fresh one.
  if (deaditemstack != (void *) NULL) {
    newitem = deaditemstack;                     // Take first item in list.
    deaditemstack = * (void **) deaditemstack;
  } else {
    // Check if there are any free items left in the current block.
    if (unallocateditems == 0) {
      // Check if another block must be allocated.
      if (*nowblock == (void *) NULL) {
        // Allocate a new block of items, pointed to by the previous block.
        newblock = (void **) malloc(itemsperblock * itembytes + sizeof(void *) 
                                    + alignbytes);
        if (newblock == (void **) NULL) {
          printf("Error:  Out of memory.\n");
          terminatetetgen(1);
        }
        *nowblock = (void *) newblock;
        // The next block pointer is NULL.
        *newblock = (void *) NULL;
      }
      // Move to the new block.
      nowblock = (void **) *nowblock;
      // Find the first item in the block.
      //   Increment by the size of (void *).
      alignptr = (unsigned long) (nowblock + 1);
      // Align the item on an `alignbytes'-byte boundary.
      nextitem = (void *)
        (alignptr + (unsigned long) alignbytes -
         (alignptr % (unsigned long) alignbytes));
      // There are lots of unallocated items left in this block.
      unallocateditems = itemsperblock;
    }
    // Allocate a new item.
    newitem = nextitem;
    // Advance `nextitem' pointer to next free item in block.
    if (itemwordtype == POINTER) {
      nextitem = (void *) ((void **) nextitem + itemwords);
    } else {
      nextitem = (void *) ((REAL *) nextitem + itemwords);
    }
    unallocateditems--;
    maxitems++;
  }
  items++;
  return newitem;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// dealloc()   Deallocate space for an item.                                 //
//                                                                           //
// The deallocated space is stored in a queue for later reuse.               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::memorypool::dealloc(void *dyingitem)
{
  // Push freshly killed item onto stack.
  *((void **) dyingitem) = deaditemstack;
  deaditemstack = dyingitem;
  items--;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// traversalinit()   Prepare to traverse the entire list of items.           //
//                                                                           //
// This routine is used in conjunction with traverse().                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::memorypool::traversalinit()
{
  unsigned long alignptr;

  // Begin the traversal in the first block.
  pathblock = firstblock;
  // Find the first item in the block.  Increment by the size of (void *).
  alignptr = (unsigned long) (pathblock + 1);
  // Align with item on an `alignbytes'-byte boundary.
  pathitem = (void *)
    (alignptr + (unsigned long) alignbytes -
     (alignptr % (unsigned long) alignbytes));
  // Set the number of items left in the current block.
  pathitemsleft = itemsperblock;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// traverse()   Find the next item in the list.                              //
//                                                                           //
// This routine is used in conjunction with traversalinit().  Be forewarned  //
// that this routine successively returns all items in the list, including   //
// deallocated ones on the deaditemqueue. It's up to you to figure out which //
// ones are actually dead.  It can usually be done more space-efficiently by //
// a routine that knows something about the structure of the item.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::memorypool::traverse()
{
  void *newitem;
  unsigned long alignptr;

  // Stop upon exhausting the list of items.
  if (pathitem == nextitem) {
    return (void *) NULL;
  }
  // Check whether any untraversed items remain in the current block.
  if (pathitemsleft == 0) {
    // Find the next block.
    pathblock = (void **) *pathblock;
    // Find the first item in the block.  Increment by the size of (void *).
    alignptr = (unsigned long) (pathblock + 1);
    // Align with item on an `alignbytes'-byte boundary.
    pathitem = (void *)
      (alignptr + (unsigned long) alignbytes -
       (alignptr % (unsigned long) alignbytes));
    // Set the number of items left in the current block.
    pathitemsleft = itemsperblock;
  }
  newitem = pathitem;
  // Find the next item in the block.
  if (itemwordtype == POINTER) {
    pathitem = (void *) ((void **) pathitem + itemwords);
  } else {
    pathitem = (void *) ((REAL *) pathitem + itemwords);
  }
  pathitemsleft--;
  return newitem;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// poolinit()    Initialize a pool of memory for allocation of items.        //
//                                                                           //
// A 'pool' is created whose records have size at least 'bytecount'.  Items  //
// will be allocated in `itemcount'-item blocks.  Each item is assumed to be //
// a collection of words, and either pointers or floating-point values are   //
// assumed to be the "primary" word type.  (The "primary" word type is used  //
// to determine alignment of items.)  If 'alignment' isn't zero, all items   //
// will be `alignment'-byte aligned in memory.  'alignment' must be either a //
// multiple or a factor of the primary word size;  powers of two are safe.   //
// 'alignment' is normally used to create a few unused bits at the bottom of //
// each item's pointer, in which information may be stored.                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::memorypool::memorypool(int bytecount, int itemcount,
  enum wordtype wtype, int alignment)
{
  int wordsize;

  // Initialize values in the pool.
  itemwordtype = wtype;
  wordsize = (itemwordtype == POINTER) ? sizeof(void *) : sizeof(REAL);
  // Find the proper alignment, which must be at least as large as:
  //   - The parameter `alignment'.
  //   - The primary word type, to avoid unaligned accesses.
  //   - sizeof(void *), so the stack of dead items can be maintained
  //       without unaligned accesses.
  if (alignment > wordsize) {
    alignbytes = alignment;
  } else {
    alignbytes = wordsize;
  }
  if ((int) sizeof(void *) > alignbytes) {
    alignbytes = (int) sizeof(void *);
  }
  itemwords = ((bytecount + alignbytes - 1) /  alignbytes)
            * (alignbytes / wordsize);
  itembytes = itemwords * wordsize;
  itemsperblock = itemcount;

  // Allocate a block of items.  Space for `itemsperblock' items and one
  //   pointer (to point to the next block) are allocated, as well as space
  //   to ensure alignment of the items. 
  firstblock = (void **) malloc(itemsperblock * itembytes + sizeof(void *)
                                + alignbytes); 
  if (firstblock == (void **) NULL) {
    printf("Error:  Out of memory.\n");
    terminatetetgen(1);
  }
  // Set the next block pointer to NULL.
  *(firstblock) = (void *) NULL;
  restart();
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// ~memorypool()   Free to the operating system all memory taken by a pool.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::memorypool::~memorypool()
{
  while (firstblock != (void **) NULL) {
    nowblock = (void **) *(firstblock);
    free(firstblock);
    firstblock = nowblock;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// link()    Initialize a link for storing items.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::link::link(int bytecount, int itemcount)
{
  // assert(bytecount > 0 && itemcount > 0);
  // Remember the real size of each item.
  linkitembytes = bytecount;
  // Set the linear order function for this link.

  // Call the constructor of 'memorypool' to initialize its variables.
  //   like: itembytes, itemwords, items, ... Each node has size
  //   bytecount + 2 * sizeof(void **), and total 'itemcount + 2' (because
  //   link has additional two nodes 'head' and 'tail').
  memorypool(bytecount + 2 * sizeof(void **), itemcount + 2, POINTER, 0);
  
  // Initial state of this link.
  head = (void **) alloc();
  tail = (void **) alloc();
  *head = (void *) tail;
  *(head + 1) = NULL;
  *tail = NULL;
  *(tail + 1) = (void *) head;
  nextlinkitem = *head;
  curpos = 1;
  linkitems = 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// clear()   Deallocate all nodes in this link.                              //
//                                                                           //
// The link is returned to its starting state, except that no memory is      //
// freed to the operating system.  Rather, the previously allocated blocks   //
// are ready to be reused.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::link::clear()
{
  // Reset the pool.
  restart();

  // Initial state of this link.
  head = (void **) alloc();
  tail = (void **) alloc();
  *head = (void *) tail;
  *(head + 1) = NULL;
  *tail = NULL;
  *(tail + 1) = (void *) head;
  nextlinkitem = *head;
  curpos = 1;
  linkitems = 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// move()    Causes 'nextlinkitem' to traverse the specified number of nodes,//
//           updates 'curpos' to be the node to which 'nextlinkitem' points. //
//                                                                           //
// 'numberofnodes' is a number indicating how many nodes need be traversed   //
// (not counter the current node) need be traversed. It may be positive(move //
// forward) or negative (move backward).  Return TRUE if it is successful.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::link::move(int numberofnodes)
{
  void **nownode;
  int i;

  nownode = (void **) nextlinkitem;
  if (numberofnodes > 0) {
    // Move forward.
    i = 0;
    while ((i < numberofnodes) && *nownode) {
      nownode = (void **) *nownode;
      i++;
    }
    if (*nownode == NULL) return false;
    nextlinkitem = (void *) nownode;
    curpos += numberofnodes;
  } else if (numberofnodes < 0) {
    // Move backward.
    i = 0;
    numberofnodes = -numberofnodes;
    while ((i < numberofnodes) && *(nownode + 1)) {
      nownode = (void **) *(nownode + 1);
      i++;
    }
    if (*(nownode + 1) == NULL) return false;
    nextlinkitem = (void *) nownode;
    curpos -= numberofnodes;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// locate()    Locates the node at the specified position.                   //
//                                                                           //
// The number 'pos' (between 1 and 'linkitems') indicates the location. This //
// routine first decides the shortest path traversing from 'curpos' to 'pos',//
// i.e., from head, tail or 'curpos'.   Routine 'move()' is called to really //
// traverse the link. If success, 'nextlinkitem' points to the node, 'curpos'//
// and 'pos' are equal. Otherwise, return FALSE.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::link::locate(int pos)
{
  int headdist, taildist, curdist;
  int abscurdist, mindist;

  if (pos < 1 || pos > linkitems) return false;

  headdist = pos - 1;
  taildist = linkitems - pos;
  curdist = pos - curpos;
  abscurdist = curdist >= 0 ? curdist : -curdist;

  if (headdist > taildist) {
    if (taildist > abscurdist) {
      mindist = curdist;
    } else {
      // taildist <= abs(curdist)
      mindist = -taildist;
      goend();
    }
  } else {
    // headdist <= taildist
    if (headdist > abscurdist) {
      mindist = curdist;
    } else {
      // headdist <= abs(curdist)
      mindist = headdist;
      rewind();
    }
  }

  return move(mindist);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// add()    Add a node at the end of this link.                              //
//                                                                           //
// A new node is appended to the end of the link.  If 'newitem' is not NULL, //
// its conents will be copied to the data slot of the new node. Returns the  //
// pointer to the newest added node.                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::link::add(void* newitem)
{
  void **newnode = tail;
  if (newitem != (void *) NULL) {
    memcpy((void *)(newnode + 2), newitem, linkitembytes);
  }
  tail = (void **) alloc();
  *tail = NULL;
  *newnode = (void*) tail;
  *(tail + 1) = (void*) newnode;
  linkitems++;
  return (void *)(newnode + 2);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insert()    Inserts a node before the specified position.                 //
//                                                                           //
// 'pos' (between 1 and 'linkitems') indicates the inserting position.  This //
// routine inserts a new node before the node of 'pos'.  If 'newitem' is not //
// NULL,  its conents will be copied into the data slot of the new node.  If //
// 'pos' is larger than 'linkitems', it is equal as 'add()'.  A pointer to   //
// the newest inserted item is returned.                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::link::insert(int pos, void* insitem)
{
  if (!locate(pos)) {
    return add(insitem);
  }

  void **nownode = (void **) nextlinkitem;

  // Insert a node before 'nownode'.
  void **newnode = (void **) alloc();
  if (insitem != (void *) NULL) {
    memcpy((void *)(newnode + 2), insitem, linkitembytes);
  }

  *(void **)(*(nownode + 1)) = (void *) newnode;
  *newnode = (void *) nownode;
  *(newnode + 1) = *(nownode + 1);
  *(nownode + 1) = (void *) newnode;

  linkitems++;

  nextlinkitem = (void *) newnode;
  return (void *)(newnode + 2);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// del()    Delete a node.                                                   //
//                                                                           //
// Returns a pointer of the deleted data. If you try to delete a non-existed //
// node (e.g. link is empty or a wrong index is given) return NULL.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::link::deletenode(void** deadnode)
{
  void **nextnode = (void **) *deadnode;
  void **prevnode = (void **) *(deadnode + 1);
  *prevnode = (void *) nextnode;
  *(nextnode + 1) = (void *) prevnode;

  dealloc((void *) deadnode);
  linkitems--;

  nextlinkitem = (void *) nextnode;
  return (void *)(deadnode + 2);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// del()    Delete a node at the specified position.                         //
//                                                                           //
// 'pos' between 1 and 'linkitems'.  Returns a pointer of the deleted data.  //
// If you try to delete a non-existed node (e.g. link is empty or a wrong    //
// index is given) return NULL.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::link::del(int pos)
{
  if (!locate(pos) || (linkitems == 0)) {
    return (void *) NULL;
  }
  return deletenode((void **) nextlinkitem);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// getitem()    The link traversal routine.                                  //
//                                                                           //
// Returns the node to which 'nextlinkitem' points. Returns a 'NULL' if the  //
// end of the link is reaching.  Both 'nextlinkitem' and 'curpos' will be    //
// updated after this operation.                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::link::getitem()
{
  if (nextlinkitem == (void *) tail) return NULL;
  void **nownode = (void **) nextlinkitem;
  nextlinkitem = *nownode;
  curpos += 1;
  return (void *)(nownode + 2);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// getnitem()    Returns the node at a specified position.                   //
//                                                                           //
// 'pos' between 1 and 'linkitems'. After this operation, 'nextlinkitem' and //
// 'curpos' will be updated to indicate this node.                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::link::getnitem(int pos)
{
  if (!locate(pos)) return NULL;
  return (void *)((void **) nextlinkitem + 2);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initializepools()    Initialize pools of mesh elements.                   //
//                                                                           //
// The sizes of the tetrahedron, shellface, and point will be calculated.    //
// Some class variables, such as 'pointmarkindex', 'elemmarkindex', 'volume- //
// boundindex', 'dummypoint', etc, are initialized. The pools of tetrahedra, //
// points, subfaces, and segments are allocated.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::initializepools()
{
  enum wordtype wtype;
  int ptsize, elesize;

  // A point contains 3 coordinates, 1 weight, plus 'n' attributes in REALs,
  //   and other fields, such as pointers, boundary markers, etc. The total
  //   size (in byte) of a point is calcualted below.
  ptsize = (4 + in->numberofpointattributes) * sizeof(REAL);
  // The index within each point at which an element pointer is found, where
  //   the index is measured in pointers. Ensure the index is aligned to a
  //   sizeof(tetrahedron)-byte address.
  point2tetindex = (ptsize + sizeof(tetrahedron) - 1) / sizeof(tetrahedron);
  // Increase the point size by two pointers, which are:
  //   - a pointer to a tet, read by point2tet();
  //   - a pointer to a parent point, read by point2ppt().
  ptsize = (point2tetindex + 2) * sizeof(tetrahedron);
  // The index within each point at which the boundary marker is found,
  //   Ensure the marker is aligned to a sizeof(int)-byte address.
  pointmarkindex = (ptsize + sizeof(int) - 1) / sizeof(int);
  // Increase the point size by two integers, which are:
  //   - an integer for boundary marker;
  //   - an integer for vertex type;
  ptsize = (pointmarkindex + 2) * sizeof(int);
  // Decide the wordtype used in point pool.
  wtype = (sizeof(REAL) >= sizeof(tetrahedron)) ? FLOATINGPOINT : POINTER;
  // Initialize the pool of vertices.
  pointpool = new memorypool(ptsize, VERPERBLOCK, wtype, 0);
  
  // Initialize spaces for 'dummypoint'.
  dummypoint = new REAL[(ptsize + sizeof(REAL) - 1) / sizeof(REAL)];
  dummypoint[0] = dummypoint[1] = dummypoint[2] = dummypoint[3] = 0.0;
  pointmark(dummypoint) = -1;

  // The number of bytes occupied by a tetrahedron.  There are 4 pointers
  //   to other tetrahedra, 4 pointers to corners, and possibly 4 pointers
  //   to subfaces.
  elesize = (8 + b->useshelles * 4) * sizeof(tetrahedron);
  // The index within each element at which its attributes are found, where
  //   the index is measured in REALs. 
  elemattribindex = (elesize + sizeof(REAL) - 1) / sizeof(REAL);
  // The index within each element at which the maximum voulme bound is
  //   found, where the index is measured in REALs.  Note that if the
  //   `b->regionattrib' flag is set, an additional attribute will be added.
  volumeboundindex = elemattribindex + in->numberoftetrahedronattributes
                   + (b->regionattrib > 0);
  // If element attributes or an constraint are needed, increase the number
  //   of bytes occupied by an element.
  if (b->varvolume) {
    elesize = (volumeboundindex + 1) * sizeof(REAL);
  }
  // The index within each element at which its marker is found, where the
  //   index is measured in ints.
  elemmarkerindex = (elesize + sizeof(int) - 1) / sizeof(int);
  // Increase the size by one interger.
  elesize = (elemmarkerindex + 1) * sizeof(int);
  // If -o2 switch is used, an additional pointer pointed to the list of
  //   higher order nodes is allocated for each element.
  highorderindex = (elesize + sizeof(tetrahedron) - 1) / sizeof(tetrahedron);
  if (b->order == 2) {
    elesize = (highorderindex + 1) * sizeof(tetrahedron);
  }
  // Having determined the memory size of an element, initialize the pools.
  tetrahedronpool = new memorypool(elesize, ELEPERBLOCK, POINTER, 16);
  hulltetrahedronpool = new memorypool(elesize, ELEPERBLOCK, POINTER, 16);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedrondealloc()    Deallocate space for a tet, marking it dead.      //
//                                                                           //
// Set the first vertex of 'dyingtet' to NULL. So we can detect dead tets    //
// when when traversing the list of all tetrahedra.                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::tetrahedrondealloc(memorypool *pool, tetrahedron *dyingtet)
{
  dyingtet[4] = (tetrahedron) NULL;
  pool->dealloc((void *) dyingtet);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedrontraverse()    Traverse the tetrahedra, skipping dead ones.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::tetrahedron* tetgenmesh::tetrahedrontraverse(memorypool *pool)
{
  tetrahedron *newtetrahedron;

  do {
    newtetrahedron = (tetrahedron *) pool->traverse();
    if (newtetrahedron == (tetrahedron *) NULL) {
      return (tetrahedron *) NULL;
    }
  } while (newtetrahedron[4] == (tetrahedron) NULL); // Skip dead ones.
  return newtetrahedron;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// shellfacedealloc()    Deallocate space for a shellface, marking it dead.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::shellfacedealloc(memorypool *pool, shellface *dyingsh)
{
  dyingsh[3] = (shellface) NULL;
  pool->dealloc((void *) dyingsh);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// shellfacetraverse()    Traverse the subfaces, skipping dead ones.         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::shellface* tetgenmesh::shellfacetraverse(memorypool *pool)
{
  shellface *newshellface;

  do {
    newshellface = (shellface *) pool->traverse();
    if (newshellface == (shellface *) NULL) {
      return (shellface *) NULL;
    }
  } while (newshellface[3] == (shellface) NULL); // Skip dead ones.
  return newshellface;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// badfacedealloc()    Deallocate space for a badface, marking it dead.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::badfacedealloc(memorypool *pool, badface *dying)
{
  dying->forg = (point) NULL;
  pool->dealloc((void *) dying);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// badfacetraverse()    Traverse the pools, skipping dead ones.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::badface* tetgenmesh::badfacetraverse(memorypool *pool)
{
  badface *newsh;

  do {
    newsh = (badface *) pool->traverse();
    if (newsh == (badface *) NULL) {
      return (badface *) NULL;
    }
  } while (newsh->forg == (point) NULL); // Skip dead ones.
  return newsh;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// pointdealloc()    Deallocate space for a point, marking it dead.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::pointdealloc(point dyingpoint)
{
  pointtype(dyingpoint) = DEADVERTEX;
  pointpool->dealloc((void *) dyingpoint);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// pointtraverse()    Traverse the points, skipping dead ones.               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::point tetgenmesh::pointtraverse()
{
  point newpoint;

  do {
    newpoint = (point) pointpool->traverse();
    if (newpoint == (point) NULL) {
      return (point) NULL;
    }
  } while (pointtype(newpoint) == DEADVERTEX); // Skip dead ones.
  return newpoint;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// maketetrahedron()    Create a new tetrahedron.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::maketetrahedron(memorypool* pool, triface *newtet)
{
  int i;

  newtet->tet = (tetrahedron *) pool->alloc();
  for (i = 0; i < 12; i++) {
    newtet->tet[i] = (tetrahedron) NULL;
  }
  for (i = 0; i < in->numberoftetrahedronattributes; i++) {
    elemattribute(newtet->tet, i) = 0.0;
  }
  if (b->varvolume) {
    volumebound(newtet->tet) = -1.0;
  }
  elemmarker(newtet->tet) = 0;
  newtet->loc = 0;
  newtet->ver = 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// makepoint()    Create a new point.                                        //
//                                                                           //
// The new point is indexed (starting from 'in->firstnumber'). It's type is  //
// initialized as UNUSEDVERTEX.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::makepoint(point* pnewpoint)
{
  int ptmark, i;

  *pnewpoint = (point) pointpool->alloc();
  // Initialize the list of coordinates and user-defined attributes.
  for (i = 0; i < 4 + in->numberofpointattributes; i++) {
    (*pnewpoint)[i] = 0.0;
  }
  // Initialize the point marker (starting from in->firstnumber).
  ptmark = (int) pointpool->items - (in->firstnumber == 1 ? 0 : 1);
  pointmark(*pnewpoint) = ptmark;
  pointtype(*pnewpoint) = UNUSEDVERTEX;
}

#endif // #ifndef memorypoolCXX