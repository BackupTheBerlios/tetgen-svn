###############################################################################
#                                                                             #
# makefile for TetGen                                                         #
#                                                                             #
# Type "make" to compile TetGen into an executable program (tetgen).          #
# Type "make tetlib" to compile TetGen into a library (libtet.a).             #
# Type "make clean" to delete all object (*.o) files.                         #
#                                                                             #
###############################################################################

# CXX should be set to the name of your favorite C++ compiler.
# ===========================================================

CXX = g++
#CXX = icpc
#CXX = CC

# CXXFLAGS is the level of optimiztion, default is -O. One should try
# -O2, -O3 ... to find the best optimization level.
# ===================================================================

CXXFLAGS = -g

# PREDCXXFLAGS is for compiling J. Shewchuk's predicates. It should
# always be equal to -O0 (no optimization). Otherwise, TetGen may not
# work properly.

PREDCXXFLAGS = -O0

# SWITCHES is a list of switches to compile TetGen.
# =================================================
#
# By default, TetGen uses double precision floating point numbers.  If you
#   prefer single precision, use the -DSINGLE switch. 
#
# The source code of TetGen includes a lot of assertions, which are mainly
#   used for catching bugs at that places.  These assertions somewhat slow
#   down the speed of TetGen.  They can be skipped by define the -DNDEBUG
#   switch.

SWITCHES = -Wall -DSELF_CHECK

# SWITCHES = -Wall -Wabi -Wctor-dtor-privacy \
#            -Woverloaded-virtual -Wno-pmf-conversions -Wsign-promo \
#            -Wsynth  -Wchar-subscripts -Wconversion -Wsign-compare \
#            -Wcomment  -Wimplicit -Wmissing-braces -Wparentheses \
#            -Wreturn-type -Wswitch -Wswitch-default \
#            -Wswitch-enum -Wtrigraphs -W -DSELF_CHECK

# RM should be set to the name of your favorite rm (file deletion program).

RM = /bin/rm

# The action starts here.

behavior.o:	behavior.cxx
	$(CXX) -I./ $(CXXFLAGS) -c behavior.cxx

io.o:	io.cxx
	$(CXX) -I./ $(CXXFLAGS) -c io.cxx

dbg.o:	dbg.cxx
	$(CXX) -I./ $(CXXFLAGS) -c dbg.cxx

memorypool.o:	memorypool.cxx
	$(CXX) -I./ $(CXXFLAGS) -c memorypool.cxx

geom.o:	geom.cxx
	$(CXX) -I./ $(CXXFLAGS) -c geom.cxx

flip.o:	flip.cxx
	$(CXX) -I./ $(CXXFLAGS) -c flip.cxx

delaunay.o:	delaunay.cxx
	$(CXX) -I./ $(CXXFLAGS) -c delaunay.cxx

surface.o:	surface.cxx
	$(CXX) -I./ $(CXXFLAGS) -c surface.cxx

constrain.o:	constrain.cxx
	$(CXX) -I./ $(CXXFLAGS) -c constrain.cxx

refine.o:	refine.cxx
	$(CXX) -I./ $(CXXFLAGS) -c refine.cxx

meshio.o:	meshio.cxx
	$(CXX) -I./ $(CXXFLAGS) -c meshio.cxx

meshstat.o:	meshstat.cxx
	$(CXX) -I./ $(CXXFLAGS) -c meshstat.cxx

predicates.o: predicates.cxx
	$(CXX) $(PREDCXXFLAGS) -c predicates.cxx

tetlib:	behavior.o io.o memorypool.o geom.o flip.o delaunay.o surface.o constrain.o refine.o meshio.o meshstat.o predicates.o
	ar r libtet.a behavior.o io.o memorypool.o geom.o flip.o delaunay.o surface.o constrain.o refine.o meshio.o meshstat.o predicates.o

tetgen:	main.cxx
	$(CXX) -I./ $(CXXFLAGS) -o tetgen main.cxx -L./ -ltet

clean:
	$(RM) *.o *.a *~ tetgen
