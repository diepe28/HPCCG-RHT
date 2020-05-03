#@HEADER
# ************************************************************************
#
#               HPCCG: Simple Conjugate Gradient Benchmark Code
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
#@HEADER


# Simple hand-tuned makefile.  Modify as necessary for your environment.
# Questions? Contact Mike Heroux (maherou@sandia.gov).
#

#
# 0) Specify compiler and linker:

#CXX=/usr/bin/g++
#LINKER=/usr/bin/g++
CXX=mpicxx
LINKER=mpicxx


# 1) Build with MPI or not?
#    If you want to run the program with MPI, make sure USE_MPI is set
#    to -DUSING_MPI

USE_MPI =
#USE_MPI = -DUSING_MPI


# 2) MPI headers:
#    If you:
#    - Are building MPI mode (-DUSING_MPI is set above).
#    - Do not have the MPI headers installed in a default search directory and
#    - Are not using MPI compiler wrappers
#    Then specify the path to your MPI header file (include a -I)

MPI_INC = -I/usr/MPICH/SDK.gcc/include


# 3) Specify C++ compiler optimization flags (if any)
#    Typically some reasonably high level of optimization should be used to
#    enhance performance.

#IA32 with GCC:
#CPP_OPT_FLAGS = -O3 -funroll-all-loops -malign-double
#CPP_OPT_FLAGS = -O3 -ftree-vectorize -ftree-vectorizer-verbose=2
#CPP_OPT_FLAGS = -O3 -std=c++11
# No compiler optimizations
CPP_OPT_FLAGS = -std=c++11 -g
#
# 4) MPI library:
#    If you:
#    - Are building MPI mode (-DUSING_MPI is set above).
#    - Do not have the MPI library installed a default search directory and
#    - Are not using MPI compiler wrappers for linking
#    Then specify the path to your MPI library (include -L and -l directives)

MPI_LIB = -L/usr/MPICH/SDK.gcc/lib -lmpich

#
# 5) Build with OpenMP or not?
#    If you want to run the program with OpenMP, make sure USING_OMP is set
#    to -DUSING_OMP

USE_OMP =
#USE_OMP = -DUSING_OMP

#
# 6) OpenMP Compiler argument
#    GCC and Intel compilers require -fopenmp and -openmp, resp.  Other compilers may differ.

#OMP_FLAGS = -fopenmp
#OMP_FLAGS = -openmp

#
# 7) System libraries: (May need to add -lg2c before -lm)

SYS_LIB =-lm -lpthread

#
# 6) Specify name if executable (optional):
TARGET = Wang

# other compilation flags
#COMP_FLAGS = -DAPPROACH_WANG=1 -DPRINT_OUTPUT=1

COMP_FLAGS = -DAPPROACH_WANG=1
#COMP_FLAGS = -DAPPROACH_WANG=1 -DVAR_GROUPING=1
#COMP_FLAGS = -DAPPROACH_WANG=1 -DJUST_VOLATILES=1

################### Derived Quantities (no modification required) ##############

CXXFLAGS= $(CPP_OPT_FLAGS) $(USE_OMP) $(USE_MPI) $(MPI_INC) $(COMP_FLAGS)
LIB_PATHS= $(SYS_LIB)

# Every *.cpp file
TEST_CPP=$(wildcard *.cpp)

# Every value of TEST_CPP change it from .cpp -> .o
TEST_OBJ=$(TEST_CPP:.cpp=.o)

#	$(LINKER) $(CPP_OPT_FLAGS) $(TEST_OBJ) $(LIB_PATHS) -o $(TARGET)

#flipIp variables
CFLAGS = -g -I$(FLIPIT_PATH)/include $(CXXFLAGS)
FILIB  = -L$(FLIPIT_PATH)/lib -lcorrupt
FIPASS = $(FLIPIT_PATH)/lib/libFlipItPass.so
LFLAGS = $(FILIB)
flipit-cxx = $(FLIPIT_PATH)/scripts/flipit-c++ $(CFLAGS)
#CC=$(flipit-cxx) #this one is needed for flipit-c++ to work
#OMPI_CXX=clang++ mpicxx -show:command

#$@, the name of the TARGET
#$<, the name of the first prerequisite

all : Wang WangVG WangJV

Wang: $(TEST_OBJ)
	$(LINKER) -o $(TARGET) $(TEST_OBJ) $(LFLAGS) $(SYS_LIB)

WangVG : COMP_FLAGS += -DVAR_GROUPING=1
WangVG: $(TEST_OBJ)
	$(LINKER) -o WangVG $(TEST_OBJ) $(LFLAGS) $(SYS_LIB)

WangJV : COMP_FLAGS += -DJUST_VOLATILES=1
WangJV: $(TEST_OBJ)
	$(LINKER) -o WangJV $(TEST_OBJ) $(LFLAGS) $(SYS_LIB)

#%.o: %.c
#	$(flipit-cxx) -o $@ -c $<

main.o: main.cpp
	$(flipit-cxx) -o $@ -c $<

RHT.o: RHT.cpp
	$(flipit-cxx) -o $@ -c $<

HPCCG.o: HPCCG.cpp
	$(flipit-cxx) -o $@ -c $<

QueueStressTest.o: QueueStressTest.cpp
	$(flipit-cxx) -o $@ -c $<

YAML_Doc.o: YAML_Doc.cpp
	$(flipit-cxx) -o $@ -c $<

YAML_Element.o: YAML_Element.cpp
	$(flipit-cxx) -o $@ -c $<

exchange_externals.o: exchange_externals.cpp
	$(flipit-cxx) -o $@ -c $<

make_local_matrix.o: make_local_matrix.cpp
	$(flipit-cxx) -o $@ -c $<

ddot.o: ddot.cpp
	$(flipit-cxx) -o $@ -c $<

waxpby.o: waxpby.cpp
	$(flipit-cxx) -o $@ -c $<

HPC_sparsemv.o: HPC_sparsemv.cpp
	$(flipit-cxx) -o $@ -c $<

dump_matlab_matrix.o: dump_matlab_matrix.cpp
	$(flipit-cxx) -o $@ -c $<

mytimer.o: mytimer.cpp
	$(flipit-cxx) -o $@ -c $<

compute_residual.o: compute_residual.cpp
	$(flipit-cxx) -o $@ -c $<

read_HPC_row.o: read_HPC_row.cpp
	$(flipit-cxx) -o $@ -c $<

generate_matrix.o: generate_matrix.cpp
	$(flipit-cxx) -o $@ -c $<

test:
	@echo "Not implemented yet..."

.PHONY : all

#clean:
#	@rm -f *.o *.bc *.bin  *~ $(TARGET) $(TARGET).exe

clean:
	@rm -f *.o *.bc *.bin *.pyc *.LLVM.txt  *~

cleanBinaries:
	@rm -f Wang* *.log *~
