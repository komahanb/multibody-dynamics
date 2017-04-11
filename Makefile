.SUFFIXES: .f90 .f .o .mod .c .cpp

#------------------------------
# Define the compiler to use
#------------------------------
CC = gcc
CX = g++
FC = gfortran-6

#------------------------------
# Define any compile-time flags
#------------------------------
CC_FLAGS =  #-g #-Wall
CX_FLAGS =  #-g #-Wall
FC_FLAGS =  -g -cpp -fbounds-check #-fbounds-check -ffree-form -Wall -cpp -dM -Wno-unused

TARGET = $(BIN_DIR)/test

default: $(OBJ)
	$(FC) $(FC_FLAGS) $(INCLUDES) -o $(TARGET) $(OBJ) $(LIB_FLAGS) $(LIBS)

complex: $(OBJ)
	$(FC) $(FC_FLAGS) -DUSE_COMPLEX $(INCLUDES) -o $(TARGET) $(OBJ) $(LIB_FLAGS) $(LIBS)

#------------------------------
# Define the suffixes in use
#------------------------------
SRC_DIR=src
OBJ_DIR=obj
BIN_DIR=bin

#-----------------------------------------------------------------------
# Define any directories containing header files other than /usr/include
#-----------------------------------------------------------------------
INCLUDES = -I/usr/local/include -I./src

#-----------------------------------------------------------------------
# Define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
#-----------------------------------------------------------------------
LIB_FLAGS = -L./lib -L${LD_LIBRARY_PATH}/lib

#-----------------------------------------------------------------------
# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
#-----------------------------------------------------------------------
LIBS = -lm -llapack

#---------------------------------------------------------------------#
# Define all the source files to compile here
#---------------------------------------------------------------------#
SRC  :=	 src/utils.f90 src/lapack.f90 src/linear_algebra.f90 \
	 src/function.f90 src/physics.f90 src/nonlinear.f90 src/integrator.f90 \
         src/adams_bashforth_moulton.f90 \
         src/newmark_beta_gamma.f90 \
         src/runge_kutta.f90 src/backward_difference.f90 \
         src/smd.f90 src/smd_functions.f90 \
         src/aero_elastic_oscillator.f90 src/oscillator_functions.f90 \
	 src/vanderpol.f90 \
         src/test_adjoint.f90
#	 src/dae.f90 \
#         src/main.f90
#	 src/utils.f90 src/rotation.f90 src/dynamics.f90 \

#-----------------------------------------------------------------------
# define the C,C++, Fortran object files 
#------------------------------------------------------------------------

OBJ = $(patsubst src/%.f90 src/interface/%.f90 src/integrator/%.f90 src/dynamics/%.f90 src/test/%.f90,obj/%.o,$(SRC))

#------------------------------
# Executable
#------------------------------

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FC_FLAGS) -c  $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CC_FLAGS) -c  $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CX) $(CX_FLAGS) -c  $< -o $@

%.o : %.mod

clean:
	$(RM) $(SRC_DIR)/*~ ${OBJ_DIR}/*.o $(TARGET) *.mod
