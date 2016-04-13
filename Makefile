.SUFFIXES: .f90 .f .o .mod .c .cpp

#------------------------------
# Define the compiler to use
#------------------------------

CC = mpicc
CX = mpicxx
FC = mpif90

#------------------------------
# Define any compile-time flags
#------------------------------

CC_FLAGS =  -g -Wall
CX_FLAGS =  -g -Wall
FC_FLAGS =  -g -cpp -dM  -Wno-unused -fbounds-check -Wall


#------------------------------
# Define the suffixes in use
#------------------------------

SRC_DIR=src
OBJ_DIR=obj
BIN_DIR=bin

#-----------------------------------------------------------------------
# Define any directories containing header files other than /usr/include
#-----------------------------------------------------------------------

INCLUDES = -I/usr/local/include -I${PETSC_DIR}/${PETSC_ARCH}/include \
	-I${PETSC_DIR}/include #-I${SLEPC_DIR}/include

#export PETSC_DIR=/home/balay/petsc-3.2-p0
#export PETSC_ARCH=linux-gnu-c-debug      
  
#-----------------------------------------------------------------------
# Define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
#-----------------------------------------------------------------------

LIB_FLAGS = -L./lib -L${LD_LIBRARY_PATH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib #-L../lib #-L/home/newhall/lib  -L../lib

#-----------------------------------------------------------------------
# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
#-----------------------------------------------------------------------

LIBS =  -ldl -lm -lstdc++ -llapack -lblas -lpetsc

SRC  :=	 src/utils.f90 src/main.f90

#-----------------------------------------------------------------------
# define the C,C++, Fortran object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .f90,.c,.cpp of all words in the macro SRCS
# with the .o suffix
#------------------------------------------------------------------------

OBJ = $(patsubst src/%.f90,obj/%.o,$(SRC))

#------------------------------
# define the executable file 
#------------------------------

TARGET = $(BIN_DIR)/test

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

all:    $(TARGET)
	@echo "\nCompilation and linking success...\n"

$(TARGET): $(OBJ)
	$(FC) $(FC_FLAGS) $(INCLUDES) -o $(TARGET) $(OBJ) $(LIB_FLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)

# Compile
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FC_FLAGS) -c  $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CC_FLAGS) -c  $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CX) $(CX_FLAGS) -c  $< -o $@

%.o : %.mod

clean:
	$(RM) $(SRC_DIR)/*~ ${OBJ_DIR}/*.o $(TARGET) *.mod

######################################
#          Compilation
######################################

#$(OBJ_DIR)/.$(SU_FC).o:
#	$(FC) $(FC_FLAGS) -c  $< -o $@

#.$(SU_CC).o:
#	$(CC) $(CC_FLAGS) -c  $< -o $@

#.$(SU_CX).o:
#	$(CX) $(CX_FLAGS) -c  $< -o $@
