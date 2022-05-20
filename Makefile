#key parameters:
OPENMP=y
COMPILER=default
GPU=n
CMD := test_petsc
ALL:${CMD}
LIB :=
LD_LIB :=

###########################################
##### Machine specific paths

CMP := mpif90

export LD_LIBRARY_PATH=/home/jkumar/Softwares/petsc-3.8.4/arch-linux2-c-debug/lib

LIB += -I/home/jkumar/Softwares/petsc-3.8.4/include/ -I/home/jkumar/Softwares/petsc-3.8.4/arch-linux2-c-debug/include/ -L/home/jkumar/Softwares/petsc-3.8.4/lib/ -L/home/jkumar/Softwares/petsc-3.8.4/arch-linux2-c-debug/lib/

OPT += -lpetsc


###########################################

OBJ := test_petsc.o

$(CMD): $(OBJ)
	        $(CMP) -o $(CMD) $(OBJ) $(LIB) $(OPT)

test_petsc.o: test_petsc.F90
	$(CMP) -c test_petsc.F90 $(LIB) $(OPT)

clean:
	rm -f *.o *.mod core test_petsc
