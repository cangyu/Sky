FC=gfortran
# #########################################################
# For pgf90 also add the flag: -Mbackslash
FCFLAGS=-g -O3
# #########################################################
# #########################################################
# Use the correct flag for your compiler
# gfortran: -fdefault-real-8
# ifort and pgf90: -r8
XFOIL_REAL8_FLAG=-fdefault-real-8
# #########################################################
OBJ=vardef.o util.o xfoil_deps.o math_deps.o optimization.o edge_grid.o memory.o surface_util.o elliptic_surface_grid.o hyperbolic_surface_grid.o surface_grid.o menu.o main.o
EXE=construct2d
DIR=src

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) -o $(EXE) $(OBJ)

allclean: 
	rm -f *.mod *.o $(EXE)

clean: 
	rm -f *.mod *.o

vardef.o: $(DIR)/vardef.f90
	$(FC) -c $(FCFLAGS) $(DIR)/vardef.f90

util.o: $(DIR)/util.f90
	$(FC) -c $(FCFLAGS) $(DIR)/util.f90

xfoil_deps.o: $(DIR)/xfoil_deps.f
	$(FC) -c $(FCFLAGS) $(XFOIL_REAL8_FLAG) $(DIR)/xfoil_deps.f

math_deps.o: $(DIR)/math_deps.f90
	$(FC) -c $(FCFLAGS) $(DIR)/math_deps.f90

optimization.o: $(DIR)/optimization.f90
	$(FC) -c $(FCFLAGS) $(DIR)/optimization.f90

edge_grid.o: $(DIR)/edge_grid.f90
	$(FC) -c $(FCFLAGS) $(DIR)/edge_grid.f90

memory.o: $(DIR)/memory.f90
	$(FC) -c $(FCFLAGS) $(DIR)/memory.f90

surface_util.o: $(DIR)/surface_util.f90
	$(FC) -c $(FCFLAGS) $(DIR)/surface_util.f90

elliptic_surface_grid.o: $(DIR)/elliptic_surface_grid.f90
	$(FC) -c $(FCFLAGS) $(DIR)/elliptic_surface_grid.f90

hyperbolic_surface_grid.o: $(DIR)/hyperbolic_surface_grid.f90
	$(FC) -c $(FCFLAGS) $(DIR)/hyperbolic_surface_grid.f90

surface_grid.o: $(DIR)/surface_grid.f90
	$(FC) -c $(FCFLAGS) $(DIR)/surface_grid.f90

menu.o: $(DIR)/menu.f90
	$(FC) -c $(FCFLAGS) $(DIR)/menu.f90

main.o: $(DIR)/main.f90
	$(FC) -c $(FCFLAGS) $(DIR)/main.f90
