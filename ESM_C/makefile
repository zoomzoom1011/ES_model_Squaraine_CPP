src     = sub/esm_c.cpp
exes    = esm_c.exe
ICPC    = icpc
INC     = -I$(MKLROOT)/include
FLAGS	= -DMKL_ILP64 -DADD_
LIBS    = -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lgomp -lpthread -lm -ldl
LIBDIRS = -L/opt/intel/mkl/lib
INCDIRS = -I/opt/intel/mkl/include


all: ${exes}

${exes}: ${src}
	$(ICPC) $(src) -o $(exes) $(FLAGS) $(LIBDIRS) $(LIBS) $(INCDIRS)

clean:
	rm esm_c.exe
