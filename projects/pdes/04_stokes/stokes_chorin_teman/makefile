include ../../makefile_header

compile: frac_step

frac_step: frac_step.cpp 
	$(MAKE) -o frac_step.exe frac_step.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi_cxx -headerpad_max_install_names -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools
	install_name_tool -change libcfd.dylib /Users/fgalarce/research/4d-flow-ultrasound/cfd/libcfd.dylib ./frac_step.exe
	install_name_tool -change libassimilation.dylib /Users/fgalarce/research/4d-flow-ultrasound/assimilation/libassimilation.dylib ./frac_step.exe
	install_name_tool -change libfem.dylib /Users/fgalarce/research/4d-flow-ultrasound/fem/libfem.dylib ./frac_step.exe
	install_name_tool -change libio.dylib /Users/fgalarce/research/4d-flow-ultrasound/io/libio.dylib ./frac_step.exe
	install_name_tool -change libtools.dylib /Users/fgalarce/research/4d-flow-ultrasound/tools/libtools.dylib ./frac_step.exe

clean:
	rm *.exe
