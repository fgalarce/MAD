include ../../makefile_header

compile_exe: error

error: error.cpp 
	$(MAKE) -o error.exe error.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi_cxx -headerpad_max_install_names -L$(CFD) -lcfd -L$(CORE) -lcore -L$(FEM) -lfem -L$(LINEAR_ALGEBRA) -llinear_algebra -L$(IO) -lio -L$(TOOLS) -ltools -L$(ASSIMILATION) -lassimilation
	install_name_tool -change libcfd.dylib /Users/fgalarce/research/4d-flow-ultrasound/cfd/libcfd.dylib ./error.exe
	install_name_tool -change libassimilation.dylib /Users/fgalarce/research/4d-flow-ultrasound/assimilation/libassimilation.dylib ./error.exe
	install_name_tool -change liblinear_algebra.dylib /Users/fgalarce/research/4d-flow-ultrasound/la/liblinear_algebra.dylib ./error.exe
	install_name_tool -change libfem.dylib /Users/fgalarce/research/4d-flow-ultrasound/fem/libfem.dylib ./error.exe
	install_name_tool -change libio.dylib /Users/fgalarce/research/4d-flow-ultrasound/io/libio.dylib ./error.exe
	install_name_tool -change libtools.dylib /Users/fgalarce/research/4d-flow-ultrasound/tools/libtools.dylib ./error.exe
clean:
	rm *.exe
