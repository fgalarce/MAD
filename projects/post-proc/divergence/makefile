include ../../makefile_header

divergence: divergence.cpp 
	$(MAKE) -o divergence.exe divergence.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi_cxx -headerpad_max_install_names -L$(CFD) -lcfd -L$(CORE) -lcore -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools -L$(LINEAR_ALGEBRA) -llinear_algebra
	install_name_tool -change libcfd.dylib /Users/fgalarce/research/4d-flow-ultrasound/cfd/libcfd.dylib ./divergence.exe
	install_name_tool -change libcore.dylib /Users/fgalarce/research/4d-flow-ultrasound/core/libcore.dylib ./divergence.exe
	install_name_tool -change liblinear_algebra.dylib /Users/fgalarce/research/4d-flow-ultrasound/la/liblinear_algebra.dylib ./divergence.exe
	install_name_tool -change libfem.dylib /Users/fgalarce/research/4d-flow-ultrasound/fem/libfem.dylib ./divergence.exe
	install_name_tool -change libio.dylib /Users/fgalarce/research/4d-flow-ultrasound/io/libio.dylib ./divergence.exe
	install_name_tool -change libtools.dylib /Users/fgalarce/research/4d-flow-ultrasound/tools/libtools.dylib ./divergence.exe
clean:
	rm *.exe
