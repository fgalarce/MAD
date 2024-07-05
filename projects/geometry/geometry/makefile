include ../../makefile_header

exe: interpolate

interpolate: interpolate.cpp 
	$(MAKE) -o interpolate.exe interpolate.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools -L$(LINEAR_ALGEBRA) -llinear_algebra
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./interpolate.exe
	install_name_tool -change libassimilation.dylib $(MAD_ROOT)/assimilation/libassimilation.dylib ./interpolate.exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./interpolate.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./interpolate.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./interpolate.exe
	install_name_tool -change liblinear_algebra.dylib $(MAD_ROOT)/la/liblinear_algebra.dylib ./interpolate.exe
#
clean:
	rm *.exe
