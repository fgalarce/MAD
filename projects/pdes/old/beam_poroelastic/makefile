include ../../../makefile_header

exe: poro

poro: poro.cpp 
	$(MAKE) -o poro.exe poro.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools -L$(LINEAR_ALGEBRA) -llinear_algebra
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./poro.exe
	install_name_tool -change libassimilation.dylib $(MAD_ROOT)/assimilation/libassimilation.dylib ./poro.exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./poro.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./poro.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./poro.exe
	install_name_tool -change liblinear_algebra.dylib $(MAD_ROOT)/la/liblinear_algebra.dylib ./poro.exe
#
clean:
	rm *.exe
