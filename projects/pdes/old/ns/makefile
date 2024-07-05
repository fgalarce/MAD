include ../../makefile_header

exe: ns

ns: ns.cpp 
	$(MAKE) -o ns.exe ns.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools -L$(LINEAR_ALGEBRA) -llinear_algebra
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./ns.exe
	install_name_tool -change libassimilation.dylib $(MAD_ROOT)/assimilation/libassimilation.dylib ./ns.exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./ns.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./ns.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./ns.exe
	install_name_tool -change liblinear_algebra.dylib $(MAD_ROOT)/la/liblinear_algebra.dylib ./ns.exe
#
clean:
	rm *.exe
