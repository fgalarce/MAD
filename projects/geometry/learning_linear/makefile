include ../../makefile_header

exe: learn

learn: learn.cpp 
	$(MAKE) -o learn.exe learn.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools -L$(LINEAR_ALGEBRA) -llinear_algebra
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./learn.exe
	install_name_tool -change libassimilation.dylib $(MAD_ROOT)/assimilation/libassimilation.dylib ./learn.exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./learn.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./learn.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./learn.exe
	install_name_tool -change liblinear_algebra.dylib $(MAD_ROOT)/la/liblinear_algebra.dylib ./learn.exe
#
clean:
	rm *.exe
