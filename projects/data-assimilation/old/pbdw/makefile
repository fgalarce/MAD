include ../../makefile_header

pbdw: pbdw.cpp 
	$(MAKE) -o pbdw.exe pbdw.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(CORE) -lcore -L$(FEM) -lfem -L$(LINEAR_ALGEBRA) -llinear_algebra -L$(IO) -lio -L$(TOOLS) -ltools -L$(ASSIMILATION) -lassimilation
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./pbdw.exe
	install_name_tool -change libcore.dylib $(MAD_ROOT)/core/libcore.dylib ./pbdw.exe
	install_name_tool -change liblinear_algebra.dylib $(MAD_ROOT)/la/liblinear_algebra.dylib ./pbdw.exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./pbdw.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./pbdw.exe
	install_name_tool -change libassimilation.dylib $(MAD_ROOT)/assimilation/libassimilation.dylib ./pbdw.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./pbdw.exe

clean:
	rm *.exe
