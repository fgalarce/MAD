include ../../makefile_header

exe: MDS

MDS: MDS.cpp 
	$(MAKE) -o MDS.exe MDS.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./MDS.exe
	install_name_tool -change libassimilation.dylib $(MAD_ROOT)/assimilation/libassimilation.dylib ./MDS.exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./MDS.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./MDS.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./MDS.exe
	install_name_tool -change liblinear_algebra.dylib $(MAD_ROOT)/la/liblinear_algebra.dylib ./MDS.exe

clean:
	rm *.exe
