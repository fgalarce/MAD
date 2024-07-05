include ../../makefile_header

exe: cMDS

cMDS: cMDS.cpp 
	$(MAKE) -o cMDS.exe cMDS.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(CORE) -lcore -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools -L$(LINEAR_ALGEBRA) -llinear_algebra
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./cMDS.exe
	install_name_tool -change libcore.dylib $(MAD_ROOT)/core/libcore.dylib ./cMDS.exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./cMDS.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./cMDS.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./cMDS.exe
	install_name_tool -change liblinear_algebra.dylib $(MAD_ROOT)/la/liblinear_algebra.dylib ./cMDS.exe

clean:
	rm *.exe
