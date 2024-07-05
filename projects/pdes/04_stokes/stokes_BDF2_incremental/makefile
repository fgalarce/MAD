include ../../makefile_header

compile: indicatrix

indicatrix: indicatrix.cpp 
	$(MAKE) -o indicatrix.exe indicatrix.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./indicatrix.exe
	install_name_tool -change libassimilation.dylib $(MAD_ROOT)/assimilation/libassimilation.dylib ./indicatrix.exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./indicatrix.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./indicatrix.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./indicatrix.exe

clean:
	rm *.exe
