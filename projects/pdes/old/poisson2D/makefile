include ../makefile_header

poisson: poisson.cpp 
	$(MAKE) -o poisson.exe poisson.cpp -L$(PETSC_LIB) -lpetsc -L$(MPI_LIB) -lmpi_cxx -headerpad_max_install_names -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools 
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./poisson.exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./poisson.exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./poisson.exe
clean:
	rm *.exe
