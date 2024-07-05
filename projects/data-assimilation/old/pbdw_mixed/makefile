include ../../makefile_header

compile_exe: least_squares.cpp 
	$(MAKE) -o least_squares.exe least_squares.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi_cxx -headerpad_max_install_names -L$(CFD) -lcfd -L$(CORE) -lcore -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools
	install_name_tool -change libcfd.dylib /Users/fgalarce/research/4d-flow-ultrasound/cfd/libcfd.dylib ./least_squares.exe
	install_name_tool -change libcore.dylib /Users/fgalarce/research/4d-flow-ultrasound/core/libcore.dylib ./least_squares.exe
	install_name_tool -change libfem.dylib /Users/fgalarce/research/4d-flow-ultrasound/fem/libfem.dylib ./least_squares.exe
	install_name_tool -change libio.dylib /Users/fgalarce/research/4d-flow-ultrasound/io/libio.dylib ./least_squares.exe
	install_name_tool -change libtools.dylib /Users/fgalarce/research/4d-flow-ultrasound/tools/libtools.dylib ./least_squares.exe

clean:
	rm *.exe
