include ../../makefile_header

imrp: imrp.cpp 
	$(MAKE) -o imrp.exe imrp.cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi_cxx -headerpad_max_install_names -L$(CFD) -lcfd -L$(CORE) -lcore -L$(FEM) -lfem -L$(IO) -lio -L$(TOOLS) -ltools
	install_name_tool -change libcfd.dylib /Users/fgalarce/research/4d-flow-ultrasound/cfd/libcfd.dylib ./imrp.exe
	install_name_tool -change libcore.dylib /Users/fgalarce/research/4d-flow-ultrasound/core/libcore.dylib ./imrp.exe
	install_name_tool -change libfem.dylib /Users/fgalarce/research/4d-flow-ultrasound/fem/libfem.dylib ./imrp.exe
	install_name_tool -change libio.dylib /Users/fgalarce/research/4d-flow-ultrasound/io/libio.dylib ./imrp.exe
	install_name_tool -change libtools.dylib /Users/fgalarce/research/4d-flow-ultrasound/tools/libtools.dylib ./imrp.exe

clean:
	rm *.exe
