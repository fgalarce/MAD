#include <mad.hpp>
#include "wss_assemble.hpp"
#include "../../../tools/petscWrapper.cpp"

int main(int argc, char *argv[]){

	/*
	 * DEBE RECIBIR COMO INPUT LOS PERFILES DE VELOCIDAD Y PRESIÓN, Y Vstokes
	 * DEBE ENSAMBLAR K y M
	 * 
	 * Debería cambiarle el nombre a wss_assemble? 
	 * Por hacer:
	 * -UI: Graficar, guardar, detectar exe, detectar nucleos y permitir cambiarlos
	 * -wss: Calcular energias*/
	Parameters par = MADinitialize(argc, argv);

	/* Initialize MAD objects */
	IO io;
	io.initialize(par);
    
	Geometry geo;
	geo.initialize(par, io);

	Boundary boundary;
	boundary.initialize(par, geo);

	Calculus calculus;
	calculus.initialize(par, geo, boundary);

	NS ns;
	ns.initialize(par, geo, boundary);
	//int TermToAssemble = 1; //1 -> viscous
	/* Loads the class that assembles the static terms */
	wss_assemble wss(par, geo, boundary, calculus);
	/* Get the static matrices */
	Mat K = wss.assemble_static(term_visc);
	//Mat K = wss.assemble_static(par, geo, boundary, 1); //9582*9582
	MADprint("+++MATRIX K ASSEMBLED+++\n");
	Mat M = wss.assemble_static(term_kin);
	//Mat M = wss.assemble_static(par, geo, boundary, 0); //9582*9582
	Mat A = wss.assemble_static(LHSmat);
	MADprint("+++MATRIX M ASSEMBLED+++\n");

	int nbDofs = io.nbVertices()*par.nbVariables();
  
	vector<Vec> v(geo.dimension()+1); //Momentum Test function
	vector<Vec> u(geo.dimension()+1); //Velocity
	vector<Vec> u1(geo.dimension()+1);

	//Falta inventar un parmetro para esto
//	string dir_vstokes="./results_vstokes";
	string dir_vstokes=par.dirSyntheticField();

	ofstream ofile_v("./energy_v.txt");
	ofstream ofile_k("./energy_k.txt");
	ofstream ofile_c("./energy_c.txt");
	ofstream ofile_drag("./energy_drag.txt");
  
	for (int t : range(par.start(), par.end(), par.jump())){
		MADprint("\n- - - - - - - - - - - - - - - - - - - - - - - - -\n");
		MADprint(" Computing Energy: ", t);
		
		//Load vectors from file
		Vec u = vec(nbDofs); //Current vector
		Vec u1 = vec(nbDofs); // Next vector
		Vec vstokes = vec(nbDofs); //Test functions
		Vec& v = vstokes;

		load(par.dirSyntheticField(), u, t, par.nbVariables(), io);
		load(par.dirSyntheticField(), u1, t+1, par.nbVariables(), io);
		load(dir_vstokes, v, t, par.nbVariables(), io);
		
		//--Estimates Ev
		double Ev = quadraticForm(u, K, v); 
		cout << "Ev = " << Ev << endl;
		ofile_v << Ev << endl;
		
		//--Estimates Ek
		VecAXPY(u1, -1.0, u);
		VecScale(u1, 1.0/par.timeStep());
		double Ek = quadraticForm(u1, M, v);
		cout << "Ek = " << Ek << endl;
		ofile_k << Ek << endl;
		
		//--Estimates Econv
		Mat C = ns.assembleLHS(A, u);
		double Ec = quadraticForm(u, C, v);
		cout << "Ec = " << Ec << endl;
		ofile_c << Ec << endl;
		
		//--Estimates the drag
		double drag = Ev+Ek+Ec;
		cout << "drag = " << drag << endl;
		ofile_drag << drag << endl;
	}

	ofile_v.close();
	ofile_k.close();
	ofile_c.close();
	ofile_drag.close();

	MADfinalize(par);
}
