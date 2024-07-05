#ifndef WSS_ASSEMBLE
#define WSS_ASSEMBLE

#include <mad.hpp>
//For some reason, I just CANNOT SPLIT THIS

using namespace std;
//Can assemble and return an specific matrix from the ns class
PetscErrorCode code;
enum terms {term_kin, term_visc, term_pres, LHSmat};
/*
double quadratic(Vec u, Mat K, Vec v) {
	Vec Ku = vec(nbDofs); 
	MatMult(K, v, Ku); 
	double E; 
	code = VecDot(Ku, u, &E); CHKERRQ(code);
	VecDestroy(&Ku);
	return E;
}*/

void load(string source, Vec& target, int t, int nbVariables, IO& io) {
	vector<string> to_load(nbVariables);
	to_load[0] = source + "/ux." + wildcard(t) + ".scl";
	to_load[1] = source + "/uy." + wildcard(t) + ".scl";
	if (nbVariables>3) {
		//3D case
		to_load[2] = source + "/uz." + wildcard(t) + ".scl";
		to_load[3] = source + "/p." + wildcard(t) + ".scl";
	} else {
		to_load[2] = source + "/p." + wildcard(t) + ".scl";
	}
	
	io.loadVector(target, to_load);
}

class wss_assemble{
public:
	wss_assemble(Parameters _par, Geometry _geo, Boundary _bd, Calculus _calculus) {
		par = _par;
		geo = _geo;
		bd = _bd;
		calculus = _calculus;
		nbDofs = geo.nbVertices*par.nbVariables();
		//u0_iter = vec(geo.nbVertices*par.nbVariables());
	}
	//Vec u0_iter;
	//enum terms {term_kin, term_visc, term_pres, LHSmat};
	//Parameters par, Geometry geo, Boundary bd, 
	Mat assemble_static(int TermToAssemble) {
		
		int nbDofs = geo.nbVertices*(geo.dimension()+1);
		//Declares a matrix
		Mat K = mat(nbDofs, nbDofs);
		//mat(K, nbDofs, nbDofs);
		//Declare and initialize a FEM object
		FEM fem;
		fem.initialize(par, geo, bd, K);
		
		Mat M = mat(nbDofs, nbDofs);
		//mat(M, nbDofs, nbDofs);
		FEM femMass;
		femMass.initialize(par, geo, bd, M);
		
		switch(TermToAssemble) {
			
			case term_kin: 
			{
				for (int partId = 0; partId < geo.elements().size(); partId++){
					cout << "Navier Stokes: Assembling discretization matrix. Part " << partId << endl;
					for (int feId = 0; feId < geo.elements()[partId].size(); feId++){
						vector<int> simplex = geo.elements()[partId][feId];
						femMass.setSimplex(simplex);
						if (feId % (geo.elements()[partId].size()) == 0){
							cout << "  Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
						}
						for (int var = 0; var < par.nbVariables()-1; var++) {
							femMass.u_dot_v_scalar(par.density(), var, var);
						}
						
					}
				}
				code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
				code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
				return M;
			}
			
			case term_visc:
			{
				//Set the matrix values
				for (int partId = 0; partId < geo.elements().size(); partId++){
					cout << "Navier Stokes: Assembling discretization matrix. Part " << partId << endl;
					for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */
						vector<int> simplex = geo.elements()[partId][feId];
						fem.setSimplex(simplex); //Saves the finite elements coordinates

						if (feId % (geo.elements()[partId].size()) == 0){
							cout << "  Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
						}

						for (int var = 0; var < par.nbVariables()-1; var++){
							/* grad u : grad v */
							for (int comp = 0; comp < geo.dimension(); comp++){
								fem.du_x_dv_y(par.viscosity(), comp, comp, var, var);
							}
						}
					}
				}
				
				code = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); CHKERR(code); 
				code = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY); CHKERR(code);
			
				return K;
			}
			
			case LHSmat:
			{
			Mat M = assemble_static(term_kin);
			Mat S = mat(nbDofs, nbDofs);
			//mat(S, nbDofs, nbDofs);
			Mat A = mat(nbDofs, nbDofs);
			FEM femStat;
			femStat.initialize(par, geo, bd, A);
			
			code = MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY); CHKERR(code); 
			code = MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY); CHKERR(code);
			code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
			code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);
			
			code = MatDuplicate(M, MAT_COPY_VALUES, &S); CHKERR(code);
			if (par.timeIntegration() == "BDF1"){
				code = MatScale(S, 1.0/par.timeStep());
			} else if (par.timeIntegration() == "BDF2"){
				code = MatScale(S, 3.0/2.0/par.timeStep());
			}
			code = MatAXPY(S, 1.0, A, DIFFERENT_NONZERO_PATTERN); CHKERR(code);

			double normLHS = norm(S);
			cout << "Navier Stokes: Norm LHS = " << normLHS << endl;

			return S;
			}
		}
	}
	
	Mat assemble_convec(Mat LHS_static, Vec u0_custom) {
		Vec u0 = zeros(nbDofs);
		VecAssemblyBegin(u0);
		VecAssemblyEnd(u0);
		Mat C = mat(nbDofs, nbDofs);
		FEM femConv;
		femConv.initialize(par, geo, bd, C);
		
		vector<Vec> u0_split;
		if (u0_custom == NULL){
			u0_split = calculus.split(u0);
		} else {
			u0_split = calculus.split(u0_custom);
		}
		vector<Vec> u0_seq(geo.dimension());
		for (int i = 0; i < geo.dimension(); i++){
			u0_seq[i] = getSequential(u0_split[i]);
		}
		code = MatZeroEntries(C); CHKERR(code);

		for (int partId = 0; partId < geo.elements().size(); partId++){
			cout << "Navier Stokes: Assembling convection matrix. Part " << partId << endl;
			for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */
				vector<int> simplex = geo.elements()[partId][feId];
				femConv.setSimplex(simplex);
				// Get solution at previous time step per element
				vector<vector<double>> u_element(geo.dimension());
				for (int comp = 0; comp < geo.dimension(); comp++){
					u_element[comp].resize(geo.dimension()+1);
					code = VecGetValues(u0_seq[comp], geo.dimension()+1, &simplex[0], &u_element[comp][0]); CHKERR(code); CHKERR(code);
				}
				double sigma = 1.0;
				double Ck = 30.0;
				double hk = femConv.feSize();
				if (par.timeIntegration() == "BDF2"){
					sigma = 2.0;
				}
				double tau = std::pow(std::pow(sigma/par.timeStep(), 2.0) + Ck*std::pow(par.viscosity(), 2.0)/std::pow(hk, 4.0), -0.5);

				// Assembling
				for (int var = 0; var < par.nbVariables()-1; var++){
					for (int comp1 = 0; comp1 < geo.dimension(); comp1++){
						/* Convection (newtonian, no stab.)*/
						femConv.a_du_x_v(par.density(), comp1, u_element[comp1], var, var);
					}
				}
			}
		}
	for (int i = 0; i < geo.dimension(); i++){
	VecDestroy(&u0_seq[i]);
	VecDestroy(&u0_split[i]);
	}
	code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERR(code);
	code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERR(code);

	/* assemble full LHS matrix */ 
	MatAXPY(C, 1.0, LHS_static, DIFFERENT_NONZERO_PATTERN);
	return C;
	}
	
	Mat assemble_all() {}
	
private:
	Parameters par;
	Geometry geo;
	Boundary bd;
	Calculus calculus;
	int nbDofs;
};  

#endif
