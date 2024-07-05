#include "wss_assemble.hpp"
#include <mad.hpp>
#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <masterElement.hpp>
#include <geometry.hpp>
#include <boundaries.hpp>
#include <calculus.hpp>
#include <fem.hpp>
#include <ns.hpp>
#include "wss_assemble.hpp"

wss_assemble::Mat assemble_static(NS ns, int TermToAssemble) {
	// Check if stabilization is required
	bool stabilize = par.cfd_stab() != "no_stab";

	for (int partId = 0; partId < geo.elements().size(); partId++){
		if (m_world_rank == 0) cout << "Navier Stokes: Assembling discretization matrix. Part " << partId << endl;
		for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

		vector<int> simplex = geo.elements()[partId][feId];

		ns.femStat.setSimplex(simplex);
		ns.femMass.copySimplex(femStat); 

		if (feId % (geo.elements()[partId].size() / m_verbose) == 0){
			if (m_world_rank == 0) cout << "  Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
		}

		for (int var = 0; var < par.nbVariables()-1; var++){
			switch(TermToAssemble)
			{
			case 0:
      		{
			/* u \cdot v */
			femMass.u_dot_v_scalar(par.density(), var, var);
			break;}
		
			case 1:
			{
			/* grad u : grad v */
			if (!m_non_newtonian){
				for (int comp = 0; comp < geo.dimension(); comp++){
				// Add laplacian terms
				if (comp == var){
					femStat.du_x_dv_y(2.0*par.viscosity(), comp, comp, var, var);
				} else {
					femStat.du_x_dv_y(par.viscosity(), comp, comp, var, var);
				}
				// Add remaining terms
				if (par.viscousTerm() == "symmetric" and comp != var){
					femStat.du_x_dv_y(par.viscosity(), comp, var, var, comp);
				}
			}
		  
			}
			break;}
			}
		}
		switch(TermToAssemble)
		{case 2:
			{
			/* -p div(v) */
			for (int comp = 0; comp < geo.dimension(); comp++){
			femStat.u_dv_y(-1.0, comp, comp, geo.dimension());
			}
	    break;}
		case 3:
			{
			/* div(u) q */
			for (int comp = 0; comp < geo.dimension(); comp++){
			femStat.du_x_v(1.0, comp, geo.dimension(), comp);
			}
		break;}
		}
		}
	}
	code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
	code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
	code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
	code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);
	//code = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); CHKERR(code); 
	//code = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY); CHKERR(code);
	//code = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERR(code); 
	//code = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERR(code);

	/* LHS static matrix */
	Mat S = mat(nbDofs, nbDofs);
	code = MatDuplicate(M, MAT_COPY_VALUES, &S); CHKERR(code);
	if (par.timeIntegration() == "BDF1"){
		code = MatScale(S, 1.0/par.timeStep());
	} else if (par.timeIntegration() == "BDF2"){
		code = MatScale(S, 3.0/2.0/par.timeStep());
	}
	code = MatAXPY(S, 1.0, A, DIFFERENT_NONZERO_PATTERN); CHKERR(code);

	double normLHS = norm(S);
	if (m_world_rank == 0) cout << "Navier Stokes: Norm LHS = " << normLHS << endl;

	return ns.K;
}

