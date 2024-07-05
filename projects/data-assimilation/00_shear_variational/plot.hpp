/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2021, 2023
    
     Felipe Galarce at INRIA/WIAS/PUCV

  MDMA is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MDMA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <mad.hpp>

class Ploter {

  public:

    Ploter(){}
    ~Ploter(){
      fPostProc.close();
    };

    void initialize(Boundary & boundary, Geometry & geometry, Parameters parameters){
      par = parameters;
      geo = geometry;
      bd = boundary;
      fPostProc.open(par.dirResults() + "/norm_solution.txt");
      m_flow.resize(geo.bdLabels().size());
    }

    void postProcess(Vec u){
      m_normSol = norm(u);
      for (int i = 0; i < geo.bdLabels().size(); i++){
        m_flow[i]= bd.flow(u, geo.bdLabels()[i]);
      }
    }

    void plot(Vec u){
      postProcess(u);
      fPostProc << m_normSol  << " ";
      for (int i = 0; i < geo.bdLabels().size(); i++){
        fPostProc << m_flow[i] << " ";
      }
      fPostProc << endl;
      MADconsole("python3 plot.py");
      if (!pdf_is_open){
        MADconsole("evince --presentation norm.pdf &");
        pdf_is_open = true;
      }
    }

  private:
    bool pdf_is_open = false;
    Boundary bd;
    Geometry geo;
    Parameters par;
    ofstream fPostProc; 
    double m_normSol;
    vector<double> m_flow;

};
