/*------------------------------------------------------------------------\
|   catalyticFOAM                                                         |
|   http://www.catalyticfoam.polimi.it/                                   |
|                                                                         |
|   Copyright (C) 2013 by Matteo Maestri and Alberto Cuoci                |
|                                                                         |
|   Contributors (2010-2013)                                              |
|   Alberto Cuoci, Matteo Maestri, Tiziano Maffei                         |
|   Sandro Goisis, Alessandra Osio, Matteo Calonaci, Federica Furnari     |
|   Giancarlo Gentile, Filippo Manelli, Stefano Rebughini                 |
|                                                                         |
|   Last modified: 16 Apr 2013                                            |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of catalyticFOAM.                                   |
|                                                                         |
|   catalyticFOAM is free software: you can redistribute it and/or modify |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   catalyticFOAM is distributed in the hope that it will be useful,      |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with catalyticFOAM. If not, see <http://www.gnu.org/licenses/>. |
|                                                                         |
\*-----------------------------------------------------------------------*/

{
	//- Initial conditions
	scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	const scalarField& rhoCells = rho.internalField();
	const scalarField& vCells = mesh.V();

	if (homogeneousReactions == true)
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
		const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies()+1;
		
		// Min and max values
		Eigen::VectorXd yMin(NEQ); for(unsigned int i=0;i<NEQ;i++) yMin(i) = 0.;  yMin(NC) = 200.;
		Eigen::VectorXd yMax(NEQ); for(unsigned int i=0;i<NEQ;i++) yMax(i) = 1.;  yMax(NC) = 6000.;
		Eigen::VectorXd y0(NEQ);
		Eigen::VectorXd yf(NEQ);

		Info <<" * Solving homogeneous chemistry... "<<endl;
		{			
			unsigned int counter = 0;
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			forAll(TCells, celli)
			{
			    //- Allocating initial conditions
				if (isNeighbourCell[celli] == 0)
				{
					//- Solving for celli:	
					{
						for(unsigned int i=0;i<NC;i++)
							y0(i) = Y[i].internalField()[celli];
						y0(NC) = TCells[celli];
						
						//batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], pCells[celli], rhoCells[celli]);
						batchReactorHomogeneousConstantPressure.SetReactor(pCells[celli]);
						OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_OpenSMOKE odebatch(batchReactorHomogeneousConstantPressure);
						OpenSMOKE::StiffOdeSolverObject_Dense o( y0, t0, &odebatch, odeParameterBatchReactorHomogeneous.linear_algebra() );
			
						o.SetMinimumConstraints(yMin);
						o.SetMaximumConstraints(yMax);
						o.SetTollAbs(odeParameterBatchReactorHomogeneous.absolute_tolerance());
						o.SetTollRel(odeParameterBatchReactorHomogeneous.relative_tolerance());
						yf = o(tf,tf);
					}
	
	
					// Check mass fractions
					normalizeFractions(yf, celli, massFractionsTol, runTime.value());

					// Assign mass fractions
					for(int i=0;i<NC;i++)
				  		Y[i].internalField()[celli] = yf(i);

					//- Allocating final values: temperature
					if (energyEquation == true)
						TCells[celli] = yf(NC);

					if (counter%(int(0.20*totalNumberOfNonCatalyticReactors)+1) == 0)
						Info <<"   Accomplished: " << counter << "/" << totalNumberOfNonCatalyticReactors << endl;

					counter++;
				}
			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(totalNumberOfNonCatalyticReactors)*1000. << " ms per reactor)" << endl;
		}
	}
    
    /*
	if (totalNumberOfCatalyticReactors != 0 && heterogeneousReactions == true)
	{
		Info <<" * Solving heterogeneous chemistry... "<<endl;
		
		int counter = 0;
		
		int NC = thermodynamicsMapXML->NumberOfSpecies();
		int NS = surface_chemistry->NumberSiteSpecies();
		double* yAdditional = new double[1+NC];

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		forAll(TCells, celli)
		{
		    //- Allocating initial conditions
			if (isNeighbourCell[celli] > 0)
			{		
				//int NEQ = NC+1+NS*isNeighbourCell[celli];
				int NEQ = NC+1+NS;	
				double* yF = new double[NEQ];
				double* y0 = new double[NEQ];
				
				for(unsigned int i=0;i<NC;i++)
					y0[i] = Y[i].internalField()[celli];
				y0[NC] = TCells[celli];

				//- Cycle on number of edges
				double totalCatalyticAreaCell=0.0;
				for (int edge=1; edge<=catalyticReactors[counter].numberOfFaces; edge++)
					totalCatalyticAreaCell = totalCatalyticAreaCell + catalyticReactors[counter].faceArea[edge];
					
					
													

				int edge = 1;
				int baseCounter = NC+1+NS*(edge-1);
					for(unsigned int i=0; i<NS;i++)
			    			y0[baseCounter+(i+1)-1] = siteFraction[i].boundaryField()[catalyticReactors[counter].patchIndex[edge]][catalyticReactors[counter].faceIndex[edge]];
	

				odeSystemObject_Heterogeneous->ode.SetPressure(pCells[celli]);
				odeSystemObject_Heterogeneous->ode.SetVolume(catalyticReactors[counter].volume);
				odeSystemObject_Heterogeneous->ode.SetFaceAreas(totalCatalyticAreaCell);

				for(int i=0;i<NEQ;i++)
					odeSystemObject_Heterogeneous->ode.Y0[i+1] = y0[i];

				ode_Heterogeneous.SetMinimumConstraints(odeSystemObject_Heterogeneous->ode.YMin);
				ode_Heterogeneous.SetMaximumConstraints(odeSystemObject_Heterogeneous->ode.YMax);
				
				//ode_Heterogeneous.SetTollAbs(1e-12);
				//ode_Heterogeneous.SetTollRel(1e-7);
				

				ode_Heterogeneous(odeSystemObject_Heterogeneous->ode.Y0, t0, odeSystemObject_Heterogeneous);
				odeSystemObject_Heterogeneous->ode.YF = ode_Heterogeneous(tf,tf);
				
				FromBzzToC(odeSystemObject_Heterogeneous->ode.YF, &yF[0]);
				
				yAdditional[0] = 0. ; //  cal/m2*s
				for(int i=1;i<=NC;i++)
					yAdditional[i] = 0.;    //  moli/m2*s
				
				//- Allocating final values: mass fractions
				for(int i=1;i<=NC;i++)
					omega[i] = yF[i-1];
	
				// Check mass fractions
				normalizeFractions(omega, celli, massFractionsTol, runTime.value());

				// Assign mass fractions
				for(int i=0;i<NC;i++)
				    Y[i].internalField()[celli] = omega[i+1];

				//- Allocating final values: temperature
				if (energyEquation == true)
					TCells[celli] = yF[NC];
				
				//- Allocating final values: surface
				for (int edge=1; edge<=catalyticReactors[counter].numberOfFaces; edge++)
				{
					int iPatch = catalyticReactors[counter].patchIndex[edge];
					int iFace  = catalyticReactors[counter].faceIndex[edge];

					//- Allocating final values: temperature
					T.boundaryField()[iPatch][iFace] = yF[NC];
					grad_T.boundaryField()[iPatch][iFace] = 0.0; //yAdditional[0];

					//- Allocating final values: mass fractions on boundary
					for(int i=0;i<NC;i++)
					{
						Y[i].boundaryField()[iPatch][iFace] = yF[i];
						Ygrad[i].boundaryField()[iPatch][iFace] = 0.0;//yAdditional[i+1];
					}

					for(int i=0; i<NS;i++)
						siteFraction[i].boundaryField()[iPatch][iFace] = yF[NC+2+i-1];  
						//siteFraction[i].boundaryField()[iPatch][iFace] = yF[i+NC+(edge-1)*NS+2-1];
				}				

				if (counter%(int(0.20*totalNumberOfCatalyticReactors)+1) == 0)
					Info <<"   Accomplished: " << counter << "/" << totalNumberOfCatalyticReactors << endl;

				counter++;

				delete y0;
				delete yF;
			}
			
		}

		delete yAdditional;
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		Info << "   Heterogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(totalNumberOfCatalyticReactors)*1000. << " ms per reactor)" << endl;
	}
	else 
	{
		Info <<" * No heterogeneous chemistry since no catalytic cells are present... "<<endl;
	}
	*/

	//- Spark
	if((spark == true)&&(runTime.value() >= ignitionTime)&&(runTime.value() <= (ignitionTime+duration))) 
	{
		for(int i = 0; i<nIgnCells; i++)
		{
		    T[ignCellsList[i]] = ignitionTemperature;
		}
	}
	
	Info << endl;
    	Info <<" * Chemistry solved with success!"<<endl;
    	Info<< " * T gas min/max   = " << min(T).value() << ", " << max(T).value() << endl;
}







				//- Solving for celli:	
			/*	{
				odeSystemObject_Heterogeneous->ode.SetVolume(catalyticReactors[counter].volume);
				odeSystemObject_Heterogeneous->ode.SetFaceAreas(catalyticReactors[counter].faceArea);
				odeSystemObject_Heterogeneous->ode.SetPressure(pCells[celli]);

				#if BZZMATH_SUPPORT == 1

					for(unsigned int i=0;i<NC;i++)
			    			odeSystemObject_Heterogeneous->ode.Y0[i+1] = Y[i].internalField()[celli];
					odeSystemObject_Heterogeneous->ode.Y0[NC+1] = TCells[celli];

					//- Cycle on number of edges
					for (int edge=1; edge<=catalyticReactors[counter].numberOfFaces; edge++)
					{
						int baseCounter = NC+1 + NS*(edge-1);

						for(int i=0; i<NS;i++)
						    odeSystemObject_Heterogeneous->ode.Y0[baseCounter+(i+1)] = 
						siteFraction[i].boundaryField()[catalyticReactors[counter].patchIndex[edge]][catalyticReactors[counter].faceIndex[edge]];
					}

					ode_Heterogeneous_1.SetMinimumConstraints(odeSystemObject_Heterogeneous->ode.YMin);
					ode_Heterogeneous_1.SetMaximumConstraints(odeSystemObject_Heterogeneous->ode.YMax);

					//FromOSToBzz(X0Catalytic[isNeighbourCell[celli]], bzz_X0Catalytic[isNeighbourCell[celli]]);
					ode_Heterogeneous_1(odeSystemObject_Heterogeneous->ode.Y0, t0, odeSystemObject_Heterogeneous);
					odeSystemObject_Heterogeneous->ode.YF = ode_Heterogeneous_1(tf,tf);
					FromBzzToC(odeSystemObject_Heterogeneous->ode.YF, &yF[0]);
					
				#else
					for(int i=0;i<NC;i++)
			     			odeSystemObject_Heterogeneous->ode.Y0[i] = Y[i].internalField()[celli];
		       		 	odeSystemObject_Heterogeneous->ode.Y0[NC] = TCells[celli];

					//- Cycle on number of edges
					for (int edge=1; edge<=catalyticReactors[counter].numberOfFaces; edge++)
					{
						int baseCounter = NC+1+NS*(edge-1);

						for(int i=0; i<NS;i++)
				    			odeSystemObject_Heterogeneous->ode.Y0[baseCounter+(i+1)-1] = siteFraction[i].boundaryField()[catalyticReactors[counter].patchIndex[edge]][catalyticReactors[counter].faceIndex[edge]];
					}	

					ode_Heterogeneous_1.SetMaximumNumberOfSteps(10000);
					ode_Heterogeneous_1.SetAnalyticalJacobian(false);
					//ode_Heterogeneous_1.SetAbsoluteTolerance(aTol);
					//ode_Heterogeneous_1.SetRelativeTolerance(rTol);
					ode_Heterogeneous_1.SetInitialValues(t0, odeSystemObject_Heterogeneous->ode.Y0);
					ode_Heterogeneous_1.Solve(tf);
					ode_Heterogeneous_1.Solution(yF);
				#endif

				yAdditional[0] = -odeSystemObject_Heterogeneous_1->ode.QReaction();
				for(int i=1;i<=NC;i++)
					yAdditional[i] = odeSystemObject_Heterogeneous_1->ode.Rgas()[i];
				}
*/
				//	for (int i=0;i<N;i++)
				//		std::cout << y0[i] << " " << yF[i] << std::endl;
				
