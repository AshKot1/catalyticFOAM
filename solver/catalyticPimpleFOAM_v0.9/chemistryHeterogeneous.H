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
