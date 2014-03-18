	{
		Info << "Properties..." << endl;
		const scalarField& pCells = p.internalField();
		const scalarField& TCells = T.internalField();

		scalarField& CpCells = cp.internalField();
		scalarField& mwCells = MWmix.internalField();
		scalarField& hCells = h.internalField();
		
		scalarField& muCells = mu.internalField();
		scalarField& lambdaCells = lambda.internalField();
		//scalarField& QCells = this->Q_.internalField();

		OpenSMOKE::OpenSMOKEVectorDouble y(Y.size());
		OpenSMOKE::OpenSMOKEVectorDouble x(Y.size());
		OpenSMOKE::OpenSMOKEVectorDouble c(Y.size());
		OpenSMOKE::OpenSMOKEVectorDouble gammamix(Y.size());
		OpenSMOKE::OpenSMOKEVectorDouble tetamix(Y.size());
		//OpenSMOKE::OpenSMOKEVectorDouble R(Y.size());
		OpenSMOKE::OpenSMOKEVectorDouble cp_species(Y.size());

		forAll(TCells, celli)
		{
			if (TCells[celli] < 0. || pCells[celli]<0.)
				Info << "Wrong T and P: " << celli << " " << TCells[celli] << " " << pCells[celli] << endl;
					
			// Sets the maps
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			thermodynamicsMapXML->SetPressure(pCells[celli]);

			// Extract the mass fractions
			for(unsigned int i=1;i<=Y.size();i++)
				y[i] = Y[i-1].internalField()[celli];

			// Calculate the molecular weight and the mole fractions
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mwCells[celli], y);	
			
			// Calculate the density
			double const cTot = pCells[celli] / (PhysicalConstants::R_J_kmol * TCells[celli]);
			double const rho = cTot * mwCells[celli];

			// Calculate constant pressure and constant volume specific heats
			double CpMixMolar; 
			thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(CpMixMolar,x);
			CpCells[celli] = CpMixMolar / mwCells[celli];
			
			if (massDiffusionInEnergyEquation_ == true)
			{
				thermodynamicsMapXML->cpMolar_Species_From_MoleFractions(cp_species);
				ElementByElementDivision(cp_species, thermodynamicsMapXML->MW(), &cp_species);
				for(unsigned int i=1;i<=Y.size();i++)
					CpSpecies_[i-1].internalField()[celli] = cp_species[i];
			}			
			

			// Calculate the enthalpy
			thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(hCells[celli], x);
			hCells[celli] /= mwCells[celli];
			
	
			// Accurate transport properties
			{
				// Effective Thermal conductivity [kg/m/s]
				transportMapXML->ThermalConductivity(lambdaCells[celli], TCells[celli], x);
				//alphaCells[celli] = lambdaCells[celli]/CpCells[celli];
				
				// Dynamic viscosity [kg/m/s]
				transportMapXML->DynamicViscosity(muCells[celli], TCells[celli], x);
				
				// Mass diffusion coefficients [m2/s]
				transportMapXML->MassDiffusionCoefficients(gammamix, TCells[celli], pCells[celli], x);
				
				// Thermal diffusion coefficients [-]
				if (soretEffect_ == true)
					transportMapXML->ThermalDiffusionRatios(tetamix, TCells[celli], x);
					
				for(unsigned int i=1;i<=Y.size();i++)
					Dmix_[i-1].internalField()[celli] = gammamix[i];
			}
			
			// Calculate formation rates
			/*
			if (calculate_reactions == true && homogeneous_reactions_ == true)
			{
				Product(cTot,x,&c);
				
				kineticsMapXML->SetPressure(pCells[celli]);
				kineticsMapXML->SetTemperature(TCells[celli]);
				
				kineticsMapXML->ReactionEnthalpiesAndEntropies();
				kineticsMapXML->ArrheniusKineticConstants();
				kineticsMapXML->ReactionRates(c);
				kineticsMapXML->FormationRates(&R);
				QCells[celli] = kineticsMapXML->HeatRelease(R);
				ElementByElementProduct(R,thermodynamicsMapXML->MW(),&R);
				
				for(unsigned int i=1;i<=Y.size();i++)
					RHomo_[i-1].internalField()[celli] = R[i];
			}
			*/
		}

		forAll(T.boundaryField(), patchi)
		{
			fvPatchScalarField& pp = p.boundaryField()[patchi];
			fvPatchScalarField& pT = T.boundaryField()[patchi];
			
			fvPatchScalarField& pmw = MWmix.boundaryField()[patchi];
			fvPatchScalarField& pcp = cp.boundaryField()[patchi];
			fvPatchScalarField& ph = h.boundaryField()[patchi];
			fvPatchScalarField& pmu = mu.boundaryField()[patchi];
			//fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];
			fvPatchScalarField& plambda = lambda.boundaryField()[patchi];
			//fvPatchScalarField& pq = this->Q_.boundaryField()[patchi];

			{
				forAll(pT, facei)
				{
					// Sets the maps
					thermodynamicsMapXML->SetTemperature(pT[facei]);
					thermodynamicsMapXML->SetPressure(pp[facei]);

					// Extract the mass fractions
					for(unsigned int i=1;i<=Y.size();i++)
						y[i] = Y[i-1].boundaryField()[patchi][facei];

					// Calculate the molecular weight and the mole fractions
					thermodynamicsMapXML->MoleFractions_From_MassFractions(x, pmw[facei], y);


					// Calculate the density
					double const cTot = pp[facei]  / (PhysicalConstants::R_J_kmol * pT[facei]);
					double const rho = cTot * pmw[facei];
					//ppsi[facei] = rho / pp[facei];

					// Calculate constant pressure and constant volume specific heats
					double CpMixMolar; 
					thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(CpMixMolar,x);
					pcp[facei]= CpMixMolar / pmw[facei];
					
					if (massDiffusionInEnergyEquation_ == true)
					{
						thermodynamicsMapXML->cpMolar_Species_From_MoleFractions(cp_species);
						ElementByElementDivision(cp_species, thermodynamicsMapXML->MW(), &cp_species);
						for(unsigned int i=1;i<=Y.size();i++)
							CpSpecies_[i-1].boundaryField()[patchi][facei] = cp_species[i];
					}

					// Calculate the enthalpy
					thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(ph[facei], x);
					ph[facei] /= pmw[facei];
					
					// Calculate accurate transport properties
					{
						// Effective Thermal conductivity [kg/m/s]
						transportMapXML->ThermalConductivity(plambda[facei], pT[facei], x);
						//palpha[facei] = plambda[facei]/pcp[facei];
						
						// Dynamic viscosity [kg/m/s]
						transportMapXML->DynamicViscosity(pmu[facei], pT[facei], x);
						
						// Mass diffusion coefficients [m2/s]
						transportMapXML->MassDiffusionCoefficients(gammamix, pT[facei], pp[facei], x);
						
						// Thermal diffusion coefficients [-]
						if (soretEffect_ == true)
							transportMapXML->ThermalDiffusionRatios(tetamix, pT[facei], x);
							
						for(unsigned int i=1;i<=Y.size();i++)
							Dmix_[i-1].boundaryField()[patchi][facei] = gammamix[i];
					}	
					
					// Calculate formation rates
					/*
					if (calculate_reactions == true && homogeneous_reactions_ == true)
					{
						Product(cTot,x,&c);
						
						kineticsMapXML->SetPressure(pp[facei]);
						kineticsMapXML->SetTemperature(pT[facei]);
						
						kineticsMapXML->ReactionEnthalpiesAndEntropies();
						kineticsMapXML->ArrheniusKineticConstants();
						kineticsMapXML->ReactionRates(c);
						kineticsMapXML->FormationRates(&R);
						pq[facei] = kineticsMapXML->HeatRelease(R);
						ElementByElementProduct(R,thermodynamicsMapXML->MW(),&R);
						
						for(unsigned int i=1;i<=Y.size();i++)
							RHomo_[i-1].boundaryField()[patchi][facei] = R[i];
					}	
					*/						
				}
			}
		}
	}
