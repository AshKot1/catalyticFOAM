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

#ifndef BatchReactorHeterogeneousConstantPressure_H
#define	BatchReactorHeterogeneousConstantPressure_H

class BatchReactorHeterogeneousConstantPressure
{
public:

	BatchReactorHeterogeneousConstantPressure(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
												OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
												OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>& thermodynamicsSurfaceMap, 
												OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>& kineticsSurfaceMap);

	void SetReactor( const double P0, const double A, const double alfaCatalyst);
	
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

	unsigned int NumberOfEquations() { return NE_; }

	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

private:

	double T_;
	double P0_;
	double A_;
	double alfaCatalyst_;
	double MW_;
	double mass_;

	unsigned int NC_;
	unsigned int SURF_NP_;
	unsigned int SURF_NC_;
	unsigned int NE_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&			thermodynamicsMap_;				//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN<double>&					kineticsMap_;					//!< kinetic map
	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&	thermodynamicsSurfaceMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&			kineticsSurfaceMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;
	OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;
	OpenSMOKE::OpenSMOKEVectorDouble Z_;
	OpenSMOKE::OpenSMOKEVectorDouble Gamma_;
	OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;
	OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;
	OpenSMOKE::OpenSMOKEVectorDouble dummy;
	
	double QRGas_;
	double QRSurface_;
	
	bool checkMassFractions_;
	bool energyEquation_;
};

BatchReactorHeterogeneousConstantPressure::BatchReactorHeterogeneousConstantPressure(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
																						OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
																						OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>& thermodynamicsSurfaceMap, 
																						OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>& kineticsSurfaceMap):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap),
	thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap), 
	kineticsSurfaceMap_(kineticsSurfaceMap)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
		SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
		NE_ = NC_ + SURF_NC_ + SURF_NP_ + 1 + 1;
		
		QRGas_     = 0.;
		QRSurface_ = 0.;
				
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &RfromGas_, true);
		ChangeDimensions(NC_, &RfromSurface_, true);
		
		ChangeDimensions(SURF_NC_, &Z_, true);
		ChangeDimensions(SURF_NP_, &Gamma_, true);
		ChangeDimensions(SURF_NC_, &Rsurface_, true);
		ChangeDimensions(SURF_NP_, &RsurfacePhases_, true);
		
		checkMassFractions_ = false;
		energyEquation_ = true;
	}

void BatchReactorHeterogeneousConstantPressure::SetReactor( const double P0, const double A, const double alfaCatalyst)
{
	P0_    = P0;
	A_     = A;
	alfaCatalyst_ = alfaCatalyst;
}

int BatchReactorHeterogeneousConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	// Recover unknowns
	{
		unsigned int k=1;
		if (checkMassFractions_ == true)
		{
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = std::max(y[k++], 0.);
			mass_ = std::max(y[k++], 0.);
			for(unsigned int i=1;i<=SURF_NP_;++i)
				Gamma_[i] = std::max(y[k++], 0.);
			for(unsigned int i=1;i<=SURF_NC_;++i)
				Z_[i] = std::max(y[k++], 0.);
		}
		else
		{
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = y[k++];
			mass_ = y[k++];
			for(unsigned int i=1;i<=SURF_NP_;++i)
				Gamma_[i] = y[k++];
			for(unsigned int i=1;i<=SURF_NC_;++i)
				Z_[i] = y[k++];
		}
		
		T_ = y[k++];
	}
	
	// Calculates the volume and the concentrations of species
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
	const double cTot_ = P0_/(PhysicalConstants::R_J_kmol * T_);
	Product(cTot_, x_, &c_);
	const double rho_ = cTot_*MW_;
	const double V_ = mass_/rho_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P0_);
		
	// Calculates homogeneous kinetics
	{
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P0_);

		kineticsMap_.ReactionRates(c_);
		kineticsMap_.FormationRates(&RfromGas_);
	}		
		
	// Calculates heterogeneous kinetics
	{
		thermodynamicsSurfaceMap_.SetPressure(P0_);
		thermodynamicsSurfaceMap_.SetTemperature(T_);
		kineticsSurfaceMap_.SetPressure(P0_);
		kineticsSurfaceMap_.SetTemperature(T_);
		kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
		kineticsSurfaceMap_.ArrheniusKineticConstants();
		kineticsSurfaceMap_.ReactionRates(c_, Z_, dummy, Gamma_);
		kineticsSurfaceMap_.FormationRates(&RfromSurface_, &Rsurface_, &dummy, &RsurfacePhases_);
	}
	
	double dm_over_dt = A_*alfaCatalyst_*Dot(RfromSurface_, thermodynamicsMap_.MW());

	// Recovering residuals
	{
		unsigned int k=1;
		
		// Gas phase species
		for (unsigned int i=1;i<=NC_;++i)	
			dy[k++] = thermodynamicsMap_.MW()[i]*RfromGas_[i]/rho_ +
					 (-omega_[i]*dm_over_dt + A_*alfaCatalyst_*RfromSurface_[i]*thermodynamicsMap_.MW()[i])/mass_;
		
		// Total mass
		dy[k++] = dm_over_dt;
		
		// Phases
		for (unsigned int i=1;i<=SURF_NP_;++i)	
			dy[k++] = RsurfacePhases_[i];
		
		// Surface site species
		for (unsigned int i=1;i<=SURF_NC_;++i)	
		{
			const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i-1]+1;
			dy[k++] = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[i-1]*Rsurface_[i] - 
									Z_[i]*dy[NC_+1+index_phase])/Gamma_[index_phase];
		}
		
		// Energy equation
		if (energyEquation_ == true)
		{
			double CpMixMolar; 
			thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
			const double CpMixMass_ = CpMixMolar / MW_;
			
			QRGas_ = kineticsMap_.HeatRelease(RfromGas_);
			QRSurface_ = kineticsSurfaceMap_.HeatRelease(RfromSurface_, Rsurface_, dummy);
			
			dy[k++] = (V_*QRGas_+ A_*alfaCatalyst_*QRSurface_) / (mass_*CpMixMass_);
		}
		else
		{
			dy[k++] = 0.;
		}
	}
		
	return 0;
}

int BatchReactorHeterogeneousConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
}

#endif // BatchReactorHeterogeneousConstantPressure_H

