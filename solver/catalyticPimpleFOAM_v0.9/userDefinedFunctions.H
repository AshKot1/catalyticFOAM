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

tmp<volScalarField> Rho(volScalarField& p, volScalarField& psi)
{
	return p*psi;
}

void normalizeMassFractions(Eigen::VectorXd& omega_plus_temperature, const label celli, const double massFractionsTol)
{
	double sumFractions = 0.;
    for(int i=0; i < omega_plus_temperature.size()-1; i++)
            sumFractions += omega_plus_temperature(i);

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of mass fractions = " << sumFractions << " in cell " << celli <<endl;
      
        if (time == 0)
        {
				Info << "\nFATAL ERROR: sum of mass fractions is not 1" << endl;
				Info << "Check internal field on cell: " << celli << endl;
				abort();
        }    
    }
	
	 for(int i=0; i < omega_plus_temperature.size()-1; i++)
            omega_plus_temperature(i) /= sumFractions;
}

void normalizeMassAndSiteFractions(	Eigen::VectorXd& yf, const label celli, 
					const unsigned int NC, const unsigned int SURF_NP, const unsigned int SURF_NC, 
					const double massFractionsTol, const double siteFractionsTol)
{
	// Gas phase
	{
		double sum = 0.;
		for(unsigned int i=0;i<NC;i++)
			sum += yf(i);

		if (sum > (1.+massFractionsTol) || sum < (1.-massFractionsTol) )
		{
			Info << "WARNING: sum of mass fractions = " << sum << " in cell " << celli <<endl;
			
			if (time == 0)
			{
				Info << "\nFATAL ERROR: sum of mass fractions is not 1" << endl;
				Info << "Check internal field on cell: " << celli << endl;
				abort();
			}    
		}
		
		for(unsigned int i=0;i<NC;i++)
			yf(i) /= sum;
	}
	
	// Surface site species
	{
		const unsigned int count = NC+1+SURF_NP;
		
		double sum = 0.;
		for(unsigned int i=0;i<SURF_NC;i++)
			sum += yf(count+i);

		if (sum > (1.+siteFractionsTol) || sum < (1.-siteFractionsTol) )
		{
			Info << "WARNING: sum of site fractions = " << sum << " in cell " << celli <<endl;
			
			if (time == 0)
			{
				Info << "\nFATAL ERROR: sum of site fractions is not 1" << endl;
				Info << "Check internal field on cell: " << celli << endl;
				abort();
			}    
		}
		
		for(unsigned int i=0;i<SURF_NC;i++)
			yf(count+i) /= sum;
	}
}

void normalizeMassAndSiteFractions(	OpenSMOKE::OpenSMOKEVectorDouble& yf, const label celli, 
					const unsigned int NC, const unsigned int SURF_NP, const unsigned int SURF_NC, 
					const double massFractionsTol, const double siteFractionsTol)
{
	// Gas phase
	{
		double sum = 0.;
		for(unsigned int i=1;i<=NC;i++)
			sum += yf[i];

		if (sum > (1.+massFractionsTol) || sum < (1.-massFractionsTol) )
		{
			Info << "WARNING: sum of mass fractions = " << sum << " in cell " << celli <<endl;
			
			if (time == 0)
			{
				Info << "\nFATAL ERROR: sum of mass fractions is not 1" << endl;
				Info << "Check internal field on cell: " << celli << endl;
				abort();
			}    
		}
		
		for(unsigned int i=1;i<=NC;i++)
			yf[i] /= sum;
	}
	
	// Surface site species
	{
		const unsigned int count = NC+1+SURF_NP;
		
		double sum = 0.;
		for(unsigned int i=1;i<=SURF_NC;i++)
			sum += yf[count+i];

		if (sum > (1.+siteFractionsTol) || sum < (1.-siteFractionsTol) )
		{
			Info << "WARNING: sum of site fractions = " << sum << " in cell " << celli <<endl;
			
			if (time == 0)
			{
				Info << "\nFATAL ERROR: sum of site fractions is not 1" << endl;
				Info << "Check internal field on cell: " << celli << endl;
				abort();
			}    
		}
		
		for(unsigned int i=1;i<=SURF_NC;i++)
			yf[count+i] /= sum;
	}
}

void normalizeMassFractions(OpenSMOKE::OpenSMOKEVectorDouble& omega_plus_temperature, const label celli, const double massFractionsTol)
{
	double sumFractions = 0.;
    for(int i=1; i <= omega_plus_temperature.Size()-1; i++)
            sumFractions += omega_plus_temperature[i];

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of mass fractions = " << sumFractions << " in cell " << celli <<endl;
        
        if (time == 0)
        {
				Info << "\nFATAL ERROR: sum of mass fractions is not 1" << endl;
				Info << "Check internal field on cell: " << celli << endl;
				abort();
        }    
    }
	
	 for(int i=1; i <= omega_plus_temperature.Size()-1; i++)
            omega_plus_temperature[i] /= sumFractions;
}
