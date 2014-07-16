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

    if(energyEquation == true)
    {
		sumCpSpecies *= 0.;
		if (iMassDiffusionContributionInEnergyEquation == true)
		{

			for (label i=0; i<Y.size(); i++)
			{
				volScalarField& Yi 		= Y[i];			
				volScalarField& Dmixi 	= Dmix[i];
				volScalarField& Cpi 	= CpSpecies[i];
		
				sumCpSpecies += rho*Cpi*Dmixi*(fvc::grad(Yi)&fvc::grad(T));
			}
		}

		fvScalarMatrix TEqn
		(
			cp*fvm::ddt(rho,T) + 
			cp*fvm::div(phi, T)
		 ==
		    fvm::laplacian(lambda,T) + 
		    sumCpSpecies +
			cp*fvOptions(rho, T)
		);

		TEqn.relax();
		fvOptions.constrain(TEqn);
		TEqn.solve(mesh.solver("T"));
		fvOptions.correct(T);			

    }

    Info<< " * T gas min/max   = " << min(T).value() << ", " << max(T).value() << endl;

    if( (spark == true) && (runTime.value() >= ignitionTime) && (runTime.value() <= (ignitionTime+duration)) ) 
    {
        for(int i = 0; i<nIgnCells; i++)
        {
            T[ignCellsList[i]] = ignitionTemperature;
        }
    }
}
