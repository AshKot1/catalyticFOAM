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

tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    
    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {

        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
	    	volScalarField& Dmixi = Dmix[i];

	    	if (diffusivityCorrection != "none")
			{
				dimensionedScalar Mi( "Mi", dimensionSet(1,0,0,0,-1,0,0), thermodynamicsMapXML->MW()[i+1] ); 

				sumDiffusionCorrections *= 0.;
				for (label j=0; j<Y.size(); j++)
				{
					volScalarField& Yj = Y[j];	
					dimensionedScalar Mj( "Mj", dimensionSet(1,0,0,0,-1,0,0), thermodynamicsMapXML->MW()[j+1] );
					sumDiffusionCorrections += fvc::laplacian(rho*PMmix/Mj*Dmixi*Yi, Yj);
				
					if ( (i==j) && (diffusivityCorrection == "semi-implicit") )
						sumDiffusionCorrections -= fvc::laplacian(rho*PMmix/Mj*Dmixi*Yi, Yj);
				}

				if (diffusivityCorrection == "semi-implicit")
				{
					fvScalarMatrix YiEqn
					(
						fvm::ddt(rho, Yi)
					  + mvConvection->fvmDiv(phi, Yi)
					  - fvm::laplacian(rho*Dmixi, Yi)
					  + fvm::laplacian(rho*Dmixi*PMmix*Yi/Mi, Yi)
					  + sumDiffusionCorrections
					 ==
						fvOptions(rho, Yi)
					);

					YiEqn.relax();
					fvOptions.constrain(YiEqn);
					YiEqn.solve(mesh.solver("Yi"));
					fvOptions.correct(Yi);						
				}
				else if (diffusivityCorrection == "explicit")
				{					
					fvScalarMatrix YiEqn
					(
						fvm::ddt(rho, Yi)
					  + mvConvection->fvmDiv(phi, Yi)
					  - fvm::laplacian(rho*Dmixi, Yi)
					  + sumDiffusionCorrections
					 ==
						fvOptions(rho, Yi)
					);

					YiEqn.relax();
					fvOptions.constrain(YiEqn);
					YiEqn.solve(mesh.solver("Yi"));
					fvOptions.correct(Yi);			
				}
			}
			else
			{
				fvScalarMatrix YiEqn
				(
					fvm::ddt(rho, Yi)
				  + mvConvection->fvmDiv(phi, Yi)
				  - fvm::laplacian(rho*Dmixi, Yi)
				 ==
					fvOptions(rho, Yi)
				);

				YiEqn.relax();
				fvOptions.constrain(YiEqn);
				YiEqn.solve(mesh.solver("Yi"));
				fvOptions.correct(Yi);			
			}

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
		    inertIndex = i;
	    }
    }

    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	Info << "   Transported solved in " << tEnd - tStart << " s " << endl;

}
