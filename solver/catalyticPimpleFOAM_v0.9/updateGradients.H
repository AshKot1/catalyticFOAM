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
/*	OpenSMOKE::OpenSMOKEVectorDouble MWi = gas_mixture->M();

	forAll(T.boundaryField(), patchi)
	{
		for(int j=0; j<indexCatalyst.size(); j++)
            	{
                	if (patchi == indexCatalyst[j])
			{
				fvPatchScalarField& pp = p.boundaryField()[patchi];
				fvPatchScalarField& pT = T.boundaryField()[patchi];
				fvPatchScalarField& ppsi  = psi.boundaryField()[patchi];
				fvPatchScalarField& plambda = lambda.boundaryField()[patchi];

				forAll(pT, facei)
				{
                    			for(label i=0; i<Ygrad.size(); i++)
						Ygrad[i].boundaryField()[patchi][facei] *=  
                          	  		MWi(i+1)*AlfaCatalyst/Dmix[i].boundaryField()[patchi][facei]/(pp[facei]*ppsi[facei]); 
			
                    			if(energyEquation == true)
                     	  			grad_T.boundaryField()[patchi][facei] *= AlfaCatalyst/plambda[facei];
				}
			}
                }
	}
*/
}
