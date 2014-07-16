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

	Info << endl;
	Info<< " * Reading chemical species mole fractions or mass fractions..." << endl;

	PtrList<volScalarField> Y(thermodynamicsMapXML->NumberOfSpecies());
	PtrList<volScalarField> Dmix(thermodynamicsMapXML->NumberOfSpecies());
	PtrList<volScalarField> Ygrad(thermodynamicsMapXML->NumberOfSpecies());
	PtrList<volScalarField> CpSpecies(thermodynamicsMapXML->NumberOfSpecies());

	for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
	{

		Info << "   Reading species: " << thermodynamicsMapXML->NamesOfSpecies()[i] << " > ";

		IOobject header
		(
			thermodynamicsMapXML->NamesOfSpecies()[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		if (header.headerOk())
		{
			Info << thermodynamicsMapXML->NamesOfSpecies()[i] << endl;
			Y.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);

			Ygrad.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"grad_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, -1, 0, 0, 0)
				)
			);
		}
		else
		{
			Info << "Ydefault" << endl;

			volScalarField Ydefault
			(
				IOobject
				(
					"Ydefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);

			Y.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					Ydefault
				)
			);

			Ygrad.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"grad_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, 0, -1, 0, 0)
				)
			);
		}

		// Sets of diffusivities and specific heats of species
		{
			Dmix.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"Dmix_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, 2, -1, 0, 0)
				)
			);
		
		    CpSpecies.set
		    (
				i,
				new volScalarField
				(
					IOobject
					(
					"CpSpecies_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, 2, -2, -1, 0)
				)
			);
		}

	}

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
	//	Viscosity and mass diffusion coefficients initialization
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
	// The fields of viscosity and mass diffusion coefficients are initialized according to
	// the values suggested by the user

	mu = dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 1.8e-5);
	for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		Dmix[i] = dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 1.8e-5);

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

