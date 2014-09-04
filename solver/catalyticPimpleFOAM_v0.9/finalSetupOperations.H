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

	//- Creating density [kg/m3]
	Info<< "Creating field rho\n" << endl;
	volScalarField rho
	(
		IOobject
		(
			"rho",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		Rho(p,psi)
	);

    #include "compressibleCreatePhi.H"

	Info<< "Creating field DpDt\n" << endl;
	volScalarField DpDt
	(
		IOobject
		(
			"DpDt",
			runTime.timeName(),
			mesh
		),
		mesh,
		dimensionedScalar("DpDt", p.dimensions()/dimTime, 0)
	);

	// Kinetic energy
	Info<< "Creating field kinetic energy K\n" << endl;
	volScalarField K("K", 0.5*magSqr(U));


	multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
	forAll(Y, i)
	{
		fields.add(Y[i]);
	}

	// Spark
	spark.Setup(mesh);

