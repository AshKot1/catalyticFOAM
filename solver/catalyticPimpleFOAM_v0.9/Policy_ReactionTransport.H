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
	#include "readTimeControls.H"
	#include "readPISOControls.H"
	#include "compressibleCourantNo.H"
	#include "setDeltaT.H"

	runTime++;
	Info<< "Time = " << runTime.timeName() << nl << endl;
	
	// Chemical step
	//#include "properties.H"
	//#include "updateGradients.H"

	{
		double t0 = runTime.value() - runTime.deltaT().value();
	   	double tf = runTime.value();

		#include "chemistry.H"
		#include "properties.H"
		#include "updateGradients.H"
	}

	#include "rhoEqn.H"

	for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
	{

		#include "UEqn.H"
		#include "YEqn.H"
		#include "TEqn.H"

		// --- PISO loop
		for (int corr=1; corr<=nCorr; corr++)
		{
			#include "pEqn.H"
		}
	}

	runTime.write();

	#include "conversion.H"


	Info	<< " ExecutionTime = " 	<< runTime.elapsedCpuTime() 	<< " s"
		<< " ClockTime = " 	<< runTime.elapsedClockTime() 	<< " s"	<< nl << endl;
}
