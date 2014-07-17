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

	//- Memory allocation: gas-phase properties
	OpenSMOKE::OpenSMOKEVectorDouble massFractions(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble moleFractions(thermodynamicsMapXML->NumberOfSpecies());

	//- Memory allocation: gas-phase chemistry
	OpenSMOKE::OpenSMOKEVectorDouble omega(thermodynamicsMapXML->NumberOfSpecies());
	
	// Batch reactor homogeneous
	BatchReactorHomogeneousConstantPressure batchReactorHomogeneousConstantPressure(*thermodynamicsMapXML, *kineticsMapXML);
	BatchReactorHomogeneousConstantVolume batchReactorHomogeneousConstantVolume(*thermodynamicsMapXML, *kineticsMapXML);
	
	
	OpenSMOKE::ODE_Parameters odeParameterBatchReactorHomogeneous;
	//odeParameterBatchReactorHomogeneousSetRelativeTolerance(1e-7);
	//odeParameterBatchReactorHomogeneousSetAbsoluteTolerance(1e-12);
	


