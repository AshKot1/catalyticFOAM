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
	
	Info<< "Reading Solver options\n" << endl;
	IOdictionary solverOptions
	(
		IOobject
		(
			"solverOptions",
			U.time().constant(),
			U.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
	
	// Read the kinetic scheme in XML format

	label inertIndex;
	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>*			thermodynamicsMapXML; 
	OpenSMOKE::KineticsMap_CHEMKIN<double>* 			kineticsMapXML;
	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>*		surfaceThermodynamicsMapXML; 
	OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>* 		surfaceKineticsMapXML;
	OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>* 		transportMapXML;
	
	const dictionary& kineticsDictionary = solverOptions.subDict("Kinetics");
	{
		//- Kinetics folder
		Foam::string kinetics_folder= kineticsDictionary.lookup("folder");
		
		// Reading homogeneous data
		Info<< " * loading homogeneous kinetics...\n" << endl;
		{	
			word kinetics("opensmoke");
			boost::filesystem::path path_kinetics = kinetics_folder;

			rapidxml::xml_document<> doc;
			std::vector<char> xml_string;
			OpenSMOKE::OpenInputFileXML(doc, xml_string,path_kinetics / "kinetics.xml");

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
			transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
			kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc); 
								
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
		}
		
		// Reading heterogeneous data
		Info<< " * loading heterogeneous kinetics...\n" << endl;
		{	
			word kinetics("opensmoke");
			boost::filesystem::path path_kinetics = kinetics_folder;

			rapidxml::xml_document<> doc;
			std::vector<char> xml_string;
			OpenSMOKE::OpenInputFileXML(doc, xml_string,path_kinetics / "kinetics.surface.xml");

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			surfaceThermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>(doc); 
			surfaceKineticsMapXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>(*surfaceThermodynamicsMapXML, doc); 
								
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
		}	
		
		//- Inert species
		word inertSpecies(kineticsDictionary.lookup("inertSpecies"));
		inertIndex = thermodynamicsMapXML->IndexOfSpecies(inertSpecies)-1;
	}
	
	
	//- Mass fractions tolerance
	scalar massFractionsTol = 1e-10;
	scalar siteFractionsTol = 1e-6;

	//- Internal models
	enum { DIFFUSIVITY_CORRECTION_NONE, DIFFUSIVITY_CORRECTION_SEMI_IMPLICIT, DIFFUSIVITY_CORRECTION_EXPLICIT} diffusivityCorrection;
	enum { STRANG_REACTION_TRANSPORT_REACTION, STRANG_REACTION_TRANSPORT, STRANG_TRANSPORT_REACTION, STRANG_REACTION_TRANSPORT_HYBRID, STRANG_TRANSPORT_REACTION_MOMENTUM} strangAlgorithm;
	enum { DIFFUSIVITY_MODEL_MULTICOMPONENT, DIFFUSIVITY_MODEL_LEWIS_NUMBERS} diffusivityModel;

	// Physical model
	Switch energyEquation;
	Switch homogeneousReactions;
	Switch heterogeneousReactions;
	Switch iMassDiffusionInEnergyEquation;
	Switch constPressureBatchReactor;
	Switch includeDpDt;
	scalar AlfaCatalyst;
	Eigen::VectorXd LewisNumbers;
	const dictionary& physicalModelDictionary = solverOptions.subDict("PhysicalModel");
	{
		energyEquation = Switch(physicalModelDictionary.lookup(word("energyEquation")));
		homogeneousReactions = Switch(physicalModelDictionary.lookup(word("homogeneousReactions")));
		heterogeneousReactions = Switch(physicalModelDictionary.lookup(word("heterogeneousReactions")));
		includeDpDt = Switch(physicalModelDictionary.lookup(word("includeDpDt")));
		iMassDiffusionInEnergyEquation = Switch(physicalModelDictionary.lookup(word("massDiffusionInEnergyEquation")));

		// Diffusivity 
		{
			word diffusivity(physicalModelDictionary.lookup("diffusivityCorrection"));
			if (diffusivity == "none")					diffusivityCorrection = DIFFUSIVITY_CORRECTION_NONE;
			else if (diffusivity == "semi-implicit")	diffusivityCorrection = DIFFUSIVITY_CORRECTION_SEMI_IMPLICIT;
			else if (diffusivity == "explicit")			diffusivityCorrection = DIFFUSIVITY_CORRECTION_EXPLICIT;
			else
			{
				Info << "Wrong diffusivityCorrection option: none || semi-implicit || explicit" << endl;
				abort();
			}
		}
		
		// Algorithm
		{
			word strang(physicalModelDictionary.lookup("strangAlgorithm"));
			     if (strang == "ReactionTransportReaction")	strangAlgorithm = STRANG_REACTION_TRANSPORT_REACTION;
			else if (strang == "ReactionTransport")			strangAlgorithm = STRANG_REACTION_TRANSPORT;
			else if (strang == "TransportReaction")			strangAlgorithm = STRANG_TRANSPORT_REACTION;
			else if (strang == "ReactionTransportHybrid")	strangAlgorithm = STRANG_REACTION_TRANSPORT_HYBRID;
			else if (strang == "TransportReactionMomentum")	strangAlgorithm = STRANG_TRANSPORT_REACTION_MOMENTUM;
			
			else
			{
				Info << "Wrong strangAlgorithm option: ReactionTransportReaction || TransportReaction || ReactionTransport || ReactionTransportHybrid || TransportReactionMomentum" << endl;
				abort();
			}
		}
		
		// Diffusivity of species
		{
			word diffusivity(physicalModelDictionary.lookup("diffusivityModel"));
			if (diffusivity == "multi-component")		diffusivityModel = DIFFUSIVITY_MODEL_MULTICOMPONENT;
			else if (diffusivity == "lewis-numbers")	diffusivityModel = DIFFUSIVITY_MODEL_LEWIS_NUMBERS;
			else
			{
				Info << "Wrong diffusivityModel option: multi-component || lewis-numbers" << endl;
				abort();
			}
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				LewisNumbers.resize(thermodynamicsMapXML->NumberOfSpecies());
				
				const dictionary& lewisNumberDictionary = solverOptions.subDict("LewisNumbers");
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					LewisNumbers(i) = readScalar(lewisNumberDictionary.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
				
				Info << endl;
				Info << "Diffusivities of species will be calculated using the following Lewis numbers: " << endl;
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Info << " " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t" << LewisNumbers(i) << endl;
				Info << endl;
			}
		}
		
		//- Constant pressure reactors
		constPressureBatchReactor = Switch(physicalModelDictionary.lookup(word("constPressureBatchReactor")));
		
		// Catalyst properties
		AlfaCatalyst = readScalar(physicalModelDictionary.lookup("alfaCatalyst"));
	}

	//- Detect spark
	const dictionary& sparkDictionary = solverOptions.subDict("Spark");
	sparkModel spark;
	spark.Read(sparkDictionary);

	// Batch reactor homogeneous: ode parameters
	const dictionary& odeHomogeneousDictionary = solverOptions.subDict("OdeHomogeneous");
	OpenSMOKE::ODE_Parameters odeParameterBatchReactorHomogeneous;
	{
		//- Mass fractions tolerance
		scalar relTolerance = readScalar(odeHomogeneousDictionary.lookup("relTolerance"));
		scalar absTolerance = readScalar(odeHomogeneousDictionary.lookup("absTolerance"));
		odeParameterBatchReactorHomogeneous.SetRelativeTolerance(relTolerance);
		odeParameterBatchReactorHomogeneous.SetAbsoluteTolerance(absTolerance);
	
		//- Full pivoting (only for OpenSMOKE solver)
		Switch fullPivoting(odeHomogeneousDictionary.lookup("fullPivoting"));
		odeParameterBatchReactorHomogeneous.SetFullPivoting(fullPivoting);
	
		//- Maximum order of integration (only for OpenSMOKE solver)
		label maximumOrder = readLabel(odeHomogeneousDictionary.lookup("maximumOrder"));
		odeParameterBatchReactorHomogeneous.SetMaximumOrder(maximumOrder);
	
		//- 
	
		// Type
		word homogeneousODESolverString(odeHomogeneousDictionary.lookup("odeSolver"));
		if (	homogeneousODESolverString != "OpenSMOKE" 	&& homogeneousODESolverString != "DVODE"  && 
				homogeneousODESolverString != "DLSODE" 		&& homogeneousODESolverString != "DLSODA" && 
				homogeneousODESolverString != "CVODE" 		&& homogeneousODESolverString != "DASPK"  &&
				homogeneousODESolverString != "MEBDF" 		&& homogeneousODESolverString != "RADAU5"  
		   )
		{
			Info << "Wrong homogeneous ODE Solver: OpenSMOKE || DVODE || DLSODE || DLSODA || CVODE || DASPK || MEBDF || RADAU5" << endl;
			abort();
		}

		if (homogeneousODESolverString == "OpenSMOKE")	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE);
		if (homogeneousODESolverString == "DVODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE);
		if (homogeneousODESolverString == "DLSODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE);
		if (homogeneousODESolverString == "DLSODA") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA);
		if (homogeneousODESolverString == "CVODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE);
		if (homogeneousODESolverString == "DASPK") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK);
		if (homogeneousODESolverString == "MEBDF") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF);
		if (homogeneousODESolverString == "RADAU5") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5);	
	
	}

	// Batch reactor heterogeneous: ode parameters
	const dictionary& odeHeterogeneousDictionary = solverOptions.subDict("OdeHeterogeneous");
	OpenSMOKE::ODE_Parameters odeParameterBatchReactorHeterogeneous;
	{
		//- Mass fractions tolerance
		scalar relTolerance = readScalar(odeHeterogeneousDictionary.lookup("relTolerance"));
		scalar absTolerance = readScalar(odeHeterogeneousDictionary.lookup("absTolerance"));
		odeParameterBatchReactorHeterogeneous.SetRelativeTolerance(relTolerance);
		odeParameterBatchReactorHeterogeneous.SetAbsoluteTolerance(absTolerance);
	
		//- Full pivoting (only for OpenSMOKE solver)
		Switch fullPivoting(odeHeterogeneousDictionary.lookup("fullPivoting"));
		odeParameterBatchReactorHeterogeneous.SetFullPivoting(fullPivoting);
	
		//- Maximum order of integration (only for OpenSMOKE solver)
		label maximumOrder = readLabel(odeHeterogeneousDictionary.lookup("maximumOrder"));
		odeParameterBatchReactorHeterogeneous.SetMaximumOrder(maximumOrder);
	
		//- 
	
		// Type
		word heterogeneousODESolverString(odeHeterogeneousDictionary.lookup("odeSolver"));
		if (	heterogeneousODESolverString != "OpenSMOKE" 	&& heterogeneousODESolverString != "DVODE"  && 
				heterogeneousODESolverString != "DLSODE" 		&& heterogeneousODESolverString != "DLSODA" && 
				heterogeneousODESolverString != "CVODE" 		&& heterogeneousODESolverString != "DASPK"  &&
				heterogeneousODESolverString != "MEBDF" 		&& heterogeneousODESolverString != "RADAU5"  
		   )
		{
			Info << "Wrong heterogeneous ODE Solver: OpenSMOKE || DVODE || DLSODE || DLSODA || CVODE || DASPK || MEBDF || RADAU5" << endl;
			abort();
		}

		if (heterogeneousODESolverString == "OpenSMOKE")	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE);
		if (heterogeneousODESolverString == "DVODE") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE);
		if (heterogeneousODESolverString == "DLSODE") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE);
		if (heterogeneousODESolverString == "DLSODA") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA);
		if (heterogeneousODESolverString == "CVODE") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE);
		if (heterogeneousODESolverString == "DASPK") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK);
		if (heterogeneousODESolverString == "MEBDF") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF);
		if (heterogeneousODESolverString == "RADAU5") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5);	
	
	}

	#include "checkOdeSolverType.H"

	Info<< "Solver options correctly read\n" << endl;
