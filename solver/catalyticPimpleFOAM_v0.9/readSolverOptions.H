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

	//- Inert species
	word inertSpecie(solverOptions.lookup("inertSpecie"));

	//- Kinetics folder
	Foam::string kinetics_folder = solverOptions.lookup("kineticsFolder");

	//- Energy equations: on/off
	Switch energyEquation(solverOptions.lookup("energyEquation"));

	//- Homogeneous reactions: on/off
	Switch homogeneousReactions(solverOptions.lookup("homogeneousReactions"));

	//- Heterogeneous reactions: on/off
	Switch heterogeneousReactions(solverOptions.lookup("heterogeneousReactions"));

	//- Detect spark
	Switch spark(solverOptions.lookup("spark"));
	Foam::vector position;
	scalar ignitionTime = 0.;
	scalar duration = 0.;
	scalar ignitionTemperature = 0.;
	scalar ignitionDiameter = 0.;

	//-  Read spark properties if detected
	if(spark == true)
	{
		position = (solverOptions.lookup("position"));
		ignitionTime = readScalar(solverOptions.lookup("time"));
		duration = readScalar(solverOptions.lookup("duration"));
		ignitionTemperature = readScalar(solverOptions.lookup("temperature"));
		ignitionDiameter = readScalar(solverOptions.lookup("diameter"));   
	}

 /*********************************************************************************************************************/ 

	// Read the kinetic scheme in XML format
	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>*				thermodynamicsMapXML; 
	OpenSMOKE::KineticsMap_CHEMKIN<double>* 					kineticsMapXML;
	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>*		surfaceThermodynamicsMapXML; 
	OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>* 			surfaceKineticsMapXML;
	OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>* 			transportMapXML;
	
	// Reading homogeneous data
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

	// Catalyst properties
	scalar AlfaCatalyst = readScalar(solverOptions.lookup("alfaCatalyst"));

	//- Mass fractions tolerance
	scalar massFractionsTol = readScalar(solverOptions.lookup("massFractionsTolerance"));

	// MassDiffusionContributionInEnergyEquation
	Switch iMassDiffusionContributionInEnergyEquation(solverOptions.lookup("massDiffusionContributionInEnergyEquation"));

	// Diffusivity 
	word diffusivityCorrection(solverOptions.lookup("diffusivityCorrection"));
	if (diffusivityCorrection != "none" && diffusivityCorrection != "semi-implicit" && diffusivityCorrection != "explicit")
	{
		Info << "Wrong diffusivityCorrection option" << endl;
		abort();
	}

	// Algorithm
	word strangAlgorithm(solverOptions.lookup("strangAlgorithm"));
	if (	strangAlgorithm != "ReactionTransportReaction" 	&& strangAlgorithm != "ReactionTransport" && 
			strangAlgorithm != "TransportReaction" 		&& strangAlgorithm != "ReactionTransportHybrid" &&
			strangAlgorithm != "TransportReactionMomentum" 
		)
	{
		Info << "Wrong strangAlgorithm option: ReactionTransportReaction || TransportReaction || ReactionTransport || ReactionTransportHybrid || TransportReactionMomentum" << endl;
		abort();
	}
