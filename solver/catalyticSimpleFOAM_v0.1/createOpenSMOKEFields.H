// Read the kinetic scheme in XML format
OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>* thermodynamicsMapXML; 
OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>* transportMapXML;
OpenSMOKE::KineticsMap_CHEMKIN<double>* kineticsMapXML;

{	
	word kinetics("opensmoke");
	boost::filesystem::path path_kinetics = kinetics;

	rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	OpenSMOKE::OpenXMLFile(doc, xml_string,path_kinetics / "kinetics.xml");

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
	transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
	kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc); 			
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
}

// TO ADD For Detailed Kinetics
/*
OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>* thermodynamicsSurfaceMapXML;
OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>* kineticsSurfaceMapXML;
{	
	word kinetics("opensmoke");
	boost::filesystem::path path_kinetics = kinetics;

	rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	OpenSMOKE::OpenXMLFile(doc, xml_string,path_kinetics / "kinetics.surface.xml");

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	thermodynamicsSurfaceMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>(doc); 
	kineticsSurfaceMapXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>(*thermodynamicsSurfaceMapXML, doc); 					
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	std::cout << "Time to read XML file: " << tEnd-tStart << std::endl;
}

label nSurfaceSpecies = thermodynamicsSurfaceMapXML->vector_names_site_species().size();
Info << nSurfaceSpecies << endl;
*/

// TO ADD For Detailed Kinetics
/*
PtrList<volScalarField> YSurface_;
PtrList<volScalarField> RSurface_;
*/

PtrList<volScalarField> Dmix_;
PtrList<volScalarField> CpSpecies_;
PtrList<volScalarField> R_;


// Creating fields
{
	// Check the order of species
	if (thermodynamicsMapXML->NumberOfSpecies() != Y.size())
		FatalError << "Kinetic schemes do not match..." << endl;

	for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		if (Y[i].name() != thermodynamicsMapXML->NamesOfSpecies()[i])
		{
			FatalError << "Kinetic schemes do not match..." << endl;
		}
	
	R_.resize(Y.size());
	Dmix_.resize(Y.size());
	CpSpecies_.resize(Y.size());
	
	// TO ADD For Detailed Kinetics
	/*
	YSurface_.resize(nSurfaceSpecies);
	RSurface_.resize(nSurfaceSpecies);
	*/

	forAll(Y,i)
	{
	
		R_.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"R_"+ Y[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.)
			)
		);

		Dmix_.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"Dmix_"+ Y[i].name(),
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.)
			)
		);

		CpSpecies_.set
		(
			i,
			new volScalarField
			(
			    IOobject
			    (
				"add::CpSpecies_"+ Y[i].name(),
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

	// TO ADD For Detailed Kinetics
	/*
	forAll(YSurface_,i)
	{
		std::cout << " * Processing species "<< i << " : " << thermodynamicsSurfaceMapXML->vector_names_site_species()[i] << " : ";
		IOobject header
		(
			"YSurface_" + thermodynamicsSurfaceMapXML->vector_names_site_species()[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		if (header.headerOk())
		{
			std::cout << "from file..." << std::endl;
			YSurface_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"YSurface_" + thermodynamicsSurfaceMapXML->vector_names_site_species()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			std::cout << "from default values..." << std::endl;

			volScalarField YSurface_default
			(
				IOobject
				(
					"YSurface_default",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);

			YSurface_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"YSurface_" + thermodynamicsSurfaceMapXML->vector_names_site_species()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					YSurface_default
				)
			);
		}

		RSurface_.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"RSurface_" + thermodynamicsSurfaceMapXML->vector_names_site_species()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.)
			)
		);
	}
	*/
}

	//- Mean molecular weight [kg/mol]
	Info<< "Creating field MWmix\n" << endl;
	volScalarField MWmix
	(
   	 IOobject
   	 (
     	   "add::MW",
     	   runTime.timeName(),
    	    mesh,
    	    IOobject::NO_READ,
    	    IOobject::AUTO_WRITE
   	 ),
   	 mesh,
  	  dimensionSet(1, 0, 0, 0, -1)
	);

	//- Dynamic viscosity [kg/m/s]
	Info<< "Creating field mu\n" << endl;
	volScalarField mu
	(
  	  IOobject
  	  (
      	  "add::mu",
      	  runTime.timeName(),
      	  mesh,
      	  IOobject::NO_READ,
     	   IOobject::AUTO_WRITE
    	),
   	 mesh,
   	 dimensionSet(1, -1, -1, 0, 0)
	);

	//- Thermal conductivity [W/m/K]
	Info<< "Creating field lambda\n" << endl;
	volScalarField lambda
	(
   	 IOobject
   	 (
     	   "add::lambda",
     	   runTime.timeName(),
     	   mesh,
   	     IOobject::NO_READ,
    	    IOobject::AUTO_WRITE
   	 ),
  	  mesh,
 	   dimensionSet(1, 1, -3, -1, 0)
	);

	//- Specific heat [J/kg/K]
	Info<< "Creating field cp\n" << endl;
	volScalarField cp
	(
   	 IOobject
   	 (
      	  "add::Cp",
      	  runTime.timeName(),
     	   mesh,
     	   IOobject::NO_READ,
     	   IOobject::AUTO_WRITE
   	 ),
   	 mesh,
   	 dimensionSet(0, 2, -2, -1, 0)
	);

	//- Allocating memory for total enthalpy [m2/s2]
	Info<< "Creating field enthalpy\n" << endl;
	volScalarField h
	(
  	  IOobject
   	 (
    	    "add::enthalpy",
    	    runTime.timeName(),
    	    mesh,
  	      IOobject::NO_READ,
    	    IOobject::AUTO_WRITE
  	  ),
 	   mesh,
 	  dimensionSet(0, 2, -2, 0, 0)
	);

// Homogeneous reactions on/off
Switch homogeneousReactions_(thermo.lookup("homogeneousReactions"));
Switch heterogeneousReactions_(thermo.lookup("heterogeneousReactions"));
Switch massDiffusionInEnergyEquation_(thermo.lookup("massDiffusionInEnergyEquation"));
Switch soretEffect_(thermo.lookup("soretEffect"));


// Tolerances
scalar relTolerance = readScalar(thermo.lookup("relTolerance"));	// default: 1.e-7
scalar absTolerance = readScalar(thermo.lookup("absTolerance"));	// default: 1.e-12

// Spark (to ignite the mixture)
#include "spark.H"

