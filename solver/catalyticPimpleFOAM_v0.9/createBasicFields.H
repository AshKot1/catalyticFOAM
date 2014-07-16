Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Reading pressure [Pa]
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Reading temperature [K]
Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Temperature gradient [K/s]
Info<< "Creating field dTdt\n" << endl;
volScalarField grad_T
(
    IOobject
    (
        "grad_T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet (0, -1, 0, 1, 0)
);

//- Allocating memory for mean molecular weight [kg/mol]
Info<< "Creating field MWmix\n" << endl;
volScalarField PMmix
(
    IOobject
    (
        "meanMolecularWeight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 0, 0, 0, -1)
);

//- Allocating memory for total concentration [mol/m3]
//- Not necessary
Info<< "Creating field cTot\n" << endl;
volScalarField cTot
(
    IOobject
    (
        "cTot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, -3, 0, 0, 1)
);


//- Allocating memory for p/rho [s2/m2]
Info<< "Creating field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0, -2, 2, 0, 0)
);

//- Allocating memory for dynamic viscosity [kg/m/s]
Info<< "Creating field mu\n" << endl;
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -1, 0, 0)
);

//- Allocating memory for thermal conductivity [W/m/K]
Info<< "Creating field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(1, 1, -3, -1, 0)
);

//- Allocating memory for mean specific heat [J/kg/K]
Info<< "Creating field cp\n" << endl;
volScalarField cp
(
    IOobject
    (
        "specificHeats",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
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
        "enthalpy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, 0, 0)
);


//- Allocating memory for radiative heat transfer term [W/m3]
Info<< "Creating field Qrad\n" << endl;
volScalarField Qrad
(
    IOobject
    (
        "Qrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -3, 0, 0)
);

//- Allocating memory for surface heat flux [W/m2]
Info<< "Creating field Qsurf\n" << endl;
volScalarField Qsurf
(
    IOobject
    (
        "Qsurf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 2, -3, 0, 0)
);

Info<< "Creating sumCorrections\n" << endl;
volScalarField sumDiffusionCorrections
(
	IOobject
	(
	"sumCorrections",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(1, -3, -1, 0, 0)
);

Info<< "Creating sumCpSpecies\n" << endl;
volScalarField sumCpSpecies
(
    IOobject
    (
	"sumCp",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -3, 0, 0)
);
