/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/


// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "reactors/utilities/Utilities"

// Additional classes
#include "catalyticReactorClass.H"

// ODE solvers
#include "math/multivalue-ode-solvers/MultiValueSolver"
#include "ode/ODE_Parameters.h"

// OpenFOAM
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"

// Additional include files
#include "sparkModel.H"
#include "userDefinedFunctions.H"

// Homogeneous reactors
#include "BatchReactorHomogeneousConstantPressure.H"
#include "BatchReactorHomogeneousConstantPressure_ODE_Interface.H"
#include "BatchReactorHomogeneousConstantVolume.H"
#include "BatchReactorHomogeneousConstantVolume_ODE_Interface.H"

// Heterogeneous reactors
#include "BatchReactorHeterogeneousConstantPressure.H"
#include "BatchReactorHeterogeneousConstantPressure_ODE_Interface.H"
#include "BatchReactorHeterogeneousConstantVolume.H"
#include "BatchReactorHeterogeneousConstantVolume_ODE_Interface.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    
    #include "createBasicFields.H"
    #include "readSolverOptions.H"
    #include "createAdditionalFields.H"
    #include "createCatalyticFields.H"
    #include "createFvOptions.H"
    #include "catalystTopology.H"
    #include "memoryAllocation.H"
    #include "properties.H"
    #include "finalSetupOperations.H"
    
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    
    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {
	#include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (strangAlgorithm == STRANG_REACTION_TRANSPORT)
	{
		#include "Policy_ReactionTransport.H"
	}
	else if (strangAlgorithm == STRANG_TRANSPORT_REACTION)
	{
		#include "Policy_TransportReaction.H"
	}
	else if (strangAlgorithm == STRANG_TRANSPORT_REACTION_MOMENTUM)
	{
		#include "Policy_TransportReactionMomentum.H"
	}
	else if (strangAlgorithm == STRANG_REACTION_TRANSPORT_REACTION)
	{
		#include "Policy_ReactionTransportReaction.H"
	}
	else if (strangAlgorithm == STRANG_REACTION_TRANSPORT_HYBRID)
	{
		#include "Policy_ReactionTransportHybrid.H"
	}
			
	runTime.write();

        Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
				<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
				<< nl << endl;
    }
    
    Info<< "End\n" << endl;
    
    return 0;
}

// ************************************************************************* //
