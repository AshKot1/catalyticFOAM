/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

{
    const volVectorField& centres = mesh.C();

    label ignitionCell;
    ignitionCell = mesh.findCell(position);
    
    if (ignitionCell == -1)
    {
        Info << "FOAM ERROR:\n theorical ignition position does not correspond to any cell.\n Please check position\nIn file constant/solverOptions.H\n"<<endl;
        abort();
    }
    else Info << "Ignition cell center: " << centres[ignitionCell] << ". Original position: " << position <<endl;

    scalar radius = ignitionDiameter/2.0;
    ignCellsList.setSize(1);
    ignCellsList[0] = ignitionCell;

    forAll(centres, celli)
    {
        scalar dist = mag(centres[celli] - position);

        if (dist < radius && celli != ignitionCell)
        {
            ignCellsList.setSize(nIgnCells+1);

            ignCellsList[nIgnCells] = celli;

            nIgnCells++;
        }
    }

    if (ignCellsList.size() > 0) 
	Info << "Added " << nIgnCells << " cells to ignition point.\n" <<endl;
    else
	Info << "The solver wasn't able to find more than one cell in the ignition point.\n" <<endl;
}
