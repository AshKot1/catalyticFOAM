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

void normalizeFractions(Eigen::VectorXd& omega_plus_temperature, label patchi, label facei, double massFractionsTol)
{
	double sumFractions = 0.;
	for(int i=0; i < omega_plus_temperature.size()-1; i++)
            sumFractions += omega_plus_temperature(i);

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of site-fractions = " << sumFractions << " in patch " << patchi<< " on face " << facei <<endl;

        if (time == 0)
        {
           Info << "\nFATAL ERROR: sum of Yi is not 1 (" << sumFractions;
           Info << ") on patch " << patchi;
           Info << " on face " << facei;
           Info << "\nProcess terminated" << endl;
           abort();
        }
    }

	 for(int i=0; i < omega_plus_temperature.size()-1; i++)
            omega_plus_temperature(i) /= sumFractions;
}

void normalizeFractions(Eigen::VectorXd& omega_plus_temperature, label celli, double massFractionsTol)
{
	double sumFractions = 0.;
    for(int i=0; i < omega_plus_temperature.size()-1; i++)
            sumFractions += omega_plus_temperature(i);

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of mass-fractions = " << sumFractions << " in cell " << celli <<endl;

        for(int i=0; i < omega_plus_temperature.size()-1; i++)
            cout << i << ")\t" << omega_plus_temperature(i) << endl;
        
        if (time == 0)
        {
            Info << "\nFATAL ERROR: sum of Yi is not 1" <<endl;
            Info << "Check internal field on " << celli <<endl;
            abort();
        }    
    }
	
	 for(int i=0; i < omega_plus_temperature.size()-1; i++)
            omega_plus_temperature(i) /= sumFractions;
}

tmp<volScalarField> Rho(volScalarField& p, volScalarField& psi)
{
	return p*psi;
}
