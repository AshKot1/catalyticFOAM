{
	label inertIndex = -1;
    	volScalarField Yt(0.0*Y[0]);

	forAll(Y, i)
    	{
	        if (Y[i].name() != inertSpecie)
	        {
		    volScalarField& Yi = Y[i];
		    Yi.max(0.0);
	            Yt += Yi;
	        }
	        else
	        {
	            inertIndex = i;
	        }
	}

    	Y[inertIndex] = scalar(1.0) - Yt;
    	Y[inertIndex].max(0.0);
}
