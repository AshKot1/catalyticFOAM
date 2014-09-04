// Check ode solver type
{
	if ( odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE ||
             odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
	{
		#if OPENSMOKE_USE_DVODE != 1
		{
			Info << "The solver was compiled without the DVODE support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if ( odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE ||
                  odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
	{
		#if OPENSMOKE_USE_SUNDIALS != 1
		{
			Info << "The solver was compiled without the CVODE support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if ( odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE ||
	          odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA || 
		  odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE ||
	          odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
	{
		#if OPENSMOKE_USE_ODEPACK != 1
		{
			Info << "The solver was compiled without the ODEPACK support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if ( odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5 ||
                  odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5 )
	{ 
		#if OPENSMOKE_USE_RADAU != 1
		{
			Info << "The solver was compiled without the RADAU support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}	
	else if ( odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK ||
                  odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK )
	{
		#if OPENSMOKE_USE_DASPK != 1
		{
			Info << "The solver was compiled without the DASPK support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if ( odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF ||
                  odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
	{
		#if OPENSMOKE_USE_MEBDF != 1
		{
			Info << "The solver was compiled without the MEBDF support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}			
}
