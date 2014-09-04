// ODE Solver (constant pressure)
typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_OpenSMOKE> denseHeterogeneousOdeConstantPressure;
typedef OdeSMOKE::MethodGear<denseHeterogeneousOdeConstantPressure> methodGearHeterogeneousConstantPressure;
OdeSMOKE::MultiValueSolver<methodGearHeterogeneousConstantPressure> odeHeterogeneousSolverConstantPressure;
odeHeterogeneousSolverConstantPressure.SetReactor(&batchReactorHeterogeneousConstantPressure);

// ODE Solver (constant volume)
typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_OpenSMOKE> denseHeterogeneousOdeConstantVolume;
typedef OdeSMOKE::MethodGear<denseHeterogeneousOdeConstantVolume> methodGearHeterogeneousConstantVolume;
OdeSMOKE::MultiValueSolver<methodGearHeterogeneousConstantVolume> odeHeterogeneousSolverConstantVolume;	
odeHeterogeneousSolverConstantVolume.SetReactor(&batchReactorHeterogeneousConstantVolume);

#if OPENSMOKE_USE_DVODE == 1
OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DVODE *dvodeHeterogeneousSolverConstantPressure;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeHeterogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DVODE::GetInstance();
	dvodeHeterogeneousSolverConstantPressure->SetBatchReactorHeterogeneousConstantPressure(&batchReactorHeterogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DVODE *dvodeHeterogeneousSolverConstantVolume;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeHeterogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DVODE::GetInstance();
	dvodeHeterogeneousSolverConstantVolume->SetBatchReactorHeterogeneousConstantVolume(&batchReactorHeterogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_ODEPACK == 1
OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DLSODE *dlsodeHeterogeneousSolverConstantPressure;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeHeterogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DLSODE::GetInstance();
	dlsodeHeterogeneousSolverConstantPressure->SetBatchReactorHeterogeneousConstantPressure(&batchReactorHeterogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DLSODE *dlsodeHeterogeneousSolverConstantVolume;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeHeterogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DLSODE::GetInstance();
	dlsodeHeterogeneousSolverConstantVolume->SetBatchReactorHeterogeneousConstantVolume(&batchReactorHeterogeneousConstantVolume);
}
OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DLSODA *dlsodaHeterogeneousSolverConstantPressure;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaHeterogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DLSODA::GetInstance();
	dlsodaHeterogeneousSolverConstantPressure->SetBatchReactorHeterogeneousConstantPressure(&batchReactorHeterogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DLSODA *dlsodaHeterogeneousSolverConstantVolume;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaHeterogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DLSODA::GetInstance();
	dlsodaHeterogeneousSolverConstantVolume->SetBatchReactorHeterogeneousConstantVolume(&batchReactorHeterogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_CVODE *cvodeHeterogeneousSolverConstantPressure;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeHeterogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_CVODE::GetInstance();
	cvodeHeterogeneousSolverConstantPressure->SetBatchReactorHeterogeneousConstantPressure(&batchReactorHeterogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_CVODE *cvodeHeterogeneousSolverConstantVolume;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeHeterogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_CVODE::GetInstance();
	cvodeHeterogeneousSolverConstantVolume->SetBatchReactorHeterogeneousConstantVolume(&batchReactorHeterogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_DASPK == 1
OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DASPK *daspkHeterogeneousSolverConstantPressure;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkHeterogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkHeterogeneousSolverConstantPressure->SetBatchReactorHeterogeneousConstantPressure(&batchReactorHeterogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DASPK *daspkHeterogeneousSolverConstantVolume;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkHeterogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkHeterogeneousSolverConstantVolume->SetBatchReactorHeterogeneousConstantVolume(&batchReactorHeterogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_RADAU == 1
OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_RADAU5 *radauHeterogeneousSolverConstantPressure;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radauHeterogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_RADAU5::GetInstance();
	radauHeterogeneousSolverConstantPressure->SetBatchReactorHeterogeneousConstantPressure(&batchReactorHeterogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_RADAU5 *radauHeterogeneousSolverConstantVolume;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radauHeterogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_RADAU5::GetInstance();
	radauHeterogeneousSolverConstantVolume->SetBatchReactorHeterogeneousConstantVolume(&batchReactorHeterogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_MEBDF == 1
OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_MEBDF *mebdfHeterogeneousSolverConstantPressure;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfHeterogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_MEBDF::GetInstance();
	mebdfHeterogeneousSolverConstantPressure->SetBatchReactorHeterogeneousConstantPressure(&batchReactorHeterogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_MEBDF *mebdfHeterogeneousSolverConstantVolume;
if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfHeterogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_MEBDF::GetInstance();
	mebdfHeterogeneousSolverConstantVolume->SetBatchReactorHeterogeneousConstantVolume(&batchReactorHeterogeneousConstantVolume);
}
#endif
