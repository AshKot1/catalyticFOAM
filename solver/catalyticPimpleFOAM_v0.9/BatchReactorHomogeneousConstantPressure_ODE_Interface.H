
#ifndef BatchReactorHomogeneousConstantPressure_ODE_Interface_H
#define	BatchReactorHomogeneousConstantPressure_ODE_Interface_H

#include "math/OpenSMOKEVector.h"
#include "ode/OpenSMOKE_OdeSystemObject.h"

#if OPENSMOKE_USE_DVODE == 1
#include "ode/OpenSMOKE_DVODE_Interface.h"
#include "ode/OpenSMOKE_DVODE.h"
#endif

#if OPENSMOKE_USE_ODEPACK == 1
#include "ode/OpenSMOKE_DLSODE_Interface.h"
#include "ode/OpenSMOKE_DLSODE.h"
#include "ode/OpenSMOKE_DLSODA_Interface.h"
#include "ode/OpenSMOKE_DLSODA.h"
#endif

#if OPENSMOKE_USE_DASPK == 1
#include "ode/OpenSMOKE_DASPK_Interface.h"
#include "ode/OpenSMOKE_DASPK.h"
#endif

#if OPENSMOKE_USE_RADAU == 1
#include "ode/OpenSMOKE_RADAU_Interface.h"
#include "ode/OpenSMOKE_RADAU.h"
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
#include "ode/OpenSMOKE_CVODE_Sundials_Interface.h"
#include "ode/OpenSMOKE_CVODE_Sundials.h"
#endif

#if OPENSMOKE_USE_MEBDF == 1
#include "ode/OpenSMOKE_MEBDF_Interface.h"
#include "ode/OpenSMOKE_MEBDF.h"
#endif

#include "ode/OpenSMOKE_OdeSystemObject.h"
#include "ODESystemVirtualClassWithOpenSMOKEVectors.H"
#include "BatchReactorHomogeneousConstantPressure.H"

namespace OpenSMOKE
{
	class BatchReactorHomogeneousConstantPressure_ODE_OpenSMOKE : public ODESystemVirtualClassWithOpenSMOKEVectors
	{
	public:

		BatchReactorHomogeneousConstantPressure_ODE_OpenSMOKE() {};

		void SetReactor(BatchReactorHomogeneousConstantPressure* reactor) { reactor_ = reactor; }

		virtual void GetEquations(const OpenSMOKE::OpenSMOKEVectorDouble& y, const double t, OpenSMOKE::OpenSMOKEVectorDouble& dy)
		{
			reactor_->Equations(t, y, dy);
		}
		virtual void PrintResults(const OpenSMOKE::OpenSMOKEVectorDouble& y, const double t)
		{
			reactor_->Print(t, y);
		}

	private:

		BatchReactorHomogeneousConstantPressure* reactor_;
	};
	
	#if OPENSMOKE_USE_DVODE == 1
	class BatchReactorHomogeneousConstantPressure_ODE_DVODE : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DVODE(BatchReactorHomogeneousConstantPressure_ODE_DVODE)

		BatchReactorHomogeneousConstantPressure* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactorHomogeneousConstantPressure(BatchReactorHomogeneousConstantPressure* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DVODE(BatchReactorHomogeneousConstantPressure_ODE_DVODE)
	#endif

	#if OPENSMOKE_USE_ODEPACK == 1
	class BatchReactorHomogeneousConstantPressure_ODE_DLSODE : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DLSODE(BatchReactorHomogeneousConstantPressure_ODE_DLSODE)

		BatchReactorHomogeneousConstantPressure* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactorHomogeneousConstantPressure(BatchReactorHomogeneousConstantPressure* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DLSODE(BatchReactorHomogeneousConstantPressure_ODE_DLSODE)

	class BatchReactorHomogeneousConstantPressure_ODE_DLSODA : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_DLSODA(BatchReactorHomogeneousConstantPressure_ODE_DLSODA)

		BatchReactorHomogeneousConstantPressure* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactorHomogeneousConstantPressure(BatchReactorHomogeneousConstantPressure* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DLSODA(BatchReactorHomogeneousConstantPressure_ODE_DLSODA)
	#endif

	#if OPENSMOKE_USE_DASPK == 1
	class BatchReactorHomogeneousConstantPressure_ODE_DASPK : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		BatchReactorHomogeneousConstantPressure* batch_;

		DEFINE_ODESOLVERINTERFACE_DASPK(BatchReactorHomogeneousConstantPressure_ODE_DASPK)

		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactorHomogeneousConstantPressure(BatchReactorHomogeneousConstantPressure* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_DASPK(BatchReactorHomogeneousConstantPressure_ODE_DASPK)
	#endif

	#if OPENSMOKE_USE_RADAU == 1
	class BatchReactorHomogeneousConstantPressure_ODE_RADAU5 : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_RADAU(BatchReactorHomogeneousConstantPressure_ODE_RADAU5)

		BatchReactorHomogeneousConstantPressure* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactorHomogeneousConstantPressure(BatchReactorHomogeneousConstantPressure* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_RADAU(BatchReactorHomogeneousConstantPressure_ODE_RADAU5)
	#endif

	#if OPENSMOKE_USE_SUNDIALS == 1
	class BatchReactorHomogeneousConstantPressure_ODE_CVODE : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_CVODE_Sundials(BatchReactorHomogeneousConstantPressure_ODE_CVODE)

		BatchReactorHomogeneousConstantPressure* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactorHomogeneousConstantPressure(BatchReactorHomogeneousConstantPressure* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_CVODE_Sundials(BatchReactorHomogeneousConstantPressure_ODE_CVODE)
	#endif

	#if OPENSMOKE_USE_MEBDF == 1
	class BatchReactorHomogeneousConstantPressure_ODE_MEBDF : public OpenSMOKE::OpenSMOKE_OdeSystemObject
	{
		DEFINE_ODESOLVERINTERFACE_MEBDF(BatchReactorHomogeneousConstantPressure_ODE_MEBDF)

		BatchReactorHomogeneousConstantPressure* batch_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	
	public:

		void SetBatchReactorHomogeneousConstantPressure(BatchReactorHomogeneousConstantPressure* batch)
		{
			batch_ = batch;
			ChangeDimensions(batch_->NumberOfEquations(), &y_, true);
			ChangeDimensions(batch_->NumberOfEquations(), &dy_, false);
		}
	 
		int GetSystemFunctions(const double t, double* y,  double* dy)
		{
			y_.CopyFrom(y);
			int flag = batch_->Equations(t, y_, dy_);
			dy_.CopyTo(dy);
			return(flag);
		}
	
		int GetAnalyticalJacobian(const double t,  double* y,  double* J)
		{
			return(0);
		}
	 
		int GetWriteFunction(const double t, double *y)
		{
			y_.CopyFrom(y);
			int flag = batch_->Print(t, y_);
			return 0;
 		}
	}; 
	COMPLETE_ODESOLVERINTERFACE_MEBDF(BatchReactorHomogeneousConstantPressure_ODE_MEBDF)
	#endif
}

#endif	// BatchReactorHomogeneousConstantPressure_ODE_Interface_H

