// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file CompositeOFScalarMultiply.hpp

Header file for the Mesquite:: CompositeOFScalarMultiply class

  \author Michael Brewer
  \date   2002-06-24
 */


#ifndef CompositeOFScalarMultiply_hpp
#define CompositeOFScalarMultiply_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>

namespace Mesquite
{
     /*! \class CompositeOFScalarMultiply
       \brief Scales a given an ObjectiveFunction.
     */
   class MsqMeshEntity;
   class PatchData;
   class CompositeOFScalarMultiply : public ObjectiveFunction
   {
	public:
	   CompositeOFScalarMultiply(double, ObjectiveFunction*);
	   ~CompositeOFScalarMultiply();
	  virtual double concrete_evaluate(PatchData &patch, MsqError &err);
     virtual std::list<QualityMetric*> get_quality_metric_list();
	protected:
     
	private:
     ObjectiveFunction* objFunc;
     double alpha;
   };
}//namespace
#endif //  CompositeOFScalarMultiply_hpp
