// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file CompositeOFMultiply.hpp

Header file for the Mesquite:: CompositeOFMultiply class

  \author Michael Brewer
  \date   2002-05-23
 */


#ifndef CompositeOFMultiply_hpp
#define CompositeOFMultiply_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>

namespace Mesquite
{
   /*!\class CompositeOFMultiply
       \brief Multiplies two ObjectiveFunction values together.
     */
   class MsqMeshEntity;
   class PatchData;
   class CompositeOFMultiply : public ObjectiveFunction
   {
	public:
	   CompositeOFMultiply(ObjectiveFunction*, ObjectiveFunction*);
	   ~CompositeOFMultiply();
	  virtual double concrete_evaluate(PatchData &patch, MsqError &err);
     virtual std::list<QualityMetric*> get_quality_metric_list();
	protected:
     
	private:
     ObjectiveFunction* objFunc1;
     ObjectiveFunction* objFunc2;
   };
}//namespace
#endif //  CompositeOFMultiply_hpp
