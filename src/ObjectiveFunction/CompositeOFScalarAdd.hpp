// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file CompositeOFScalarAdd.hpp

Header file for the Mesquite:: CompositeOFScalarAdd class

  \author Michael Brewer
  \date   2002-06-24
 */


#ifndef CompositeOFScalarAdd_hpp
#define CompositeOFScalarAdd_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>


namespace Mesquite
{
   /*! \class CompositeOFScalarAdd.
       \brief Adds a scalar to a given ObjectiveFunction.
     */
   class MsqMeshEntity;
   class PatchData;
   class CompositeOFScalarAdd : public ObjectiveFunction
   {
	public:
      CompositeOFScalarAdd(double, ObjectiveFunction*);
	   ~CompositeOFScalarAdd();
	  virtual double concrete_evaluate(PatchData &patch, MsqError &err);
     virtual std::list<QualityMetric*> get_quality_metric_list();
	protected:
     
	private:
     double alpha;
     ObjectiveFunction* objFunc;
   };
}//namespace
#endif //  CompositeOFScalarAdd_hpp
