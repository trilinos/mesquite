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
MSQ_USE(list);

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
     virtual ~CompositeOFScalarMultiply();
     virtual bool concrete_evaluate(PatchData &patch, double &fval,
                                    MsqError &err);
     virtual list<QualityMetric*> get_quality_metric_list();
     
   protected:
     //!Implement the scalar multiply analytic gradient
     bool compute_analytical_gradient(PatchData &patch,Vector3D *const &grad,
				      double &OF_val,MsqError &err,
				      size_t array_size);
   private:
     ObjectiveFunction* objFunc;
     double mAlpha;
   };
}//namespace
#endif //  CompositeOFScalarMultiply_hpp
