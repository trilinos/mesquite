// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file LPtoPTemplate.hpp
  \brief Header file for the Mesquite::LPtoPTemplate class
 \author Michael Brewer
 \author Thomas Leurent
  \date   2002-05-23
 */


#ifndef LPtoPTemplate_hpp
#define LPtoPTemplate_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>

namespace Mesquite
{
     /*! \class LPtoPTemplate
       \brief Calculates the L_p objective function raised to the pth
       power.  That is, sums the p_th powers of (the absolute value of)
       the quality metric values.
     */
   class PatchData;
   class MsqMeshEntity;
   class LPtoPTemplate :public ObjectiveFunction
   {
   public:
     LPtoPTemplate(QualityMetric *, short, MsqError &);
     virtual ~LPtoPTemplate();
     virtual bool concrete_evaluate(PatchData &patch, double &fval,
                                    MsqError &err);
   protected:
     virtual bool compute_analytical_gradient(PatchData &patch,
					      Vector3D *const &grad,
					      double &OF_val,
					      MsqError &err, 
					      size_t array_size);
     
     virtual bool  compute_analytical_hessian(PatchData &patch,
					      MsqHessian &hessian, 
					      Vector3D *const &grad,
					      double &OF_val,
					      MsqError &err);
     
   private:
     double compute_function(double metric_values[], size_t total_num,
                             MsqError &err);
       //! \param pVal
     short pVal;
   };
   
   inline double LPtoPTemplate::compute_function(double metric_values[],
                                                 size_t total_num,
                                                 MsqError &/*err*/)
   {
     size_t ind=0;
     short jnd=0;
     double temp_value=1;
     double total_value=0;
     for(ind=0;ind<total_num;++ind){
       temp_value=1;
       for(jnd=0;jnd<pVal;++jnd){
         temp_value*=metric_values[ind];
       }
       total_value+=temp_value;
     }
     return total_value;
   }
   
}//namespace

#endif // LPtoPTemplate_hpp

