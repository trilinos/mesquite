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
       /*!Use set_dividing_by_n to control whether this objective
         function divides it's final value by the number of
         metric values used to compute the objective function
         value.  That is, if the associated metric is element
         based, the obejctive function value is divided by
         the number of elements.  If it is vertex based, the
         objective function is divided by the number of vertices.
         If this function is passed 'true', the function value
         will be scale.  If it is passed false, the function
         value will not be scaled.*/
     void set_dividing_by_n(bool d_bool){dividingByN=d_bool;}
     
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
       //! The metric value entries are raised to the pVal power
     short pVal;
       //! dividingByN is true if we are dividing the objective function
       //! by the number of metric values.
     bool dividingByN;
   };
   
   inline double LPtoPTemplate::compute_function(double metric_values[],
                                                 size_t total_num,
                                                 MsqError &err)
   {
     double scale_factor=1.0;
       //if scaling, divid by total_num
     if(dividingByN){
       if(total_num<=0)
         err.set_msg("\nIn compute_function, attempting to divide by zero.");
       else
         scale_factor/=total_num;
     }
     size_t ind=0;
     short jnd=0;
     double temp_value=1;
     double total_value=0;
     for(ind=0;ind<total_num;++ind){
       temp_value=1;
       for(jnd=0;jnd<pVal;++jnd){
         temp_value*=metric_values[ind];
       }
       total_value+=(scale_factor*temp_value);
     }
     return total_value;
   }
   
}//namespace

#endif // LPtoPTemplate_hpp

