// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file LPTemplate.hpp
  \brief Header file for the Mesquite::LPTemplate class
  \author Michael Brewer
  \date   2002-05-23
 */


#ifndef LPTemplate_hpp
#define LPTemplate_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>

namespace Mesquite
{
     /*! \class LPTemplate
       \brief Calculates the L_p objective function.  That is, sums
       the p_th powers of the (absolute value of the) quality metric values
       and takes that value to the (1/p)_th power.
     */
   class PatchData;
   class MsqMeshEntity;
   class LPTemplate :public ObjectiveFunction
   {
	public:
	  LPTemplate(QualityMetric *, int, MsqError &);
	  ~LPTemplate();
	  virtual double concrete_evaluate(PatchData &patch, MsqError &err);
	protected:
     virtual void  compute_analytical_gradient(PatchData &patch,
                                               Vector3D *const &grad,
                                               MsqError &err, int array_size);
     
	private:
     double compute_function(double metric_values[], int total_num, MsqError &err);
       //! \param pVal
	  int pVal;
   };
   
   inline double LPTemplate::compute_function(double metric_values[],
                                              int total_num, MsqError &err)
   {
     int ind=0;
     int jnd=0;
     double temp_value=1;
     double total_value=0;
     for(ind=0;ind<total_num;++ind){
       temp_value=1;
       for(jnd=0;jnd<pVal;++jnd){
         temp_value*=metric_values[ind];
       }
       total_value+=temp_value;
     }
     return pow(total_value, 1/((double) pVal));
   }
   
}//namespace



#endif // LPTemplate_hpp

