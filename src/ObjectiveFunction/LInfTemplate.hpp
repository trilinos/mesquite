// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file LInfTemplate.hpp

Header file for the Mesquite::LInfTemplate class

  \author Michael Brewer
  \date   2002-07-3
 */


#ifndef LInfTemplate_hpp
#define LInfTemplate_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>

namespace Mesquite
{
   class PatchData;
  /*! \class LInfTemplate
    \brief Computes the L_infinity objective function for a given patch,
    i.e., LInfTemplate::concrete_evaluate returns the maximum absolute value of
    the quality metric values  on 'patch'.
  */
   class LInfTemplate :public ObjectiveFunction
   {
	public:
	  LInfTemplate(QualityMetric *);
	  ~LInfTemplate();
	  virtual double concrete_evaluate(PatchData &patch, MsqError &err);
	protected:
	private:
	  
   };
}//namespace
#endif // LInfTemplate_hpp
