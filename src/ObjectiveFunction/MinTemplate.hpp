// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file MinTemplate.hpp

Header file for the Mesquite::MinTemplate class

  \author Lori Freitag
  \date   2002-07-18
 */


#ifndef MinTemplate_hpp
#define MinTemplate_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>

namespace Mesquite
{
   class PatchData;
   class MinTemplate :public ObjectiveFunction
   {
	public:
	  MinTemplate(QualityMetric *);
	  ~MinTemplate();
	  virtual double concrete_evaluate(PatchData &patch, MsqError &err);
	protected:
	private:
	  
   };
}//namespace
#endif // MinTemplate_hpp
