// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file MaxTemplate.hpp

Header file for the Mesquite::MaxTemplate class

  \author Lori Freitag
  \date   2002-07-18
 */


#ifndef MaxTemplate_hpp
#define MaxTemplate_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>

namespace Mesquite
{
   class PatchData;
   /*! \class MaxTemplate
    \brief Computes the maximum quality metric value.

    This function is the same as the LInfTemplate except that no
    absolute values are used.
  */
   class MaxTemplate :public ObjectiveFunction
   {
   public:
     MaxTemplate(QualityMetric *);
     virtual ~MaxTemplate();
     virtual bool concrete_evaluate(PatchData &patch, double &fval,
                                    MsqError &err);
   protected:
   private:
	  
   };
}//namespace
#endif // MaxTemplate_hpp
