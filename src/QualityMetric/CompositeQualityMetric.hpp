// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file CompositeQualityMetric.hpp
\brief Header file for the Mesquite::CompositeQualityMetric class

  \author Thomas Leurent
  \date   2002-09-01
 */


#ifndef CompositeQualityMetric_hpp
#define CompositeQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityMetric.hpp"

namespace Mesquite
{
   class MsqMeshEntity;
     /*! \class CompositeQualityMetric
       \brief Parent class for the Composite Quality Metrics.  Contains
       private data members qMetric1, qMetric2, and scaleAlpha, and
       protected member functions to access these data members.
     */
   class CompositeQualityMetric : public QualityMetric
   {
  public:
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~CompositeQualityMetric()
        {};
     
  protected:
       //!Set qMetric1.
     inline void set_qmetric1(QualityMetric* qm){
       qMetric1=qm;
     }
       //!Get qMetric1
     inline QualityMetric*  get_qmetric1(){
       return qMetric1;
     }
     
       //!Set qMetric2.
     inline void set_qmetric2(QualityMetric* qm){
       qMetric2=qm;
     }
      //!Get qMetric2
     inline QualityMetric* get_qmetric2(){
       return qMetric2;
     }
       //!Set scaleAlpha.
     inline void set_scalealpha(double alpha){
       scaleAlpha=alpha;
     }
      //!Get scaleAlpha.
     inline double get_scalealpha(){
       return scaleAlpha;
     }
  private:
     QualityMetric* qMetric1;
     QualityMetric* qMetric2;
     double scaleAlpha; 
   };

} //namespace


#endif // CompositeQualityMetric_hpp
