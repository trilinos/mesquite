/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file CompositeQualityMetric.hpp
\brief Header file for the Mesquite::CompositeQualityMetric class

  \author Thomas Leurent
  \date   2002-09-01
 */


#ifndef CompositeQualityMetric_hpp
#define CompositeQualityMetric_hpp

#include "Mesquite.hpp"
#include "QualityMetric.hpp"

namespace Mesquite
{
   class MsqMeshEntity;
     /*! \class CompositeQualityMetric
       \brief Parent class for the Composite Quality Metrics.

       Contains private data members qMetric1, qMetric2,
       and scaleAlpha, and protected member functions to access
       these data members.
     */
   class CompositeQualityMetric : public QualityMetric
   {
  public:
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~CompositeQualityMetric()
        {};
     
  protected:
     //!The constructor is protected so that this class can not be
     //!instantiated.
     CompositeQualityMetric():QualityMetric(), qMetric1(NULL),
			      qMetric2(NULL), scaleAlpha(0.0)
     {};
     
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
//******************BEGIN INLINE FUNCTIONS ******************************


#endif // CompositeQualityMetric_hpp
