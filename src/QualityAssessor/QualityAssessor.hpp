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

/*! \file QualityAssessor.hpp

Header file for the Mesquite::QualityAssessor class

  \author Thomas Leurent
  \date   2002-05-01
 */


#ifndef MSQ_QUALITYASSESSOR_HPP
#define MSQ_QUALITYASSESSOR_HPP


#include "Mesquite.hpp"
#include "PatchDataUser.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <string.h>
#else
#  include <list>
#  include <string>
   using std::string;
   using std::list;
#endif


namespace Mesquite 
{

   class QualityMetric;
   class MsqError;
   class MeshSet;

  /*! \class QualityAssessor

      \brief The QualityAssessor class contains all the information needed
      to perform a single sweep over the mesh and make all relevant
      quality assessments.

      The relevant quality assessments are set by the user or
      automatically (default) by Mesquite when an InstructionQueue
      object is used.  If the mesh has been changed (improved),
      it is often useful to reuse the same QualityAssessor object
      to reassess the mesh quality.
  */
   class QualityAssessor : public PatchDataUser
  {
  public:
    /*! \enum QAFunction
      type of function used in conjunction with QualityMetric to compute mesh quality */ 
    enum QAFunction {
       AVERAGE=1,
       HISTOGRAM=2,
       MAXIMUM=4, 
       MINIMUM=8,
       RMS=16,
       STDDEV=32,
       ALL_MEASURES=63
    };
    string get_QAFunction_name(enum QualityAssessor::QAFunction);
    
    //! Constructor requires a QualityMetric and an evaluation function
    QualityAssessor(QualityMetric*, enum QAFunction  func,
                    string name = "DefaultQualAssessName");

      //!Destructor
    ~QualityAssessor();
    
      //! Provides a name to the QualityAssessor (use it for default name in constructor).
    void set_name(string name) { qualityAssessorName = name; };
      //! Retrieves the QualityAssessor name. A default name should be set in the constructor.
    virtual string get_name() { return qualityAssessorName; }

    virtual AlgorithmType get_algorithm_type() { return QUALITY_ASSESSOR; }
      //! Adds a quality metric and a wrapper function (min, max, ...).
    void add_quality_assessment(QualityMetric* qm, enum QAFunction  func, MsqError &err);
    
      //! Set the min and max values to be used for the histogram.
    void set_histogram_range(QualityMetric* qm, double min_val, double max_val, MsqError &err);
    
      //! Does one sweep over the mesh and assess the quality with the metrics previously added.
    virtual double loop_over_mesh(MeshSet &ms, MsqError &err);

      //! Do not print results of assessment.
    void disable_printing_results()
       {
         printingTurnedOff=1;
       }

      /*!Sets the QualityMetric and QAFunction combination that will
        be returned when loop_over_mesh is called.
      */
    void set_stopping_assessment(QualityMetric* qm, enum QAFunction func,
                                 MsqError &err);
    
  private:
    string qualityAssessorName;  
    struct Assessor
    {
      QualityMetric* metric;
      long unsigned int funcFlagBits;
      double minHist;
      double maxHist;

      // structure constructor
      Assessor() :
        metric(0),
        funcFlagBits(0),
        minHist(0),
        maxHist(0)
      {}
    };

      //QAVars is used in loop_over_mesh.  It gives
      //us a place to store data as we are looping over
      //the mesh
    
    struct QAVars
    {
      double avgVar;//holds sum of metric vals until print
      int histVar[MSQ_HIST_SIZE+2];//holds number per col of hist
      double histMax;
      double histMin;
      double histDelta;//holds step size for histogram
      double maxVar;//hold max metric val
      double minVar;//holds min metric val
      double rmsVar;//hold sum of squares until print
      double stdVar;//holds std dev squared until print
      int numInvalid;//counts the number of invalid metric values
    };
    
    list<Assessor*> assessList;
      //flag to turn off printing
    int printingTurnedOff;
      //pointer to qm used form return value
    QualityMetric* stoppingMetric;
    QAFunction stoppingFunction;
  };

  
} //namespace


#endif // QualityAssessor_hpp
