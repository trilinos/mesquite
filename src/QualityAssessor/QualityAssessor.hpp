// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file QualityAssessor.hpp

Header file for the Mesquite::QualityAssessor class

  \author Thomas Leurent
  \date   2002-05-01
 */


#ifndef QualityAssessor_hpp
#define QualityAssessor_hpp

#include <list>
#include <string>

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "MeshSet.hpp"
#include "PatchDataUser.hpp"

namespace Mesquite {

   class QualityMetric;

  /*! \class QualityAssessor

      \brief The QualityAssessor class contains all the information needed to perform a
      single sweep over the mesh and make all relevant quality assessments.

      The relevant quality assessments are set by the user or automatically (default)
      by Mesquite when an InstructionQueue object is used. 

      If the mesh has been changed (improved), it is often useful to reuse the same
      QualityAssessor object to reassess the mesh quality.
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
    std::string get_QAFunction_name(enum QualityAssessor::QAFunction);
    
    //! Constructor requires a QualityMetric and an evaluation function
    QualityAssessor(QualityMetric*, enum QAFunction  func,
                    std::string name = "DefaultQualAssessName");

      //!Destructor
    ~QualityAssessor();
    
      //! Provides a name to the QualityAssessor (use it for default name in constructor).
    void set_name(std::string name) { qualityAssessorName = name; };
      //! Retrieves the QualityAssessor name. A default name should be set in the constructor.
    std::string get_name() { return qualityAssessorName; }
    
      //! Adds a quality metric and a wrapper function (min, max, ...).
    void add_quality_assessment(QualityMetric* qm, enum QAFunction  func, MsqError &err);
    
      //! Set the min and max values to be used for the histogram.
    void set_histogram_range(QualityMetric* qm, double min_val, double max_val, MsqError &err);
    
      //! Does one sweep over the mesh and assess the quality with the metrics previously added.
    double assess_mesh_quality(MeshSet &ms, MsqError &err);

      //! Do not print results of assessment.
    void disable_printing_results()
       {
         printingTurnedOff=1;
       }

      /*!Sets the QualityMetric and QAFunction combination that will
        be returned when assess_mesh_quality is called.
      */
    void set_stopping_assessment(QualityMetric* qm, enum QAFunction func,
                                 MsqError &err);
    
  private:
    std::string qualityAssessorName;  
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

      //QAVars is used in assess_mesh_quality.  It gives
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
    };
    
    std::list<Assessor*> assessList;
      //flag to turn off printing
    int printingTurnedOff;
      //pointer to qm used form return value
    QualityMetric* stoppingMetric;
    QAFunction stoppingFunction;
  };

  
} //namespace


#endif // QualityAssessor_hpp
