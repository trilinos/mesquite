// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file InstructionQueue.hpp

Header file for the Mesquite::InstructionQueue class

  \author Thomas Leurent
  \date   2002-05-01
 */


#ifndef InstructionQueue_hpp
#define InstructionQueue_hpp

#include <list>

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "QualityAssessor.hpp"
#include "QualityImprover.hpp"

namespace Mesquite {

 // class QualityAssessor;
//  class QualityImprover;

  /*! \class InstructionQueue
    \brief An InstructionQueue object gathers Mesquite Instructions and ensures
           that the instruction queue is coherent for mesh improvement and/or
           mesh quality assessment purposes.

           The user can instantiate several InstructionQueue objects to be used
           with various MeshSet objects.
           
           The most commonly used functions are:
           -# add_preconditioner(...)
           -# add_quality_assessor(...)
           -# set_master_quality_improver(...)
           -# run_instructions(...)
  */
  class InstructionQueue
  {

    struct QueueEntry
    {
      enum EntryType
        {
          IMPROVER,
          ASSESSOR
        };

      enum EntryType mType;
      union
      {
        QualityImprover* mImprover;
        QualityAssessor* mAssessor;
      };

      QueueEntry(QualityImprover* entry) :
        mType(IMPROVER),
        mImprover(entry)
      {}

      QueueEntry(QualityAssessor* entry) :
        mType(ASSESSOR),
        mAssessor(entry)
      {}
    };

  public:
    InstructionQueue();
    
    void add_preconditioner(QualityImprover* instr, MsqError &err);
    void remove_preconditioner(int index, MsqError &err);
    void insert_preconditioner(QualityImprover* instr, int index, MsqError &err);
    
    void add_quality_assessor(QualityAssessor* instr, MsqError &err);
    void remove_quality_assessor(int index, MsqError &err);
    void insert_quality_assessor(QualityAssessor* instr, int index, MsqError &err);
    
    void set_master_quality_improver(QualityImprover* instr, MsqError &err);
    
    void disable_automatic_quality_assessment()
       { autoQualAssess = false; }
    void enable_automatic_quality_assessment()
       { autoQualAssess = true; }
    
    void run_instructions(MeshSet &msc, MsqError &err);
    void clear();
    
  protected:
    
  private:
    std::list<QueueEntry>::iterator clear_master(MsqError &err);

    std::list<QueueEntry> instructions;

    bool autoQualAssess;
    
    int nbPreConditionners;
    bool isMasterSet;
    int masterInstrIndex; // 0-based
    // keeping an index instead of an iterator in case list is reallocated
  };


} //namespace


#endif // InstructionQueue_hpp
