// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file StoppingCriterion.hpp

Header file for the StoppingCriterion classes.
Includes the base class as well as two composite classes that over-ride
the stop() virtual function.
All three classes can be instantiated.

  \author Thomas Leurent
  \date   2002-05-01
 */


#ifndef StoppingCriterion_hpp
#define StoppingCriterion_hpp

#include <string>

#include "Mesquite.hpp"
#include "MesquiteError.hpp"

namespace Mesquite {
   class QualityAssessor;
   class ObjectiveFunction;
   class MeshSet;
  /*! \class StoppingCriterion

      \brief The StoppingCriterion class contains
  */
  class StoppingCriterion
  {
  public:
    
     /*! */
    enum SCFunction {
       MAX_NODE_MOVEMENT = 1,
       NUMBER_OF_PASSES = 2, 
       OTHER = 4
    };

     struct StoppingCriterionEntry
     {
        enum CriterionEntryType {
           OBJ_FUNC_CRIT,
           ASSESSOR_CRIT,
           SCFUNCTION_CRIT
        };

        enum CriterionEntryType mType;

        union
        {
           QualityAssessor* mAssessor;
           ObjectiveFunction* mObjFunc;
           enum SCFunction mSCFunc;
        };

        // Structure constructors
        StoppingCriterionEntry(ObjectiveFunction* entry) :
           mType(OBJ_FUNC_CRIT),
           mObjFunc(entry)
        {}
        StoppingCriterionEntry(QualityAssessor* entry) :
           mType(ASSESSOR_CRIT),
           mAssessor(entry)
        {}
        StoppingCriterionEntry(enum SCFunction entry) :
           mType(SCFUNCTION_CRIT),
           mSCFunc(entry)
        {}
       StoppingCriterionEntry()
          {}

     };

     
     StoppingCriterion(ObjectiveFunction* of, double lower, double upper);
    StoppingCriterion(QualityAssessor* qm, double lower, double upper);
     StoppingCriterion(enum SCFunction func, double lower, double upper);
    StoppingCriterion(enum SCFunction func, int alpha);

      /*!Resets dynamic information so that a StoppingCriterion may be used
        multiple times.*/
    inline void reset_all(MsqError &err)
       {
         loopCounter=0;
       }
    
     
      //!Destructor
    virtual ~StoppingCriterion(){};

    void increment_counter()
       {
         ++loopCounter;
       }
    void reset_counter()
       {
         loopCounter=0;
       }
    
    
    void set_name(std::string name)
       { SCName = name; };
    
    std::string get_name()
       { return SCName; }

     /*! returns true if stopping criterion is met. */
    virtual bool stop(MeshSet &ms, MsqError &err); 
 protected:
    StoppingCriterion()
       {}
    
  private:
    static int loopCounter;
     std::string SCName;
     StoppingCriterionEntry criterion;
     /* bitset containing embeded functions. */
     long unsigned int standardStoppingCrit;
     double lowerBound;
     double upperBound;
    int mAlpha;
     
  };


   /*! \class CompositeAndStoppingCriterion
       \brief Composes 2 stopping criterion with the AND operator 
   */
   class CompositeAndStoppingCriterion : public StoppingCriterion
   {
   public:
      CompositeAndStoppingCriterion( StoppingCriterion* SC1,
                                     StoppingCriterion* SC2 );
     
      
      virtual bool stop(MeshSet &ms, MsqError &err);
   private:
      StoppingCriterion *stopCrit1;
      StoppingCriterion *stopCrit2;
   };
   

   /*! \class CompositeOrStoppingCriterion
       \brief Composes 2 stopping criterion with the OR operator 
   */
   class CompositeOrStoppingCriterion : public StoppingCriterion
   {
   public:
      CompositeOrStoppingCriterion( StoppingCriterion* SC1,
                                    StoppingCriterion* SC2 );
     
      
      virtual bool stop(MeshSet &ms, MsqError &err);
   private:
      StoppingCriterion *stopCrit1;
      StoppingCriterion *stopCrit2;
   }; 
   
} //namespace


#endif // StoppingCriterion_hpp
