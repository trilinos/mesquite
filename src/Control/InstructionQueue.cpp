// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file InstructionQueue.cpp

Member functions of the Mesquite::InstructionQueue class

  \author Thomas Leurent
  \date   2002-05-01
 */

#include <string>

#include "InstructionQueue.hpp"
#include "MeshSet.hpp"

using namespace Mesquite;


InstructionQueue::InstructionQueue() :
  autoQualAssess(true),
  nbPreConditionners(0),
  isMasterSet(false),
  masterInstrIndex(-1)
{
}

#undef __FUNC__
#define __FUNC__ "InstructionQueue::add_preconditioner"
/*! \fn InstructionQueue::add_preconditioner(QualityImprover* instr, MsqError &err)
    \brief adds a QualityImprover at the end of the instruction list

    This function cannot be used once the set_master_quality_improver()
    function has been used.
    
    See also insert_preconditioner().
  */
void InstructionQueue::add_preconditioner(QualityImprover* instr,
                                        MsqError &err)
{
  if (isMasterSet) {
    err.set_msg("cannot add preconditionners once the master QualityImprover has been set.");
    return;
  }
  
  // creepy shallow copy
  //  QualityImprover* instr_copy = instr->clone();  
  // QueueEntry entry(instr_copy);
  
  QueueEntry entry(instr);
  instructions.push_back(entry);
  nbPreConditionners++;
}


#undef __FUNC__
#define __FUNC__ "InstructionQueue::remove_preconditioner"
/*! \fn InstructionQueue::remove_preconditioner(int index, MsqError &err)
    \brief removes a QualityImprover* from the instruction queue

    \param index is 0-based. An error is set if the index does not correspond
           to a valid element in the queue.
*/
void InstructionQueue::remove_preconditioner(int index, MsqError &err)
{
  // checks index is valid
  if ( index == masterInstrIndex ) {
    err.set_msg("cannot remove master QualityImprover.");
    return;
  } else if (index >= instructions.size() ) {
    err.set_msg("index points beyond end of list.");
    return;
  }
  
  // position the instruction iterator over the preconditionner to delete
  std::list<QueueEntry>::iterator pos;
  pos = instructions.begin();
  std::advance(pos, index);
  if ( pos->mType == QueueEntry::IMPROVER ) {
    std::string name = pos->mImprover->get_name();
    std::cout << "  o InstructionQueue: removing QualityImprover " << name <<  ".\n"; 
    instructions.erase(pos);
    nbPreConditionners--;
  }
  else
    err.set_msg("index does not point to a QualityImprover.");
}  


#undef __FUNC__
#define __FUNC__ "InstructionQueue::insert_preconditioner"
/*! \fn InstructionQueue::insert_preconditioner(QualityImprover* instr, int index, MsqError &err)
    \brief inserts a QualityImprover* into the instruction queue.

    Pre-conditionners can only be inserted before the master QualityImprover.

    \param index is 0-based. An error is set if the index does not correspond
           to a valid position in the queue.
*/
void InstructionQueue::insert_preconditioner(QualityImprover* instr,
                                           int index, MsqError &err)
{
  // checks index is valid
  if (isMasterSet==true && index > masterInstrIndex) {
    err.set_msg("cannot add a preconditionner after the master QualityImprover.");
    return;
  }
  if (index >= instructions.size() ) {
    err.set_msg("index points beyond end of list.");
    return;
  }

  // position the instruction iterator
  std::list<QueueEntry>::iterator pos;
  pos = instructions.begin();
  std::advance(pos, index);
  // adds the preconditioner
  QueueEntry entry(instr);
  instructions.insert(pos,entry);
  nbPreConditionners++;
}


#undef __FUNC__
#define __FUNC__ "InstructionQueue::add_quality_assessor"
/*! \fn InstructionQueue::add_quality_assessor(QualityAssessor* instr, MsqError &err)
    \brief adds a QualityAssessor to the instruction queue.

    QualityAssessor pointers can be added at any time to the instruction queue.
*/
void InstructionQueue::add_quality_assessor(QualityAssessor* instr,
                                          MsqError &err)
{
  // creepy shallow copy
  //  QualityAssessor* instr_copy = instr->clone(); 
  // QueueEntry entry(instr_copy);
  
  QueueEntry entry(instr);
  instructions.push_back(entry);
}


#undef __FUNC__
#define __FUNC__ "InstructionQueue::remove_quality_assessor"
/*! \fn InstructionQueue::remove_quality_assessor(int index, MsqError &err)
    \brief removes a QualityAssessor* from the instruction queue

    \param index is 0-based. An error is set if the index does not correspond
           to a valid element in the queue.
*/
void InstructionQueue::remove_quality_assessor(int index, MsqError &err)
{
  // checks index is valid
  if (index >= instructions.size() ) {
    err.set_msg("index points beyond end of list.");
    return;
  }
  
  // position the instruction iterator over the QualityAssessor to delete
  std::list<QueueEntry>::iterator pos;
  pos = instructions.begin();
  std::advance(pos, index);
  if ( pos->mType == QueueEntry::ASSESSOR ) {
    std::string name = pos->mAssessor->get_name();
    std::cout << "  o InstructionQueue: removing QualityAssessor " << name << ".\n"; 
    instructions.erase(pos);
  }
  else
    err.set_msg("index does not point to a QualityAssessor.");
}  


#undef __FUNC__
#define __FUNC__ "InstructionQueue::insert_quality_assessor"
/*! \fn InstructionQueue::insert_quality_assessor(QualityAssessor* instr, int index, MsqError &err)
    \brief inserts a QualityAssessor* into the instruction queue.

    QualityAssessors can be inserted at any position in the instruction queue.

    \param index is 0-based. An error is set if the index is past the end of the queue.
*/
void InstructionQueue::insert_quality_assessor(QualityAssessor* instr,
                                           int index, MsqError &err)
{
  // checks index is valid
  if (index > instructions.size()) {
    err.set_msg("index points two positions beyond end of list.");
    return;
  }

  // position the instruction iterator
  std::list<QueueEntry>::iterator pos;
  pos = instructions.begin();
  std::advance(pos, index);
  // adds the QualityAssessor
  QueueEntry entry(instr);
  instructions.insert(pos,entry);
}


#undef __FUNC__
#define __FUNC__ "InstructionQueue::set_master_quality_improver"
void InstructionQueue::set_master_quality_improver(QualityImprover* instr,
                                                 MsqError &err)
{
  // creepy shallow copy
  //  QualityImprover* instr_copy = instr->clone();  
  // masterInstruction = instr_copy;

  if (isMasterSet) {
    std::cout << "WARNING: InstructionQueue::set_master_quality_improver():\n"
         << "\tOverwriting previously specified master quality improver.\n";
    // if master is already set, clears it and insert the new one at the same position.
    std::list<QueueEntry>::iterator master_pos;
    master_pos = this->clear_master(err); MSQ_CHKERR(err);
    QueueEntry entry(instr);
    instructions.insert(master_pos, entry);
    isMasterSet = true;
  } else {
    // if master is not set, add it at the end of the queue.
    QueueEntry entry(instr);
    instructions.push_back(entry);
    isMasterSet = true;
    masterInstrIndex = instructions.size()-1;
  }
}

#undef __FUNC__
#define __FUNC__ "InstructionQueue::run_instructions"
void InstructionQueue::run_instructions(MeshSet &ms, MsqError &err)
{   
  if (nbPreConditionners != 0 && isMasterSet == false ) {
    err.set_msg("no pre-conditionners allowed if master QualityImprover is not set.");
    return;
  }
  
  std::list<QueueEntry>::const_iterator instr_iter;
    //Michael
    //std::cout<<"\nFirst check of time "<<err.since_last_check();
  // For each pass QualityImprover/QualityAssessor in the preconditionner list
  for (instr_iter = instructions.begin();
       instr_iter != instructions.end(); ++instr_iter) {
    // applies the QualityImprover/QualityAssessor to the MeshSet
    if (instr_iter->mType == QueueEntry::IMPROVER) {
      instr_iter->mImprover->loop_over_mesh(ms, err); MSQ_CHKERR(err); }
    else if (instr_iter->mType == QueueEntry::ASSESSOR) {
      instr_iter->mAssessor->assess_mesh_quality(ms, err); MSQ_CHKERR(err); }
    else
      err.set_msg("Unknown instruction type.");
      //std::cout<<"\nInstruction queue time "<<err.since_last_check();
  }
    //std::cout<<"\nApproximate TOTAL time "<<err.since_birth();
}

  
#undef __FUNC__
#define __FUNC__ "InstructionQueue::clear"
void InstructionQueue::clear()
{
  instructions.clear();
  autoQualAssess = true;
  isMasterSet = false;
  masterInstrIndex = -1;
}


#undef __FUNC__
#define __FUNC__ "InstructionQueue::clear_master"
std::list<InstructionQueue::QueueEntry>::iterator InstructionQueue::clear_master(MsqError &err)
{
  std::list<QueueEntry>::iterator instr_iter;
  std::list<QueueEntry>::iterator master_pos;
  
  if (!isMasterSet) {
    err.set_msg("No master quality improver to clear.");
    return instr_iter;
  }
  
    // position the instruction iterator over the master quality improver
  master_pos = instructions.begin();
  std::advance(master_pos, masterInstrIndex);
  
    // erases the master quality improver
  instr_iter = instructions.erase(master_pos);
  isMasterSet = false;
  
    // returns the position where the Master was
  return instr_iter;
}
