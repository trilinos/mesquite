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

/*! \file MsqTag.hpp

Header file for the Mesquite::MsqTag class

  \author Thomas Leurent
  \date   2004-03-29
 */


#ifndef MsqTag_hpp
#define MsqTag_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "TargetMatrix.hpp"

namespace Mesquite
{
  
  /*! \class MsqTag
    \brief The MsqTag class manages the storage and access to information associated with
    MsqMeshEntity and MsqVertex instances.
  */
  class MsqTag 
  {
  public:
      //! Constructor initializes all pointers to 0.
    MsqTag() : targets(0), inverseTargets(0), scalars(0),
               numTargets(0), maxNumTargets(0), numScalars(0),
               maxNumScalars(0) {}

      //! copy constructor does a deep copy of all the tag data. 
    MsqTag(const MsqTag& A) { copy(A); }

      //! operator= does a deep copy of all the tag data. 
    MsqTag& operator=(const MsqTag& A)
    { if (&A != this) copy(A); return *this; }

     //! deep copy used in operators.
    void copy(const MsqTag& A) 
    { numScalars = A.numScalars;
      numTargets = A.numTargets;
      maxNumTargets = A.maxNumTargets;
      maxNumScalars = A.maxNumScalars;
      targets = new TargetMatrix[maxNumTargets];
      scalars = new double[maxNumScalars];
      short i;
      for (i=0; i<numTargets; ++i)
        targets[i] = A.targets[i];
      for (i=0; i<numScalars; ++i)
        scalars[i] = A.scalars[i]; }
     
    
      //! 
    ~MsqTag()
    {
      delete [] targets;   
      delete [] scalars;
    }

      //! This sets the targets to an existing array of TargetMatrix s. No allocation is made
      //! by this function.
    void set_targets(TargetMatrix* Ws_pt, short int num_targets, MsqError &err)
    {
      if (targets != 0) {
        MSQ_SETERR(err)( MsqError::INVALID_STATE,
                   "Targets already allocated. Why are you reallocating ?");
        return;
      }

      targets = Ws_pt;
      numTargets = num_targets;
      maxNumTargets = num_targets;
    }

      //! This function allocates target matrices within the tag.
      //! \param num_targets is the number of target matrices to be allocated.
      //! \param alloc_inv if true, allocates space also for the inverse matrices. 
    void allocate_targets(short int num_targets, MsqError &err, bool alloc_inv=false)
    {
      if (targets != 0) {
        MSQ_SETERR(err)( MsqError::INVALID_STATE,
                   "Targets already allocated. Why are you reallocating ?");
        return;
      }

      targets = new TargetMatrix[num_targets];
      if (alloc_inv) inverseTargets = new TargetMatrix[num_targets];
      numTargets = num_targets;
      maxNumTargets = num_targets;
    }
      //!This function can be called when a new tag is being re-used.
      //!  The space needed for targets is reallocated if needed.
    void reallocate_targets(short int num_targets, MsqError &err,
                             bool alloc_inv=false)
    {
      if (num_targets > maxNumTargets) {
        delete[] targets;
        targets = new TargetMatrix[num_targets];      
        if (alloc_inv){
          delete [] inverseTargets;
          inverseTargets = new TargetMatrix[num_targets];
        }
        maxNumTargets = num_targets;
      }
      numTargets = num_targets;
    }
      //! This function can be called when a new tag is created to allocate
      //!  the space needed for scalars.
    void allocate_scalars(short int num_scalars, MsqError &err)
    {
      if (scalars != 0) {
        MSQ_SETERR(err)( MsqError::INVALID_STATE,
                   "Scalars already allocated. Why are you reallocating ?");
        return;
      }

      scalars = new double[num_scalars];
      numScalars = num_scalars;
      maxNumScalars = num_scalars;
    }
      //!This function can be called when a new tag is being re-used.
      //!  The space needed for scalars is reallocated if needed.
    void reallocate_scalars(short int num_scalars, MsqError &err)
    {
      if (num_scalars > maxNumScalars) {
        delete [] scalars;
        scalars = new double[num_scalars];
        maxNumScalars = num_scalars;
      }
      numScalars = num_scalars;
    }
    
      //!
    TargetMatrix* get_targets(const short int &num_targets)
      { assert(num_targets==numTargets); return targets; }

      //!
    TargetMatrix* get_inverse_targets(const short int &num_targets)
      { assert(num_targets==numTargets); return inverseTargets; }

      //!
    TargetMatrix& target_matrix(const short int &i)
      { assert(i<numTargets); return targets[i]; }

      //!
    TargetMatrix& inverse_target_matrix(const short int &i)
      { assert(i<numTargets); return inverseTargets[i]; }

      //!
    double& scalar(const short int &i) { assert(i<numScalars); return scalars[i]; }

  private:
    TargetMatrix* targets;
    TargetMatrix* inverseTargets;
    double* scalars;
      //!number of targets we can access
    short int numTargets;
      //!maximum number of targets for which we have allocated space
    short int maxNumTargets;
      //!number of scalars we can access
    short int numScalars;
      //!maximum number of scalars for which we have allocated space
    short int maxNumScalars;
  };
} //namespace


#endif // MsqTag_hpp
