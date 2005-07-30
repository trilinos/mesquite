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
/*!
  \file   PatchDataUser.hpp
  \brief    This file contains the PatchDataUser and the PatchDataParameters classes.

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef PatchDataUser_hpp
#define PatchDataUser_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

#ifndef MSQ_USE_OLD_C_HEADERS
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#else
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#endif


namespace Mesquite
{

  /*! \class PatchDataParameters
   contains all information necessary to fill up a PatchData instance. */
  class PatchDataParameters
  {
  public:
    PatchDataParameters() :
      mType(PatchData::UNDEFINED_PATCH_TYPE),
      mParam1(0),
      mParam2(0)
    {}

    friend class PatchDataUser;
    
    //! Tells the MeshSet what kind of data the patches should include.
    /*! \param patch_type see the PatchData::PatchType enumeration.
      \param patch_param1 meaning depends on patch_type.
      \param patch_param2 meaning depends on patch_type.
    */
    inline void set_patch_type(PatchData::PatchType patch_type, MsqError &err,
                        int patch_param1=0, int patch_param2=0);

    //! Returns Patch Type (local around vertices, local around elements,  global)
    PatchData::PatchType get_patch_type()
    {return mType;}
    //! Returns numbers of layers for local patch. This might not always be a valid measure,
    //! depending on the partition algorythm.
    inline int get_nb_layers();
    
    inline void set_global_patch_type();
    inline void set_element_on_vertex_patch_type( unsigned num_layers );

  private:
    PatchData::PatchType mType; //!< see the enum ... 
    int mParam1, mParam2; //!< For general use in conjunction with PatchType. 
  };


  /*! \class PatchDataUser
    \brief This should be the parent class of all algorithms retrieving information 
    from a MeshSet object. 

    It makes sure that the Patch settings are accessed 
    uniformaly across all those algorithms. Children of PatchDataUser are, 
    among others, thw QualityImprover and QualityAssessor classes. 
    PatchDataUser delegates all settings to its PatchDataParameter member. 

    Alternatively, a PatchDataParameters object can be copied directly (see
    set_all_parameters).
  */
  class PatchDataUser
  {
  protected:
    PatchDataUser() : mParams()
      {}
  public:
    virtual ~PatchDataUser()
      {}

    //! Sets the Patch Type. 
    virtual void set_patch_type(PatchData::PatchType patch_type, MsqError &err,
                                int param1=0, int param2=0) {
      mParams.set_patch_type(patch_type, err, param1, param2); }
    //! Returns the Patch Type.
    PatchData::PatchType get_patch_type() { 
      return mParams.get_patch_type(); }
    //! Returns number of layers (if relevant for partition algorythm). 
    int get_nb_layers() { 
      return mParams.get_nb_layers(); }
    
    /*! Sets all parameters at once by copying a PatchDataParameters object */
    void set_all_parameters(PatchDataParameters &params)
    { mParams = params; }
    //! Returns the PatchDataParameters object.
    PatchDataParameters& get_all_parameters()
    { return mParams; }

    // *** functions related to the algorithms, no to the patch parameters *** 
    
    //! This is the "run" function of PatchDataUser. It can do anything really. 
    virtual double loop_over_mesh(Mesh* ms,                // may NOT be null
                                  MeshDomain* domain,      // may be null
                                  PatchData* global_patch, // may be null
                                  MsqError &err) = 0;
    
    //! Returns the algorithm name
    virtual msq_std::string get_name() = 0;
    enum AlgorithmType { QUALITY_IMPROVER, QUALITY_ASSESSOR, MESH_TRANSFORM, TARGET_CALCULATOR };
    //! Return the algorithm type (to avoid RTTI use). 
    virtual enum AlgorithmType get_algorithm_type() = 0;

  private:
    PatchDataParameters mParams; //!< Contains Patch parameters
  };
  

  /*! \fn PatchDataParameters::set_patch_type(PatchData::PatchType patch_type, MsqError &err, int patch_param1, int patch_param2)
      This function can be over-ridden by the concrete PatchDataUser, in order to cusotomize the Patches available to for the specific algorithm implemented. 

  Utimately, we might want to return an error in the Parent class implementation, 
  in order to force the concrete classes to specify the available types.*/
  inline void PatchDataParameters::set_patch_type(PatchData::PatchType patch_type,
                                                  MsqError &err,
                                                  int patch_param1,
                                                  int patch_param2)
  {
    // For now, no support for VERTICES_ON_ELEMENT_PATCH
    if ( patch_type != PatchData::ELEMENTS_ON_VERTEX_PATCH
	 && patch_type != PatchData::GLOBAL_PATCH )
      {
	MSQ_SETERR(err)("VERTICES_ON_ELEMENT_PATCH not supported yet.",
                        MsqError::NOT_IMPLEMENTED);
	return;
      }
    
    if (patch_type == PatchData::ELEMENTS_ON_VERTEX_PATCH)
      {
        if (patch_param1 < 0)
        {
	  MSQ_SETERR(err)("ELEMENTS_ON_VERTEX_PATCH not supported yet.",
                          MsqError::NOT_IMPLEMENTED);
	  return;
        }
      }
    
    mType = patch_type;
    mParam1 = patch_param1;
    mParam2 = patch_param2;
    
    return;
  }
  
  inline void PatchDataParameters::set_global_patch_type()
  {
    mType = PatchData::GLOBAL_PATCH;
  }
  
  inline void PatchDataParameters::set_element_on_vertex_patch_type( unsigned num_layers )
  {
    mType = PatchData::ELEMENTS_ON_VERTEX_PATCH;
    mParam1 = num_layers;
  }
  
  
  inline int PatchDataParameters::get_nb_layers()
  {
    return mParam1;
  }


} // namespace

#endif //  PatchDataUser_hpp
