/*!
  \file   PatchDataUser.hpp
  \brief    This file contains the PatchDataUser and the PatchDataParameters classes.

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef PatchDataUser_hpp
#define PatchDataUser_hpp

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "Mesquite.hpp"
#include "PatchData.hpp"
#include "MesquiteError.hpp"

namespace Mesquite {

  /*! \class PatchDataParameters
   contains all information necessary to fill up a PatchData instance. */
  class PatchDataParameters {
  public:
    PatchDataParameters() :
      mType(PatchData::UNDEFINED_PATCH_TYPE),
      mParam1(0),
      mParam2(0),
      cullingMethodBits(0)
    {}

    friend class PatchDataUser;
    
    //! Tells the MeshSet what kind of data the patches should include.
    /*! \param patch_type see the PatchData::PatchType enumeration.
      \param patch_param1 meaning depends on patch_type.
      \param patch_param2 meaning depends on patch_type.
    */
    void set_patch_type(PatchData::PatchType patch_type, MsqError &err,
                        int patch_param1=0, int patch_param2=0);

    PatchData::PatchType get_patch_type()
    {return mType;}
    int get_nb_layers(MsqError &err);

    //! Sets on a culling criterion.
    void add_culling_method(enum PatchData::culling_method cm);
    //! No culling performed (sets off all culling criteria).
    void no_culling_method();
    //! Sets off a certain culling criteria. 
    void remove_culling_method(enum PatchData::culling_method cm);
    //! returns the bitset.
    long unsigned int get_culling_method_bits() { return cullingMethodBits; }

  private:
    PatchData::PatchType mType;
    int mParam1, mParam2;
    long unsigned int cullingMethodBits;
  };


  /*! \class PatchDataUser
    \brief This should be the parent class of all algorithms retrieving information 
    from a MeshSet object. 

    It makes sure that the Patch settings are accessed 
    uniformaly across all those algorithms. Children of PatchDataUser are, 
    among others, thw QualityImprover and QualityAssessor classes. 
    PatchDataUser delegates all settings to its PatchDataParameter menber. 

    Alternatively, a PatchDataParameters object can be copied directly (see
    set_all_parameters).*/
  class PatchDataUser {
  protected:
    PatchDataUser() : mParams() { }
  public:
    virtual void set_patch_type(PatchData::PatchType patch_type, MsqError &err,
                                int param1=0, int param2=0) {
      mParams.set_patch_type(patch_type, err, param1, param2); }
    PatchData::PatchType get_patch_type() { 
      return mParams.get_patch_type(); }
    int get_nb_layers(MsqError &err) { 
      return mParams.get_nb_layers(err); }

    void add_culling_method(enum PatchData::culling_method cm) {
      mParams.add_culling_method(cm); }
    void no_culling_method() {
      mParams.no_culling_method(); }
    void remove_culling_method(enum PatchData::culling_method cm) {
      mParams.remove_culling_method(cm); }
    long unsigned int get_culling_method_bits() { 
      return mParams.get_culling_method_bits(); }

    /*! Sets all parameters at once by copying a PatchDataParameters object */
    void set_all_parameters(PatchDataParameters &params)
    { mParams = params; }
    PatchDataParameters& get_all_parameters()
    { return mParams; }

  private:
    PatchDataParameters mParams;
  };
  

#undef __FUNC__
#define __FUNC__ "PatchDataParameter::set_patch_type"
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
	err.set_msg("VERTICES_ON_ELEMENT_PATCH not supported yet.");
	return;
      }
    
    if (patch_type == PatchData::ELEMENTS_ON_VERTEX_PATCH)
      {
        if (patch_param1 < 1)
        {
          err.set_msg("ELEMENTS_ON_VERTEX_PATCH must have the number of layers in patch_param1.");
	return;
        }
      }
    
    mType = patch_type;
    mParam1 = patch_param1;
    mParam2 = patch_param2;
    
    return;
  }
  
  
#undef __FUNC__
#define __FUNC__ "PatchDataParameter::get_nb_layers"
  inline int PatchDataParameters::get_nb_layers(MsqError &err)
  {
    if (mType == PatchData::GLOBAL_PATCH) {
      err.set_msg("Patch Type is GLOBAL_PATCH.");
      return 0;
    }

    return mParam1;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchDataParameter::add_culling_method"
  /*! \fn PatchDataParameter::add_culling_method(enum culling_method cm)
   */
  inline void PatchDataParameters::add_culling_method(enum PatchData::culling_method cm)
  {
    cullingMethodBits |= cm;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchDataParameter::no_culling_method"
  /*! \fn PatchDataParameter::no_culling_method()
   */
  inline void PatchDataParameters::no_culling_method()
  {
    cullingMethodBits = 0;
  }
  
#undef __FUNC__
#define __FUNC__ "PatchDataParameter::remove_culling_method"
  /*! \fn PatchDataParameter::remove_culling_method(enum culling_method cm)
   */
  inline void PatchDataParameters::remove_culling_method(enum PatchData::culling_method cm)
  {
    cullingMethodBits &= ~cm;
  }
  
  
} // namespace

#endif //  PatchDataUser_hpp
