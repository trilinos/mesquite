// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE:  7-Nov-02 at 16:22:26
//  LAST-MOD:  8-Nov-02 at 10:27:00 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file QualityImprover.cpp

Implements a couple of default virtual functions of the virtual class

 */
// DESCRIP-END.
//

#include "QualityImprover.hpp"

using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "QualityImprover::set_patch_type"
/*! \fn QualityImprover::set_patch_type(MeshSet::PatchType type, MsqError &err)

    By default, this virtual function sets an error flag.
    It must be explicitly overriden in the concrete class in order to be available,
    since only the concrete algorythm class knows the types of Patch it can use.
 */
void QualityImprover::set_patch_type(MeshSet::PatchType type, MsqError &err)
{
  err.set_msg("This algorythm does not allow to redefine the Patch type.");
}
