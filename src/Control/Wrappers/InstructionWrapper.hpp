/*! \file InstructionWrapper.hpp

Header file for the Mesquite::InstructionWrapper class

  \author Darryl Melander
  \date   2002-11-14
 */


#ifndef INSTRUCTIONWRAPPER_HPP
#define INSTRUCTIONWRAPPER_HPP

namespace Mesquite
{
  /*! \class InstructionWrapper
    \brief An InstructionWrapper provides a quick, user-friendly way to
           indicate how you want a mesh's quality to be improved.  This
           class is just a placeholder - it doesn't provide anything to
           it's child classes.  At some point in the future it may serve
           a better purpose.
  */
  class InstructionWrapper
  {
      // If wrappers are ever used generically, make this virtual
      // virtual ~InstructionWrapper();
  };
}

#endif
