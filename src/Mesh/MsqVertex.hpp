// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-

/*! \file MsqVertex.hpp
  \brief Mesquite's vertex object.
*/
#ifndef MSQVERTEX_HPP
#define MSQVERTEX_HPP

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
    /*!
      \class MsqVertex
      \brief MsqVertex is the Mesquite object that stores information about
      the vertices in the mesh.*/
  class MsqVertex : public Vector3D
   {
   public:
       //!Construct vertex using three doubles.
     MsqVertex(double x, double y, double z) 
         : Vector3D(x, y, z), vertexBitFlags(0)
       {}
       //!Construct vertex using Vector3D.
     MsqVertex(const Vector3D &vec) 
         : Vector3D(vec), vertexBitFlags(0)
       {}
       //!Construct default vertex with coordinates (0.0,0.0,0.0)
     MsqVertex() 
         : Vector3D(0,0,0), vertexBitFlags(0)
       {}
     
     void operator=(const Vector3D& rhs)
       { Vector3D::operator=(rhs);
         vertexBitFlags = 0; }
     
     void operator=(const MsqVertex& rhs)
       { Vector3D::operator=(rhs);
         vertexBitFlags = rhs.vertexBitFlags; }
     
       // This allows for 8 flag bits.
       // I don't think we'll want more than that (yet).
     typedef char FlagMask;
     
       //! \enum FlagMaskID
       //!   Those are the available flags... currently only return
       //!   is_free.
       //!   Developers: The values used in that enum are used by a bitset,
       //!               so they have to be 2-based (2,4,8,16,32, ...)
     enum FlagMaskID
     {
       MSQ_NO_VTX_FLAG = 0,
       MSQ_ALGO_FLAG0 = 1, //!< vertex is "free"
       MSQ_SOFT_FIXED = 2,  //!< vertex is fixed. This flag can be set on and off. 
       MSQ_HARD_FIXED = 4,  //!< vertex is always fixed. This can only be set on and never off.
       MSQ_COORDS_CHANGED = 8,
       MSQ_FLAG_3 = 16,
       MSQ_FLAG_4 = 32,
       MSQ_ALGO_FLAG1 = 64, //!< free bit, to be used by algorithm if needed.
       MSQ_ALGO_FLAG2 = 128 //!< free bit, to be used by algorithm if needed. 
     };
       //!Returns true if vertex is ``free''.
     bool is_free_vertex()
       { return ( !bool(vertexBitFlags & MSQ_SOFT_FIXED) &&
                  !bool(vertexBitFlags & MSQ_HARD_FIXED) ); }
     
     void set_soft_fixed_flag()
       { vertexBitFlags|=MSQ_SOFT_FIXED; }
     
     void remove_soft_fixed_flag()
       { vertexBitFlags &= (~MSQ_SOFT_FIXED); }
     
     void set_hard_fixed_flag()
       { vertexBitFlags|=MSQ_HARD_FIXED; }
     
     void set_vertex_flag(FlagMaskID flag)
       { vertexBitFlags|=flag; }
     
     void remove_vertex_flag(FlagMaskID flag)
       { vertexBitFlags &= (~flag); }
     
     bool is_flag_set(FlagMaskID flag)
       { return bool(vertexBitFlags & flag); }
     
     void move_to_owner()
       { // Not yet written
       }
     
   private:
     FlagMask vertexBitFlags;
   };
  
} //namespace


#endif // MsqVertex_hpp
