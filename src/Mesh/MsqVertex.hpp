// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-

//! \file MsqVertex.hpp


#ifndef MSQVERTEX_HPP
#define MSQVERTEX_HPP

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
   class MsqVertex : public Vector3D
   {
   public:
     
     MsqVertex(double x, double y, double z) 
         : Vector3D(x, y, z), vertexBitFlags(0)
       {}
     
     MsqVertex(const Vector3D &vec) 
         : Vector3D(vec), vertexBitFlags(0)
       {}
     
     MsqVertex() 
         : Vector3D(0,0,0), vertexBitFlags(0)
       {}

     void operator=(const Vector3D& rhs)
       { Vector3D::operator=(rhs); }
     
       // This allows for 8 flag bits.
       // I don't think we'll want more than that (yet).
     typedef unsigned char FlagMask;
     
       //! \enum FlagMaskID
       //!   Those are the available flags... currently only return
       //!   is_free.
       //!   Developers: The values used in that enum are used by a bitset,
       //!               so they have to be 2-based (2,4,8,16,32, ...)
     enum FlagMaskID
     {
       MSQ_FREE_VERTEX = 1, //!< vertex is "free")
       MSQ_FLAG_2 = 2,      //!< no other flags at this time
       MSQ_FLAG_3 = 4,
       MSQ_FLAG_4 = 8
     };
     
     int is_free_vertex()
       { return (vertexBitFlags & MSQ_FREE_VERTEX); }
     
     void set_vertex_flag(FlagMask alpha)
       { vertexBitFlags|=alpha; }
     
     void remove_vertex_flag(FlagMask alpha)
       { vertexBitFlags &= (~alpha); }
     
     bool is_flag_set(FlagMask flag)
       { return (vertexBitFlags & flag); }
     
     void move_to_owner()
       { // Not yet written
       }
     
   private:
     FlagMask vertexBitFlags;
   };
  
} //namespace


#endif // MsqVertex_hpp
