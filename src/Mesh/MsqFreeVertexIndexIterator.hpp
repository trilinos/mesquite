/*!
  \file   MsqFreeVertexIndexIterator.hpp
  \brief    This file contains the MsqFreeVertexIndexIterator class

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef MsqFreeVertexIndexIterator_hpp
#define MsqFreeVertexIndexIterator_hpp

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "Mesquite.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MesquiteError.hpp"

namespace Mesquite {

  /*! \class MsqFreeVertexIndexIterator
    \brief iterates over indexes of free vetices in a PatchData.

    A free vertex is defined as not having the MSQ_SOFT_FIXED and MSQ_HARD_FIXED
    flags activated.
    
    Use the iterator as follow:
    MsqFreeVertexIndexIterator ind(&patch_data,err);
    ind.reset();
    while (ind.next()) {
      cout << ind.value();
    }  .*/
  class MsqFreeVertexIndexIterator {
  public:
    MsqFreeVertexIndexIterator(PatchData *pd, MsqError &err) :
      originator(pd), current_index(-1)
    { vertex_array = pd->get_vertex_array(err); }
    //! Resets the iterator. 
    //! The next call to next() will set the iterator on the first free vertex. 
    void reset() { current_index=-1; }
    //! Increments the iterator. returns false if there is no more free vertex.
    inline bool next();
    //! Returns an index corresponding to a free vertex.
    int value() {return current_index;}
  private:
    PatchData* originator;
    int current_index;
    MsqVertex* vertex_array;
  };
  

#undef __FUNC__
#define __FUNC__ "MsqFreeVertexIndexIterator::next"
  /*! \fn inline bool MsqFreeVertexIndexIterator::next() */
  inline bool MsqFreeVertexIndexIterator::next()
  {
    bool fixed=true;
    while ( fixed ) 
      {
	++current_index;
	if ( current_index == originator->num_vertices() ) {
	  return false; 
	}
	fixed = vertex_array[current_index].is_flag_set(MsqVertex::MSQ_SOFT_FIXED) ||
        	vertex_array[current_index].is_flag_set(MsqVertex::MSQ_HARD_FIXED) ;
      }
    return true;
  }
  


} // namespace

#endif //  MsqFreeVertexIndexIterator_hpp
