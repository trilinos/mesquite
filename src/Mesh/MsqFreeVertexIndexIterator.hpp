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

#include "Mesquite.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MesquiteError.hpp"

namespace Mesquite
{

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
      iterOriginator(pd), iterCurrentIndex(0), initialState(true)
    { iterVertexArray = pd->get_vertex_array(err); }
    //! Resets the iterator. 
    //! The next call to next() will set the iterator on the first free vertex. 
    void reset() { initialState=true; iterCurrentIndex=0; }
    //! Increments the iterator. returns false if there is no more free vertex.
    inline bool next();
    //! Returns an index corresponding to a free vertex.
    size_t value() {return iterCurrentIndex;}
  private:
    PatchData* iterOriginator;
    size_t iterCurrentIndex;
    MsqVertex* iterVertexArray;
    bool initialState;
  };
  

#undef __FUNC__
#define __FUNC__ "MsqFreeVertexIndexIterator::next"
  /*! \fn inline bool MsqFreeVertexIndexIterator::next() */
  inline bool MsqFreeVertexIndexIterator::next()
  {
    bool fixed=true;
    while ( fixed ) 
      {
        if ( initialState==true )  initialState=false;
        else  ++iterCurrentIndex;
        
        if ( iterCurrentIndex == iterOriginator->num_vertices() ) {
          return false; 
        }
        fixed = iterVertexArray[iterCurrentIndex].is_flag_set(MsqVertex::MSQ_SOFT_FIXED) ||
          iterVertexArray[iterCurrentIndex].is_flag_set(MsqVertex::MSQ_HARD_FIXED) ;
      }
    return true;
  }
  


} // namespace

#endif //  MsqFreeVertexIndexIterator_hpp
