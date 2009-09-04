
#ifndef Mesquite_SmartLaplacianSmoother_hpp 
#define Mesquite_SmartLaplacianSmoother_hpp

#include "Mesquite.hpp"
#include "RelaxationSmoother.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
#endif

namespace MESQUITE_NS
{
  /*\brief Do laplacian smooth, but don't invert elements.
   */  
  class MESQUITE_EXPORT SmartLaplacianSmoother : public RelaxationSmoother 
  {
  public:
    /**
     *\param OF ObjectiveFunction used by some termination criteria
     */
    SmartLaplacianSmoother( ObjectiveFunction* OF = NULL ) 
      : RelaxationSmoother(OF) {}
    
    ~SmartLaplacianSmoother();
    virtual msq_std::string get_name() const;
    
    static size_t num_inverted( PatchData& pd, MsqError& err );
    
  protected:
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);

  private:
    msq_std::vector<size_t> adjVtxList;    
  };

  

  
}

#endif
