#ifndef SIMPLIFIED_GEOMETRY_ENGINE_HPP
#define SIMPLIFIED_GEOMETRY_ENGINE_HPP
/*!\file SimplifiedGeometryEngine.hpp
  Limited geometry capabilities.
  \author Michael Brewer          
  \date   2003-01-27    
*/
#include <vector>
#include "Vector3D.hpp"
#include "MsqVertex.hpp"
namespace Mesquite
{
  /*!\class SimplifiedGeometryEngine
        
    \brief  Contains some simple geometry capabilities.
    
    The SimplifiedGeometryEngine class contains limited geometry
    functionality.  It designed to be able to perform two geometry
    calls needed within Mesquite on a PLANE or a SPHERE.  These two
    functions are snap_vertex_on_domain and get_surface_normal.
  */
  class SimplifiedGeometryEngine
  {
      
  public:
      /*! Simplified Geometry has two geometrical entities implemented:
        PLANE and SPHERE.  The third option, NONE, is the default value
        in SimplifiedGeeometryEngine's constructor for geomType.*/
    enum MSQ_GEOM_TYPE{
      NONE,
      PLANE,
      SPHERE
    };
      //!SimplifiedGeometryEngine defualts GeomType to be 'NONE'.
    SimplifiedGeometryEngine();

      //!Set the SimplifiedGeometryEngine to be a plane with surface normal
      //!'norm' and containing the point 'point'.       
    void set_geometry_to_plane(Vector3D norm, Vector3D point, MsqError &err);
      //!Set the SimplifiedGeoemtryEngine to be a sperhe with radius 'rad' and
      //!center 'center'.
    void set_geometry_to_sphere(Vector3D center, double rad, MsqError &err);
    
      //!destructor
    ~SimplifiedGeometryEngine()
      {}
      //!Calculates the surface normal and returns it in 'surf_norm'.
    void get_surface_normal(MsqVertex* vert, Vector3D &surf_norm,
                            MsqError &err);

      //!Moves the vertex pointed to by 'vert' to the surface.
    void snap_vertex_to_domain(MsqVertex* vert, MsqError &err);
    
  private:
    MSQ_GEOM_TYPE geomType;
    Vector3D geomVec;
    Vector3D geomPoint;
    double sphereRad;
  };
  
}

#endif
