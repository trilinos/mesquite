#include "SimplifiedGeometryEngine.hpp"
#include "MsqMessage.hpp"
using namespace Mesquite;
#undef __FUNC__
#define __ FUNC__ "SimplifiedGeometryEngine::SimplifiedGeometryEngine()"
  /*!  
  Mesquite engine for simple geometry calculations.  Currently, this
  object is designed to provide "surface normal" and "move to surface"
  functionality for spheres and planes.  The geometry is defined by
  the data passed to the constructor.  Planes are defined by the
  surface normal and a point in the plane.  Spheres are defined by
  the radius vector and the ceter point.
 */
SimplifiedGeometryEngine::SimplifiedGeometryEngine()
    :geomType(NONE)
{
  
}


//!Set the SimplifiedGeometryEngine to be a plane with surface normal
//!'norm' and containing the point 'point'.
void SimplifiedGeometryEngine::set_geometry_to_plane(Vector3D norm,
                                                     Vector3D point,
                                                     MsqError &err)
{
  geomType=PLANE;
  geomPoint=point;
  geomVec=norm;
}

//!Set the SimplifiedGeoemtryEngine to be a sperhe with radius 'rad' and
//!center 'center'.
void SimplifiedGeometryEngine::set_geometry_to_sphere(Vector3D center,
                                                      double rad,
                                                      MsqError &err)
{
  geomType=SPHERE;
  geomPoint=center;
  sphereRad=rad;
  if(sphereRad<MSQ_MIN){
    err.set_msg("Sphere's radius is either too small or negative");
  }
}


void SimplifiedGeometryEngine::get_surface_normal(MsqVertex* vert,
                                                  Vector3D &surf_norm,
                                                  MsqError &err)
{

  Vector3D temp_vec;
  double len;
  switch(geomType){
    case PLANE:
        //geomVec holds planes normal
      surf_norm=geomVec;
      break;
    case SPHERE:
        //translate sphere and point to origin
      temp_vec[0]=vert->x()-geomPoint[0];
      temp_vec[1]=vert->y()-geomPoint[1];
      temp_vec[2]=vert->z()-geomPoint[2];
      len=temp_vec.length();
      if(len<MSQ_MIN){
        err.set_msg("Vert lies at center of sphere");
      }
      else{
        surf_norm=(temp_vec/len);
      }
      break;
    default:
      err.set_msg("get_surface_normal called with geometry type NONE.");
  }       
}


void SimplifiedGeometryEngine::snap_vertex_to_domain(MsqVertex* vert,
                                                     MsqError &err)
{
  Vector3D temp_vec;
  double len;
  switch(geomType){
    case PLANE:
        //translate plane/point to pass through origin
      temp_vec[0]=vert->x()-geomPoint[0];
      temp_vec[1]=vert->y()-geomPoint[1];
      temp_vec[2]=vert->z()-geomPoint[2];
        //len=1;//temp_vec.length_squared();
      len=((temp_vec%geomVec));
      *vert=temp_vec-(len*geomVec)+geomPoint;
        //As a debug step, make sure step was 'correct'
      MSQ_DEBUG_ACTION(1,{
        if(fabs((*vert-geomPoint)%geomVec)>
           fabs(((temp_vec+(len*geomVec))%geomVec))){
          PRINT_INFO("POINT NOT IN PLANE x %f, y %f, z %f",vert->x(),vert->y(),vert->z());
          PRINT_INFO("\nFIRST dot  = %e",(*vert-geomPoint)%geomVec);
          PRINT_INFO("\nSECOND dot = %e",((temp_vec+(len*geomVec))%geomVec));
          err.set_msg("Vertex moved out of the plane.");
        }
      });
      
      break;
    case SPHERE:
        //translate sphere and point to origin
      temp_vec[0]=vert->x()-geomPoint[0];
      temp_vec[1]=vert->y()-geomPoint[1];
      temp_vec[2]=vert->z()-geomPoint[2];
      temp_vec/=temp_vec.length();
      temp_vec*=sphereRad;
      vert->x(temp_vec[0]+geomPoint[0]);
      vert->y(temp_vec[1]+geomPoint[1]);
      vert->z(temp_vec[2]+geomPoint[2]);
      break;
    default:
      err.set_msg("snap_vertex_on_domain called with geometry type NONE.");
  }
}



