#include "Mesquite.hpp"
#include "SphericalDomain.hpp"
#include "Vector3D.hpp"

// Currently only works if "entity_handle" refers to a vertex
void Mesquite::SphericalDomain::snap_to(Mesh::EntityHandle /*entity_handle*/,
                                        Vector3D &coordinate)
{
  if (!mMesh)
  {
    return;
  }
  
    //translate sphere to origin
  coordinate -= mCenter;
  
    // Move point to the sphere surface
  double len = coordinate.length();
  
    // If it's not right at the center...
  if (len > MSQ_MIN)
    coordinate *= (mRadius / len);
  
    // Move it back off of the origin
  coordinate += mCenter;
}

void Mesquite::SphericalDomain::normal_at(Mesh::EntityHandle /*entity_handle*/,
                                          Vector3D &coordinate)
{
  if (!mMesh)
  {
    coordinate.set(0, 0, 0);
    return;
  }
  
    //translate sphere to origin
  coordinate -= mCenter;
  
    // See how far it is from the center of the sphere.
  double len=coordinate.length();
  if(len<MSQ_MIN)
  {
      // If it's right at the center...
    coordinate.set(0, 0, 0);
  }
  else
  {
      // Otherwise normalize...
    coordinate /= len;
  }
}
