#include "PlanarDomain.hpp"

void Mesquite::PlanarDomain::snap_to(Mesquite::Mesh::EntityHandle entity_handle,
                                     Vector3D &coordinate)
{
  if (!mMesh)
    return;
  
    //translate plane/point to pass through origin
  coordinate -= mPoint;
  
  double len=(coordinate%mNormal);
  coordinate += mPoint-(len*mNormal);
}


void Mesquite::PlanarDomain::normal_at(
  Mesquite::Mesh::EntityHandle /*entity_handle*/,
  Vector3D &coordinate)
{
  coordinate = mNormal;
}
