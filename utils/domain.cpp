#include "MeshDomain1D.hpp"
#include "PlanarDomain.hpp"
#include "CylinderDomain.hpp"
#include "SphericalDomain.hpp"
#include "DomainClassifier.hpp"
#include "MeshImpl.hpp"
#include "CLArgs.hpp"
#include "MsqError.hpp"

#include "domain.hpp"
#include <iostream>
#include <stdlib.h>

using namespace Mesquite;

class SphereDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  SphereDomainArg( msq_std::vector<MeshDomain*>& domlist,
                   msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool SphereDomainArg::value( const msq_std::vector<double>& list )
{
  double rad = list[0];
  if (rad <= 0.0)
    return false;
  Vector3D center(0,0,0);
  if (list.size() == 4)
    center.set( list[1], list[2], list[3] );
  domList.push_back( new SphericalDomain( center, rad ) );
  dimList.push_back( 2 );
  return true;
}

class CylinderDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  CylinderDomainArg( msq_std::vector<MeshDomain*>& domlist,
                     msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool CylinderDomainArg::value( const msq_std::vector<double>& vals )
{
  double rad = vals[0];
  Vector3D normal( vals[1], vals[2], vals[3] );
  Vector3D point(0,0,0);
  if (vals.size() == 7)
    point.set( vals[4], vals[5], vals[6] );
  domList.push_back( new CylinderDomain( rad, normal, point ) );
  dimList.push_back( 2 );
  return true;
}

class PlanarDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  PlanarDomainArg( msq_std::vector<MeshDomain*>& domlist,
                   msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool PlanarDomainArg::value( const msq_std::vector<double>& list )
{
  Vector3D normal( list[0], list[1], list[2] );
  Vector3D point(0,0,0);
  if (list.size() == 6)
    point.set( list[3], list[4], list[5] );
  domList.push_back( new PlanarDomain( normal, point ) );
  dimList.push_back( 2 );
  return true;
}

class LineDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  LineDomainArg( msq_std::vector<MeshDomain*>& domlist,
                   msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool LineDomainArg::value( const msq_std::vector<double>& vals )
{
  Vector3D dir( vals[0], vals[1], vals[2] );
  Vector3D point(0,0,0);
  if (vals.size() == 6)
    point.set( vals[3], vals[4], vals[5] );
  LineDomain* pdom = new LineDomain( point, dir );
  domList.push_back( pdom );
  dimList.push_back( 1 );
  return true;
}

class CircleDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  CircleDomainArg( msq_std::vector<MeshDomain*>& domlist,
                   msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool CircleDomainArg::value( const msq_std::vector<double>& vals )
{
  double rad = vals[0];
  Vector3D normal( vals[1], vals[2], vals[3] );
  Vector3D point(0,0,0);
  if (vals.size() == 7)
    point.set( vals[4], vals[5], vals[6] );
  CircleDomain* pdom = new CircleDomain( point, normal, rad );
  domList.push_back( pdom );
  dimList.push_back( 1 );
  return true;
}

class PointDomainArg : public CLArgs::DoubleListArgI
{
  private:
  msq_std::vector<MeshDomain*>& domList;
  msq_std::vector<int>& dimList;
  public:
  PointDomainArg( msq_std::vector<MeshDomain*>& domlist,
                  msq_std::vector<int>& dims ) 
    : domList(domlist), dimList( dims ) {}
  virtual bool value( const msq_std::vector<double>& list );
};
bool PointDomainArg::value( const msq_std::vector<double>& vals )
{
  Vector3D point( vals[0], vals[1], vals[2] );
  PointDomain* pdom = new PointDomain( point );
  domList.push_back( pdom );
  dimList.push_back( 0 );
  return true;
}

msq_std::vector<MeshDomain*> domains;
msq_std::vector<int> domain_dims;
SphereDomainArg     sphere_arg( domains, domain_dims );
CylinderDomainArg cylinder_arg( domains, domain_dims );
PlanarDomainArg      plane_arg( domains, domain_dims );
CircleDomainArg     circle_arg( domains, domain_dims );
LineDomainArg         line_arg( domains, domain_dims );
PointDomainArg       point_arg( domains, domain_dims );
CLArgs::ToggleArg skin_mesh( false );

const char* SPHERE_VALUES[] = { "rad", "x", "y", "z" };
const char* CYLINDER_VALUES[] = { "rad", "i", "j", "k", "x", "y", "z" };

void add_domain_args( CLArgs& args )
{
  args.toggle_flag( SKIN_FLAG, "Mark boundary vertices as fixed (default if no domain specified)", &skin_mesh );
  args.double_list_flag( SPHERE_FLAG, "Spherical domain as center and radius", &sphere_arg );
  args.limit_list_flag( SPHERE_FLAG, 4, SPHERE_VALUES );
  args.limit_list_flag( SPHERE_FLAG, 1, SPHERE_VALUES );
  args.double_list_flag( PLANE_FLAG, "Planar domain as normal and point", &plane_arg );
  args.limit_list_flag( PLANE_FLAG, 3, CYLINDER_VALUES+1 );
  args.limit_list_flag( PLANE_FLAG, 6, CYLINDER_VALUES+1 );
  args.double_list_flag( CYLINDER_FLAG, "Cylindrical domain as point, axis, and radius", &cylinder_arg );
  args.limit_list_flag( CYLINDER_FLAG, 4, CYLINDER_VALUES );
  args.limit_list_flag( CYLINDER_FLAG, 7, CYLINDER_VALUES );
  args.double_list_flag( LINE_FLAG, "Linear domain as direction and point", &line_arg );
  args.limit_list_flag( LINE_FLAG, 3, CYLINDER_VALUES+1 );
  args.limit_list_flag( LINE_FLAG, 6, CYLINDER_VALUES+1 );
  args.double_list_flag( CIRCLE_FLAG, "Circular domain as radius, normal, and center", &circle_arg );
  args.limit_list_flag( CIRCLE_FLAG, 4, CYLINDER_VALUES );
  args.limit_list_flag( CIRCLE_FLAG, 7, CYLINDER_VALUES );
  args.double_list_flag( POINT_FLAG, "Point domain", &point_arg );
  args.limit_list_flag( POINT_FLAG, 3, SPHERE_VALUES+1 );
}

MeshDomain* process_domain_args( MeshImpl* mesh )
{
  MsqError err;
  MeshDomain* rval = 0;
  
  if (!domains.empty()) {
    DomainClassifier* result = new DomainClassifier();
    DomainClassifier::classify_geometrically( *result,
                                              mesh,
                                              1e-4,
                                              &domains[0],
                                              &domain_dims[0],
                                              domains.size(),
                                              err );
    rval = result;
  }
  else if (skin_mesh.value()) {
    mesh->mark_skin_fixed( err, false );
  }
  if (err) {
    msq_stdio::cerr << err << msq_stdio::endl;
    exit( 3 );
  }
  
  return rval;
}
