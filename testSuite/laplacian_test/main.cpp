/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
  
#include "Mesquite_all_headers.hpp"
#include <ostream>
using namespace Mesquite;
int main(int argc, char* argv[])
{
  Mesquite::MeshImpl mesh;
  MsqPrintError err(std::cout);
  mesh.read_vtk("/home/bktidwell/tmp/meshFiles-new/3D/vtk/prisms/untangled/6-wedge-prism.vtk", err);
  if (err) { std::cout << "read error" << std::endl; return 1; }
  
     //create geometry: plane z=0, normal (0,0,1)
  Vector3D pnt(0,0,0);
  Vector3D s_norm(0,0,1);
  Mesquite::PlanarDomain msq_geom(s_norm, pnt);

    // creates an intruction queue
  InstructionQueue queue1;

    // creates a mean ratio quality metric ...
  ConditionNumberQualityMetric shape_metric;

  QualityAssessor qa=QualityAssessor(&shape_metric);

  queue1.add_quality_assessor(&qa,err);
  if (err) return 1;

  queue1.run_instructions(&mesh, err);
  if (err) return 1;

  return 0;
}
