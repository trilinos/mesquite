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
#include <iostream>
//#include "TSTT_Base.hpp"
#include "MDBMesquiteUtil.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "ConditionNumberQualityMetric.hpp"
#include "LPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "ASMQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "EdgeLengthRangeQualityMetric.hpp"
#include "ConjugateGradient.hpp"
#include "SimplifiedGeometryEngine.hpp"
#include "MsqMessage.hpp"
#include "CompositeOFScalarMultiply.hpp"
#include "CompositeOFMultiply.hpp"
#include "CompositeOFAdd.hpp"
// Function declarations
int run_mesquite(TSTT::Mesh_Handle &mh,
                  TSTT::MeshError *error);
int main(int argc, char** argv)
{
  TSTT::Mesh_Handle my_mesh;
  TSTT::MeshError error=0;
  
  if(argc < 3)
  {
    std::cout << "Usage: tstt_app <inputfile> <outputfile>\n";
    return 1;
  }
  
  init_mesh(&my_mesh, argv[1], &error);

  if(error == 0)
    get_vert_info(my_mesh, &error);

  if(error == 0)
    get_edge_info(my_mesh, &error);
    // get face and region info and tag fixed nodes
  if(error == 0)
    get_face_info(my_mesh, &error);
  TSTT::Int num_regions=0;
  if(error == 0){
    num_regions=get_region_info(my_mesh, &error);
    
      //set the fixed boundary as the nodes of faces, if regions exist
    if(num_regions>0){
      set_fixed_boundary_regions(my_mesh, &error);
    }
      //else set the fixed boundary as the nodes of edges.
    else{
      set_fixed_boundary_faces(my_mesh, &error);
    }
  }
  
    // Smooth the mesh
  if (error == 0)
    error = run_mesquite(my_mesh, &error);

    // Write out the results
  if (error == 0)
    write_results(my_mesh, argv[2], &error);
  
    // clean up
  deinit_mesh(my_mesh);
  return 0;
}

int run_mesquite(TSTT::Mesh_Handle &mh, TSTT::MeshError* /*error*/)
{
  Mesquite::MeshSet mesh_set1;
  Mesquite::MsqError err;
  Mesquite::InstructionQueue queue1;
    //add the given mesh to mesquite's meshset object
  mesh_set1.add_mesh(mh, err); 
  if (err.errorOn) return 1;
    //create geometry
  Mesquite::Vector3D pnt(0, 4, 15);
  Mesquite::Vector3D s_norm(0, 0, 1);
  Mesquite::SimplifiedGeometryEngine msq_geom;
    //ADD GEOMETRY OBJECT IF 2D and NEEDED
  msq_geom.set_geometry_to_plane(s_norm,pnt,err);
  mesh_set1.set_simplified_geometry_engine(&msq_geom);
  
    //**********************************************************
    // creates quality metrics ...    
    //**********************************************************    

   Mesquite::SmoothnessQualityMetric* smooth =
     Mesquite::ASMQualityMetric::create_new();
   
   Mesquite::SmoothnessQualityMetric* edg =
     Mesquite::EdgeLengthQualityMetric::create_new();
   edg->set_averaging_method(Mesquite::QualityMetric::SUM, err);
   if (err.errorOn) return 1;   
   
   Mesquite::SmoothnessQualityMetric* edg_range =
     Mesquite::EdgeLengthRangeQualityMetric::create_new(0,0,err);
   if (err.errorOn) return 1;   
  
   Mesquite::ShapeQualityMetric* cond_no =
     Mesquite::ConditionNumberQualityMetric::create_new();
   
     //**********************************************************
     // creates objective functions using above metrics ...    
     //**********************************************************
    
   Mesquite::LPTemplate* obj_func_cond =
     new Mesquite::LPTemplate(cond_no, 5, err);
   
   Mesquite::LPTemplate* obj_func_smooth =
     new Mesquite::LPTemplate(smooth, 5, err);

   Mesquite::LPTemplate* obj_func_edg =
     new Mesquite::LPTemplate(edg_range, 5, err);

   Mesquite::CompositeOFScalarMultiply* obj_func_cond_scaled =
     new Mesquite::CompositeOFScalarMultiply(1,obj_func_cond);
   
   Mesquite::CompositeOFScalarMultiply* obj_func_smooth_scaled =
     new Mesquite::CompositeOFScalarMultiply(.001,obj_func_smooth);
   
   Mesquite::CompositeOFScalarMultiply* obj_func_edg_scaled =
     new Mesquite::CompositeOFScalarMultiply(1,obj_func_smooth);

   Mesquite::CompositeOFAdd* obj_temp =
     new Mesquite::CompositeOFAdd(obj_func_cond_scaled,obj_func_smooth_scaled);
   
   Mesquite::CompositeOFAdd* obj_final =
     new Mesquite::CompositeOFAdd(obj_temp,obj_func_edg_scaled);

     //specify which type of gradient to use for each of the objective funcs
   obj_func_cond->set_gradient_type(
     Mesquite::ObjectiveFunction::ANALYTICAL_GRADIENT);
   
   obj_func_cond_scaled->set_gradient_type(
     Mesquite::ObjectiveFunction::ANALYTICAL_GRADIENT);
   
   obj_func_smooth->set_gradient_type(
     Mesquite::ObjectiveFunction::ANALYTICAL_GRADIENT);
   
   obj_func_smooth_scaled->set_gradient_type(
     Mesquite::ObjectiveFunction::ANALYTICAL_GRADIENT);

   obj_func_edg->set_gradient_type(
     Mesquite::ObjectiveFunction::ANALYTICAL_GRADIENT);
   
   obj_func_edg_scaled->set_gradient_type(
     Mesquite::ObjectiveFunction::ANALYTICAL_GRADIENT);

   obj_temp->set_gradient_type(
     Mesquite::ObjectiveFunction::ANALYTICAL_GRADIENT);
   
   obj_final->set_gradient_type(
     Mesquite::ObjectiveFunction::ANALYTICAL_GRADIENT);
   
     //**********************************************************       
     // creates quality improver using above objective functions
     //**********************************************************
   
   Mesquite::ConjugateGradient* cg_pass =
     new Mesquite::ConjugateGradient( obj_final, err );
   if (err.errorOn) return 1;   

   cg_pass->set_patch_type( Mesquite::PatchData::GLOBAL_PATCH, err,1 ,1);
   if (err.errorOn) return 1;   
     //cg_pass->set_patch_type(Mesquite::PatchData::ELEMENTS_ON_VERTEX_PATCH,
     //                      err,1 ,1);
   
     //**********************************************************       
     // creates quality assessor (qual_assessor)
     //**********************************************************

   Mesquite::QualityAssessor qual_assessor =
     Mesquite::QualityAssessor(cond_no,
                               Mesquite::QualityAssessor::ALL_MEASURES);
   qual_assessor.add_quality_assessment(edg,Mesquite::QualityAssessor::ALL_MEASURES,err);
   if (err.errorOn) return 1;   
   qual_assessor.add_quality_assessment(edg_range,Mesquite::QualityAssessor::ALL_MEASURES,err);
   if (err.errorOn) return 1;   
   qual_assessor.add_quality_assessment(smooth,Mesquite::QualityAssessor::ALL_MEASURES,err);
   if (err.errorOn) return 1;   
     //**********************************************************       
     // creates inner and outer termination criteria
     //**********************************************************
   Mesquite::TerminationCriterion t1;
   t1.add_criterion_type_with_int(
     Mesquite::TerminationCriterion::NUMBER_OF_ITERATES,60, err);
   if (err.errorOn) return 1;   

   Mesquite::TerminationCriterion t0;
   t0.add_criterion_type_with_int(
     Mesquite::TerminationCriterion::NUMBER_OF_ITERATES,2, err);
   if (err.errorOn) return 1;   

     //add termination criterion to quality improvers
     //if a local scheme do:
     //cg_pass->set_inner_termination_criterion(&t0);
     //cg_pass->set_outer_termination_criterion(&t1);
     //if a global scheme do:
   cg_pass->set_inner_termination_criterion(&t1);

     // sets a culling method on the first QualityImprover
   cg_pass->add_culling_method(Mesquite::PatchData::NO_BOUNDARY_VTX);

     //Conjugate Gradient specific debug flag
   cg_pass->set_debugging_level(0);

     //**********************************************************       
     // Now add objects to the Instruction Queue
     //**********************************************************

   queue1.add_quality_assessor(&qual_assessor,err);
   queue1.set_master_quality_improver(cg_pass, err); 
   if (err.errorOn) return 1;
   queue1.add_quality_assessor(&qual_assessor,err);

   //writeVtkMesh("original_mesh", mesh, err);
   //if (err.errorOn) return 1;
  
     // launches optimization on mesh_set1
   queue1.run_instructions(mesh_set1, err); 
   if (err.errorOn) return 1;

     //PRINT SOME TIMING INFORMATION
   PRINT_TIMING_DIAGNOSTICS();
   return 0;
}
