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
#include "MDBInterface.hpp"
#include "MDBRange.hpp"
#include "MDBSkinner.hpp"
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorithms
#include "MeanRatioQualityMetric.hpp"
#include "LPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"

// Function declarations
void set_fixed_boundary(TSTT::Mesh_Handle &mh,
                        TSTT::MeshError *error);
int run_mesquite(TSTT::Mesh_Handle &mh,
                  TSTT::MeshError *error);
void write_results(TSTT::Mesh_Handle &mh,
                   const char* outfile,
                   TSTT::MeshError *error);
void init_mesh(TSTT::Mesh_Handle *mh,
               const char* filename,
               TSTT::MeshError *error);
void get_vert_info(TSTT::cMesh_Handle my_mesh,
                   TSTT::MeshError* error);
void get_edge_info(TSTT::cMesh_Handle my_mesh,
                   TSTT::MeshError* error);
void get_face_info(TSTT::cMesh_Handle my_mesh,
                   TSTT::MeshError* error);
void get_region_info(TSTT::cMesh_Handle my_mesh,
                     TSTT::MeshError* error);
void deinit_mesh(TSTT::Mesh_Handle my_mesh);

// Main program
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

  if(error == 0)
    get_face_info(my_mesh, &error);

  if(error == 0)
    get_region_info(my_mesh, &error);

    // Add the "fixed" tag
  if (error == 0)
    set_fixed_boundary(my_mesh, &error);
  
    // Smooth the mesh
  if (error == 0)
    error = run_mesquite(my_mesh, &error);

    // Write out the results
  if (error == 0)
    write_results(my_mesh, argv[2], &error);
  
    // clean up
  deinit_mesh(my_mesh);
  return !error;
}

void write_results(TSTT::Mesh_Handle &mh,
                   const char* outfile,
                   TSTT::MeshError *error)
{
  MDBInterface* mdb = reinterpret_cast<MDBInterface*>(mh);
  std::vector<MDBEntityHandle> output_list;
  MDBErrorCode err = mdb->write_mesh(outfile, output_list);
  *error = reinterpret_cast<TSTT::MeshError>(err);
}

void init_mesh(TSTT::Mesh_Handle *mh, const char* filename, TSTT::MeshError *error)
{
  TSTT::Mesh_Create(mh, error);
  if(*error != 0)
  {
    std::cout << "Failed to Create Mesh\n";
    return;
  }
  
  TSTT::Mesh_Load(*mh, filename, error);
  if(*error != 0)
  {
    std::cout << "Failed to Load Mesh\n";
    return;
  }
  else
  {
    std::cout << "Successfully read \"" << filename << "\"\n";
  }
}

void get_vert_info(TSTT::cMesh_Handle my_mesh, TSTT::MeshError* error)
{
  int num = 0;
  TSTT::Mesh_Services_GetInt(my_mesh, "number of vertices", &num, error);
  if(*error != 0)
  {
    std::cout << "  Failed to get number of vertices\n";
    return; 
  }
  else
  {
    int num_verts = 0;
    TSTT::Entity_Handle* handles=0;
    double* coords = (double*)malloc(sizeof(double)*3);
    int dum = 3;
    
    std::cout << "  Mesh has " << num << " vertices\n";
    TSTT::Mesh_GetEntities(my_mesh, TSTT::VERTEX, &handles, &num_verts, error);
    assert(num_verts == num);
    TSTT::Entity_GetVertexCoords(my_mesh, (TSTT::cEntity_Handle*)(&handles[0]), 1, TSTT::INTERLEAVED, 
                           &coords, &dum, error);
    std::cout << "    coordinates for vertex " << handles[0] << " are: "
              << coords[0] << ", " << coords[1] << ", " << coords[2] << "\n";

    free(coords);
    TSTT::Mesh_FreeEntityHandles(my_mesh, handles, error);
  }
}


void get_edge_info(TSTT::cMesh_Handle my_mesh, TSTT::MeshError* error)
{
  int num=0;
  TSTT::Mesh_Services_GetInt(my_mesh, "number of edges", &num, error);
  if(*error != 0)
  {
    std::cout << "  Failed to get number of edges\n";
    return;
  }
  else
  {
    std::cout << "  Mesh has " << num << " edges\n";
  }
}

void get_face_info(TSTT::cMesh_Handle my_mesh, TSTT::MeshError* error)
{
    
  int num=0;
  TSTT::Mesh_Services_GetInt(my_mesh, "number of faces", &num, error);
  if(*error != 0)
  {
    std::cout << "  Failed to get number of faces\n";
    return; 
  }
  else
  {

    std::cout << "  Mesh has " << num << " faces\n";
    if(num > 0)
    {
      TSTT::Entity_Handle *adj_ents=0;
      TSTT::Int* csr_pointer=0;
      TSTT::Int* csr_data=0;
      TSTT::Int num_adj=0;
      int i, j, k;
      TSTT::Entity_Handle *handles = 0;
      TSTT::Int num_faces=0;
      TSTT::EntityTopology* top=0;
      int num_tris=0, num_quads=0;
      
      TSTT::Mesh_GetEntities(my_mesh, TSTT::FACE, &handles, &num_faces, error);
      top = (TSTT::EntityTopology*) malloc(sizeof(TSTT::EntityTopology)*num_faces);
      TSTT::Entity_GetTopology(my_mesh, (TSTT::cEntity_Handle*)handles,
                         num_faces, top, error);
      for(i=0; i<num_faces; i++)
      {
        if(top[i] == TSTT::TRIANGLE)
          num_tris++;
        else if(top[i] == TSTT::QUADRILATERAL)
          num_quads++;
      }
      std::cout << "    " << num_tris << " of which are triangles\n";
      std::cout << "    " << num_quads << " of which are quadrilaterals\n";
      free(top);
      
      TSTT::Entity_GetAdjacencies(my_mesh, (TSTT::cEntity_Handle*)handles, num_faces,
          TSTT::VERTEX, &adj_ents, &csr_pointer,
          &csr_data, &num_adj, error);
      if(*error == 0)
      {
        std::cout << "      " << num_faces << " faces share "
             << num_adj << " vertices\n";
        std::cout << "      the last face, " << handles[num_faces] << " has "
             << csr_pointer[num_faces] - csr_pointer[num_faces-1]
             << " vertices adjacent to it\n";
        std::cout << "        vertices are -- ";
        k = csr_pointer[num_faces] - csr_pointer[num_faces-1];
        for(j=0; j<k; j++)
          std::cout << adj_ents[csr_data[csr_pointer[num_faces-1]+j]]
               <<" ";
        std::cout << "\n";
      }
      
      TSTT::Mesh_FreeEntityAdjacency(my_mesh, adj_ents,
                               csr_pointer, csr_data, error);
      TSTT::Mesh_FreeEntityHandles(my_mesh, handles, error);
    }
    
  }
}


void get_region_info(TSTT::cMesh_Handle my_mesh, TSTT::MeshError* error)
{
  int num=0;
  TSTT::Mesh_Services_GetInt(my_mesh, "number of regions", &num, error);
  if(*error != 0)
  {
    std::cout << "  Failed to get number of regions\n";
    return; 
  }
  else
  {
    TSTT::Entity_Handle * handles = 0;
    enum TSTT::EntityTopology* top=0;
    TSTT::cEntity_Handle* c_handles;
    TSTT::Int num_regions=0;
    int num_tets=0, num_hexes=0;
    int i;
    
    std::cout << "  Mesh has " << num << " regions\n";
    TSTT::Mesh_GetEntities(my_mesh, TSTT::REGION, &handles, &num_regions, error);
    top = (TSTT::EntityTopology*) malloc(sizeof(TSTT::EntityTopology)*num_regions);
    c_handles = (TSTT::cEntity_Handle*)handles;
    TSTT::Entity_GetTopology(my_mesh, c_handles, num_regions, top, error);
    for(i=0; i<num_regions; i++)
    {
      if(top[i] == TSTT::TETRAHEDRON)
        num_tets++;
      else if(top[i] == TSTT::HEXAHEDRON)
        num_hexes++;
    }
    std::cout << "    " << num_tets  << " of which are tetrahedrons\n";
    std::cout << "    " << num_hexes << " of which are hexahedrons\n";
    if(num_regions > num_hexes+num_tets)
      std::cout << "    " << num_regions-(num_hexes+num_tets)
           << " of which have another entity topology\n";
    free(top);
    TSTT::Mesh_FreeEntityHandles(my_mesh, handles, error);
  }
}


void deinit_mesh(TSTT::Mesh_Handle my_mesh)
{
  TSTT::MeshError error;
  TSTT::Mesh_Destroy(my_mesh, &error);
  if(error != 0)
    std::cout << "Failed to shutdown properly\n";
}

int run_mesquite(TSTT::Mesh_Handle &mh, TSTT::MeshError *error)
{
    // Initialize a MeshSet object
  Mesquite::MeshSet mesh_set1;
  Mesquite::MsqPrintError err;
  mesh_set1.add_mesh(mh, err); 
  if (err) return 1;
  
    // Create an intruction queue
  Mesquite::InstructionQueue queue1;

    // Create a mean ratio quality metric
  Mesquite::ShapeQualityMetric* mean_ratio = Mesquite::MeanRatioQualityMetric::create_new();
  
    // Build an objective function with it
  Mesquite::LPTemplate* obj_func = new Mesquite::LPTemplate(mean_ratio, 2, err);
  if (err) return 1;
  
    // Create the steepest descent  optimization procedures
  Mesquite::SteepestDescent* pass1 = new Mesquite::SteepestDescent( obj_func );
  Mesquite::TerminationCriterion tc2;
  tc2.add_criterion_type_with_int(Mesquite::TerminationCriterion::NUMBER_OF_ITERATES,10,err);
  if (err) return 1;
    //Mesquite::StoppingCriterion sc2(Mesquite::StoppingCriterion::NUMBER_OF_PASSES, 10);
  pass1->set_outer_termination_criterion(&tc2);
  pass1->add_culling_method(Mesquite::PatchData::NO_BOUNDARY_VTX);
  
  queue1.set_master_quality_improver(pass1, err); 
  if (err) return 1;
  queue1.run_instructions(mesh_set1, err); 
  if (err) return 1;
}

void set_fixed_boundary(TSTT::Mesh_Handle &mh, TSTT::MeshError *error)
{
  MDBInterface* mdb = reinterpret_cast<MDBInterface*>(mh);
  
    // Create the "fixed" tag with a value of 0
  MDBTag fixed_tag_handle = 0;
  int val = 0;
  mdb->tag_create_dense("fixed", sizeof(int),
                        fixed_tag_handle,
                        reinterpret_cast<void*>(&val));
  
    // Get all the quads/tris
  MDBRange quad_range;
  mdb->get_entities_by_dimension(2, quad_range);
  
    // Skin the collection of quads/tris
  MDBRange forward_edge_range, reverse_edge_range;
  MDBSkinner skinner(mdb);
  skinner.find_skin(quad_range, forward_edge_range, reverse_edge_range);

  std::cout << "There are " << forward_edge_range.size()
            << " forward edges in the model" << std::endl;
  std::cout << "There are " << reverse_edge_range.size()
            << " reverse edges in the model" << std::endl;

    // Examine each edge, mark its nodes as fixed
  std::vector<MDBEntityHandle> nodes;
  for (MDBRange::const_iterator edge_iter = forward_edge_range.begin();
       edge_iter != forward_edge_range.end();
       ++edge_iter)
  {
    mdb->get_adjacencies(*edge_iter, 0, true, nodes);
    for (std::vector<MDBEntityHandle>::const_iterator node_iter = nodes.begin();
         node_iter != nodes.end();
         ++node_iter)
    {
      int tag_val = 1;
      mdb->tag_set_data(fixed_tag_handle, *node_iter,
                        reinterpret_cast<void*>(&tag_val));
    }
  }

    // Now just for fun, let's see how many tagged nodes we get
  int tval = 1;
  std::vector<MDBEntityHandle> tagged_verts;
  mdb->get_entities_with_tag_value(MDBVertex, fixed_tag_handle,
                                   &tval, tagged_verts);
  std::cout << "Number of tagged boundary vertices: " << tagged_verts.size()
            << std::endl;
}
