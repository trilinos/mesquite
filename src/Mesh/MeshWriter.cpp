/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2005 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_MESH_WRITER_CPP
#define MSQ_MESH_WRITER_CPP

#include "MeshWriter.hpp"
#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"

#ifdef MSQ_USE_OLD_IO_HEADERS
#include <fstream.h>
#include <string.h>
#include <iomanip.h>
#else
#include <fstream>
#include <string>
#include <iomanip>
#endif

namespace Mesquite {

namespace MeshWriter {

/**\brief Write VTK file
 *
 * Copied from src/Mesh/MeshSet.cpp and adapted for removal of 
 * MeshSet class by J.Kraftcheck on 2005-7-28.
 *
 * This code is provided mainly for debugging.  A more efficient
 * and complete writer implementation is provided in the MeshImpl
 * class for saving meshes that were read from a file initially.
 */
void write_vtk( Mesh* mesh, const char* out_filename, MsqError &err)
{
    // Open the file
  msq_stdio::ofstream file(out_filename);
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }

    // loads a global patch
  PatchData pd;
  pd.set_mesh( mesh );
  pd.fill_global_patch( err ); MSQ_ERRRTN(err);
    
    // Write a header
  file << "# vtk DataFile Version 2.0\n";
  file << "Mesquite Mesh " << out_filename << " .\n";
  file << "ASCII\n";
  file << "DATASET UNSTRUCTURED_GRID\n";
  
    // Write vertex coordinates
  file << "POINTS " << pd.num_nodes() << " float\n";
  size_t i;
  for (i = 0; i < pd.num_nodes(); i++)
  {
    file << pd.vertex_by_index(i)[0] << ' '
         << pd.vertex_by_index(i)[1] << ' '
         << pd.vertex_by_index(i)[2] << '\n';
  }
  
    // Write out the connectivity table
  size_t connectivity_size = 0;
  for (i = 0; i < pd.num_elements(); ++i)
    connectivity_size += pd.element_by_index(i).node_count()+1;
    
  file << "CELLS " << pd.num_elements() << ' ' << connectivity_size << '\n';
  for (i = 0; i < pd.num_elements(); i++)
  {
    msq_std::vector<size_t> vtx_indices;
    pd.element_by_index(i).get_node_indices(vtx_indices);
    file << vtx_indices.size();
    for (msq_stdc::size_t j = 0; j < vtx_indices.size(); ++j)
    {
      file << ' ' << vtx_indices[j];
    }
    file << '\n';
  }
  
    // Write out the element types
  file << "CELL_TYPES " << pd.num_elements() << '\n';
  for (i = 0; i < pd.num_elements(); i++)
  {
    unsigned char type_id = 0;
    switch (pd.element_by_index(i).get_element_type())
    {
      case Mesquite::TRIANGLE:
        type_id = 5;
        break;
      case Mesquite::QUADRILATERAL:
        type_id = 9;
        break;
      case Mesquite::TETRAHEDRON:
        type_id = 10;
        break;
      case Mesquite::HEXAHEDRON:
        type_id = 12;
        break;
    default:
      MSQ_SETERR(err)("element type not implemented",MsqError::NOT_IMPLEMENTED);
      break;
    }
    file << (int)type_id << '\n';
  }
  
    // Write out which points are fixed.
  file << "POINT_DATA " << pd.num_nodes()
       << "\nSCALARS fixed float\nLOOKUP_TABLE default\n";
  for (i = 0; i < pd.num_nodes(); ++i)
  {
    if (pd.vertex_by_index(i).is_free_vertex())
      file << "0\n";
    else
      file << "1\n";
  }
  
    // Close the file
  file.close();
}



/** Writes a gnuplot file directly from the MeshSet.
 *  This means that any mesh imported successfully into Mesquite
 *  can be outputed in gnuplot format.
 *
 *  Within gnuplot, use \b plot 'file1.gpt' w l, 'file2.gpt' w l  
 *   
 *  This is not geared for performance, since it has to load a global Patch from
 *  the mesh to write a mesh file. 
 *
 * Copied from src/Mesh/MeshSet.cpp and adapted for removal of 
 * MeshSet class by J.Kraftcheck on 2005-7-28.
*/
void write_gnuplot( Mesh* mesh, const char* out_filebase, MsqError &err)
{
    // Open the file
  msq_std::string out_filename = out_filebase;
  out_filename += ".gpt";
  msq_stdio::ofstream file(out_filename.c_str());
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }

    // loads a global patch
  PatchData pd;
  pd.set_mesh( mesh );
  pd.fill_global_patch( err ); MSQ_ERRRTN(err);
    
    // Write a header
  file << "\n";
  
  for (size_t i=0; i<pd.num_elements(); ++i)
  {
    msq_std::vector<size_t> vtx_indices;
    pd.element_by_index(i).get_node_indices(vtx_indices);
    for (size_t j = 0; j < vtx_indices.size(); ++j)
    {
      file << pd.vertex_by_index(vtx_indices[j])[0] << ' '
           << pd.vertex_by_index(vtx_indices[j])[1] << ' '
           << pd.vertex_by_index(vtx_indices[j])[2] << '\n';
    }
      file << pd.vertex_by_index(vtx_indices[0])[0] << ' '
           << pd.vertex_by_index(vtx_indices[0])[1] << ' '
           << pd.vertex_by_index(vtx_indices[0])[2] << '\n';
    file << '\n';
  }
  
    // Close the file
  file.close();
}


} // namespace MeshWriter

} // namespace Mesquite

#endif
