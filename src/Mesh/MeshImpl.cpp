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
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 16-May-02 at 10:26:21
//  LAST-MOD: 18-Jun-04 at 16:07:03 by Thomas Leurent
//
/*! \file MeshImpl.cpp

\brief This files contains a mesh database implementation that can be used
to run mesquite by default.
  
    \author Thomas Leurent
    \author Darryl Melander
    \date 2003-05-16  
 */

#include "MeshImpl.hpp"
#include "MsqMessage.hpp"

#include "Vector3D.hpp"

#ifdef USE_STD_INCLUDES
#include <fstream>
#include <string>
#include <iomanip>
#else
#include <fstream.h>
#include <string.h>
#include <iomanip.h>
#endif

#ifdef MSQ_USING_EXODUS
#include "exodusII.h"
#endif

MSQ_USE(ifstream);
MSQ_USE(ofstream);
MSQ_USE(setprecision);
MSQ_USE(string);
MSQ_USE(cerr);

namespace Mesquite
{
  template<typename X> class MeshImpl_EntityIterator : public Mesquite::VertexIterator
  {
  public:
    MeshImpl_EntityIterator(X* array,
                            size_t element_count) 
        : mArray(array),
          mCount(element_count),
          mIndex(0)
      {}
    
      // Moves the iterator back to the first
      // entity in the list.
    virtual void restart()
      {
        mIndex = 0;
      }
    
      // *iterator.  Return the handle currently
      // being pointed at by the iterator.
    virtual Mesquite::Mesh::EntityHandle operator*() const
      {
        return is_at_end() ? NULL : reinterpret_cast<void*>(mArray+mIndex);
      }
    
      // ++iterator
    virtual void operator++()
      {
        ++mIndex;
      }
    
      // Returns false until the iterator has
      // been advanced PAST the last entity.
      // Once is_at_end() returns true, *iterator
      // returns NULL.
    virtual bool is_at_end() const
      {
        return (mIndex < mCount ? false : true);
      }
    
  private:
    X *mArray;
    size_t mCount;
    size_t mIndex;
  };
}

Mesquite::MeshImpl::MeshImpl() 
    : vertexArray(NULL),
      elementArray (NULL),
      vertexCount(0),
      elementCount(0),
      onBoundaryBits(NULL),
      vertexMesquiteByte(NULL),
      newVertIndices(NULL),
      v2eOffset(NULL),
      totalVertexUses(0),
      v2E(NULL),
      numCoords(3)
{}

Mesquite::MeshImpl::~MeshImpl() 
{
  delete [] vertexArray;
  delete [] elementArray;
  delete [] onBoundaryBits;
  delete [] vertexMesquiteByte;
  delete [] newVertIndices;
  delete [] v2eOffset;
  delete [] v2E;
}

void Mesquite::MeshImpl::read_vtk(const char* in_filename,
                                  Mesquite::MsqError &err)
{
    // Open the file
  ifstream ifs(in_filename);
  if (!ifs)
  {
    string err_msg("file ");
    err_msg += in_filename;
    err_msg += " not found";
    err.set_msg(err_msg);
    return;
  }
  
    // Skip the first 4 lines
  char line[256];
  size_t i;
  for (i = 0; i < 4; ++i)
    ifs.getline(line,256);
  
    // get the number of vertices.
  char word[81];
  ifs >> word;
  if (strcmp(word, "POINTS") != 0)
  {
    string err_msg = "MeshImpl::read_vtk: expecting word POINTS, not ";
    err_msg += word;
    err_msg += ".";
    err.set_msg(err_msg);
    return;
  }
  ifs >> vertexCount;
  ifs.getline(line,256);
  
    // Create an array big enough for all the vertices
  vertexArray = new Mesquite::MeshImpl::Vertex[vertexCount];
  
    // Also create a cache for elements_get_attached_vertices().
  newVertIndices = new size_t[vertexCount];
  memset(newVertIndices, 0, vertexCount*sizeof(size_t));
  
    // Get the coordinates for each vertex
  double* coord_ptr = reinterpret_cast<double*>(vertexArray);
  for(i = 0; i < vertexCount; ++i)
  {
    ifs >> *coord_ptr >> *(coord_ptr + 1) >> *(coord_ptr + 2);
    coord_ptr += 3;
  }
  
    // get the number of regions
  ifs >> word;
  if (strcmp(word, "CELLS") != 0)
  {
    string err_msg = "MeshImpl::read_vtk: expecting word CELLS, not ";
    (err_msg += word) += ".";
    err.set_msg(err_msg);
  }
  
  size_t table_size;
  ifs >> elementCount >> table_size;
  
    // Read the connectivity table
  size_t *connectivity_table = new size_t[table_size];
  size_t *cur_entry = connectivity_table;
  
  for (i = 0; i < table_size; ++i)
  {
    ifs >> *(cur_entry++);
  }
  
    // get the number of CELL_TYPES identifier
  ifs >> word;
  if (strcmp(word, "CELL_TYPES") != 0)
  {
    string err_msg = "MeshImpl::read_vtk: expecting word CELL_TYPES, not ";
    (err_msg += word) += ".";
    err.set_msg(err_msg);
    return;
  }
  ifs >> i;
  if (i != elementCount)
  {
    err.set_msg("MeshImpl::read_vtk: nb of CELL_TYPES != nb of CELLS.");
  }
  
    // Store each cell type
  vector<size_t> cell_type;
  cell_type.reserve(elementCount);
  for (i=0; i < elementCount; ++i)
    ifs >> cell_type[i];
  
    // Create the array of elements
  elementArray = new MeshImpl::Element[elementCount];
  
    // Create an element for each cell
  totalVertexUses = 0;
  cur_entry = connectivity_table;
  for (i=0; i<elementCount; ++i)
  {
    int j;

    switch (cell_type[i])
    {
      case 5:  // Triangle
        if (*cur_entry++ != 3)
        {
          err.set_msg("MeshImpl::read_vtk: expecting 3 vtx for a triangle.");
          return;
        }
        elementArray[i].mType = Mesquite::TRIANGLE;
        totalVertexUses += 3;
        for (j = 0; j < 3; j++)
        {
          elementArray[i].vertexIndices[j] = *cur_entry++;
        }
        break;
      case 9:  // Quad
        if (*cur_entry++ != 4)
        {
          err.set_msg("MeshImpl::read_vtk: expecting 4 vtx for a Quad.");
          return;
        }
        elementArray[i].mType = Mesquite::QUADRILATERAL;
        totalVertexUses += 4;
        for (j = 0; j < 4; j++)
        {
          elementArray[i].vertexIndices[j] = *cur_entry++;
        }
        break;
      case 10: // Tet
        if (*cur_entry++ != 4)
        {
          err.set_msg("MeshImpl::read_vtk: expecting 4 vtx for a Tet.");
          return;
        }
        elementArray[i].mType = Mesquite::TETRAHEDRON;
        totalVertexUses += 4;
        for (j = 0; j < 4; j++)
        {
          elementArray[i].vertexIndices[j] = *cur_entry++;
        }
        break;
      case 12: // Hex 
        if (*cur_entry++ != 8)
        {
          cerr << "MeshImpl::read_vtk: expecting 8 vtx for an Hex.\n";
          return;
        }
        elementArray[i].mType = Mesquite::HEXAHEDRON;
        totalVertexUses += 8;
        for (j = 0; j < 8; j++)
        {
          elementArray[i].vertexIndices[j] = *cur_entry++;
        }
        break;
    }
  }
  delete [] connectivity_table;
  
    // Create a set of bits to store whether it's on the boundary or not
  size_t num_bitset_bytes = (vertexCount / 8);
  num_bitset_bytes += (vertexCount % 8) ? 1 : 0;
  onBoundaryBits = new unsigned char[num_bitset_bytes];
  memset(onBoundaryBits, 0, sizeof(unsigned char)*num_bitset_bytes);
  
    // Get the number of BOUNDARY_POINTS,
    // Set them only if they are present.
  ifs >> word;
  if (strcmp(word, "POINT_DATA") == 0)
  {
    int num_boundary_verts;
    ifs >> num_boundary_verts;
    
      // Skip the next 5 words
    for (i = 1; i <= 5; i++)
      ifs >> word;
    
      // See if this vertex is marked as a boundary vertex
    for (int vert_index=0; vert_index < num_boundary_verts; vert_index++)
    {
      int val;
      ifs >> val;
      if (val != 0)
      {
        unsigned char bit_flag = 1 << (vert_index%8);
        onBoundaryBits[vert_index / 8] |= bit_flag;
      }
    }
  }

    // Create a Mesquite byte for each vertex
  vertexMesquiteByte = new unsigned char[vertexCount];
  memset(vertexMesquiteByte, 0, sizeof(unsigned char)*vertexCount);
  
  ifs.close();
}

void Mesquite::MeshImpl::write_vtk(const char* out_filebase,
                                   Mesquite::MsqError &err)
{
    // Open the file
  string out_filename = out_filebase;
  out_filename += ".vtk";
  ofstream file(out_filename.c_str());
  if (!file)
  {
    err.set_msg("Unable to open file");
    return;
  }
  
    // Write a header
  file << "# vtk DataFile Version 2.0\n";
  file << "Mesquite Mesh " << out_filebase << " .\n";
  file << "ASCII\n";
  file << "DATASET UNSTRUCTURED_GRID\n";
  
    // Write vertex coordinates
  file << "POINTS " << vertexCount << " float\n";
  size_t i;
  for (i = 0; i < vertexCount; i++)
  {
      //MB: Is there a way to use setprecision (or an equivalent)
      //when use_std_includes is not defined.
#ifdef USE_STD_INCLUDES   
    file <<setprecision(15)<< vertexArray[i].coords[0] << ' '
         <<setprecision(15)<< vertexArray[i].coords[1] << ' '
         <<setprecision(15)<< vertexArray[i].coords[2] << '\n';
#else
    file << vertexArray[i].coords[0] << ' '
         << vertexArray[i].coords[1] << ' '
         << vertexArray[i].coords[2] << '\n';
#endif
  }
  
    // Write out the connectivity table
  size_t connectivity_size = elementCount + totalVertexUses;
  file << "CELLS " << elementCount << ' ' << connectivity_size << '\n';
  for (i = 0; i < elementCount; i++)
  {
    size_t verts_this_elem = vertices_in_topology(elementArray[i].mType);
    file << verts_this_elem;
    for (size_t j = 0; j < verts_this_elem; j++)
    {
      file << ' ' << elementArray[i].vertexIndices[j];
    }
    file << '\n';
  }
  
    // Write out the element types
  file << "CELL_TYPES " << elementCount << '\n';
  for (i = 0; i < elementCount; i++)
  {
    unsigned char type_id = 0;
    switch (elementArray[i].mType)
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
      err.set_msg("element type not implemented");
      break;
    }
    file << (int)type_id << '\n';
  }
  
    // Write out which points are fixed.
  file << "POINT_DATA " << vertexCount
       << "\nSCALARS fixed float\nLOOKUP_TABLE default\n";
  for (i = 0; i < vertexCount; i++)
  {
    if (onBoundaryBits[i/8] & ((unsigned char)(1) << i % 8))
      file << "1\n";
    else
      file << "0\n";
  }
  
    // Close the file
  file.close();
}

void Mesquite::MeshImpl::read_exodus(const char*
#ifdef MSQ_USING_EXODUS
                                     in_filename
#endif                                     
                                     , Mesquite::MsqError &err)
{
#ifndef MSQ_USING_EXODUS
  err.set_msg("Exodus not enabled in this build of Mesquite");
  return;
#else
  if (vertexArray != NULL)
  {
    err.set_msg("Attempting to read second file into a MeshImpl");
    return;
  }
  
  int app_float_size = sizeof(double);
  int file_float_size = 0;
  float exo_version = 0;
  int exo_err = 0;
  
    // Open the file
  int file_id = ex_open(in_filename, EX_READ, &app_float_size,
                    &file_float_size, &exo_version);

    // Make sure we opened the file correctly
  if (file_id < 0)
  {
    err.set_msg("Unable to open file");
    return;
  }
  
    // make sure the file is saved as doubles
  if (file_float_size != sizeof(double))
  {
    err.set_msg("File saved with float-sized reals.  Can only read files "
                "saved with doubles.");
    return;
  }

  char title[MAX_LINE_LENGTH];
  int dim, vert_count, elem_count, block_count, ns_count, ss_count;
  
    // get info about the file
  exo_err = ex_get_init(file_id, title, &dim, &vert_count,
                        &elem_count, &block_count, &ns_count, &ss_count);
  if (exo_err < 0)
  {
    err.set_msg("Unable to get entity counts from file.");
    return;
  }
  
  vertexCount = vert_count;
  elementCount = elem_count;
  
    // Now that we know how big our arrays have to be, allocate...
    // ...vertices
  vertexArray = new Mesquite::MeshImpl::Vertex[vertexCount];
    // ...cache for elements_get_attached_vertices()
  newVertIndices = new size_t[vertexCount];
  memset(newVertIndices, 0, vertexCount*sizeof(size_t));
    // ...elements
  elementArray = new MeshImpl::Element[elementCount];
    // ...bits indicating what's on the boundary
  size_t num_bitset_bytes = (vertexCount / 8);
  num_bitset_bytes += (vertexCount % 8) ? 1 : 0;
  onBoundaryBits = new unsigned char[num_bitset_bytes];
  memset(onBoundaryBits, 0, sizeof(unsigned char)*num_bitset_bytes);
    // ...a "Mesquite" byte for each vertex
  vertexMesquiteByte = new unsigned char[vertexCount];
  memset(vertexMesquiteByte, 0, sizeof(unsigned char)*vertexCount);
  
    // Now fill in the data
  
    // Get the vertex coordinates
  double* temp_doubles;
  if (dim == 2)
  {
    numCoords = 2;
    temp_doubles = new double[vertexCount*2];
    exo_err = ex_get_coord(file_id,
                           reinterpret_cast<void*>(temp_doubles),
                           reinterpret_cast<void*>(temp_doubles + vertexCount),
                           NULL);
  }
  else
  {
    numCoords = 3;
    temp_doubles = new double[vertexCount*3];
    exo_err = ex_get_coord(file_id,
                           reinterpret_cast<void*>(temp_doubles),
                           reinterpret_cast<void*>(temp_doubles + vertexCount),
                           reinterpret_cast<void*>(temp_doubles +
                                                   vertexCount + vertexCount));
  }
    // Make sure it worked
  if (exo_err < 0)
  {
    err.set_msg("Unable to retrieve vertex coordinates from file.");
    return;
  }
    // Stuff coordinates into vertexArray
  double *cur_coord = temp_doubles;
  size_t i, j, k;
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < vertexCount; j++)
    {
      vertexArray[j].coords[i] = *cur_coord;
      cur_coord++;
    }
  }
  if (dim == 2)
  {
    for (j = 0; j < vertexCount; j++)
      vertexArray[j].coords[2] = 0.0;
  }
  else
  {
    for (j = 0; j < vertexCount; j++)
    {
      vertexArray[j].coords[2] = *cur_coord;
      cur_coord++;
    }
  }
  delete [] temp_doubles;

    // Process elements block by block
  int *block_ids = new int[block_count];
  exo_err = ex_get_elem_blk_ids(file_id, block_ids);
  if (exo_err < 0)
  {
    err.set_msg("Unable to read block IDs from file.");
    delete [] block_ids;
    return;
  }
  size_t conn_table_size = 0;
  int *connectivity_table = NULL;
  MeshImpl::Element *cur_elem = elementArray;
  totalVertexUses=0;
  for (i = 0; i < block_count; i++)
  {
      // Get info about this block's elements
    char elem_type_str[MAX_STR_LENGTH];
    int num_block_elems, verts_per_elem, num_atts;
    exo_err = ex_get_elem_block(file_id, block_ids[i], elem_type_str,
                                &num_block_elems, &verts_per_elem,
                                &num_atts);
    totalVertexUses+=(num_block_elems*verts_per_elem);
    if (exo_err < 0)
    {
      err.set_msg("Unable to read parameters for block.");
      delete [] block_ids;
      return;
    }
    
      // Figure out which type of element we're working with
    Mesquite::EntityTopology elem_type;
    for (j = 0; j < 3; j++)
      elem_type_str[j] = toupper(elem_type_str[j]);
    if (!strncmp(elem_type_str, "TRI", 3))
    {
      elem_type = Mesquite::TRIANGLE;
    }
    else if (!strncmp(elem_type_str, "QUA", 3) ||
             !strncmp(elem_type_str, "SHE", 3))
    {
      elem_type = Mesquite::QUADRILATERAL;
    }
    else if (!strncmp(elem_type_str, "HEX", 3))
    {
      elem_type = Mesquite::HEXAHEDRON;
    }
    else if (!strncmp(elem_type_str, "TET", 3))
    {
      elem_type = Mesquite::TETRAHEDRON;
    }
    else
    {
      err.set_msg("Unrecognized element type in block");
      continue;
    }
    
      // Get the connectivity
    if (conn_table_size < num_block_elems*verts_per_elem)
    {
      conn_table_size = num_block_elems*verts_per_elem;
      delete [] connectivity_table;
      connectivity_table = new int[conn_table_size];
    }
    exo_err = ex_get_elem_conn(file_id, block_ids[i], connectivity_table);
    if (exo_err < 0)
    {
      err.set_msg("Unable to read element block connectivity.");
      delete [] block_ids;
      delete [] connectivity_table;
      return;
    }

      // Put the connectivity into the elementArray
    int *cur_entry = connectivity_table;
    for (j = 0; j < num_block_elems; j++)
    {
      cur_elem->mType = elem_type;
      for (k = 0; k < verts_per_elem; k++)
      {
        cur_elem->vertexIndices[k] = (*cur_entry)-1;
        cur_entry++;
      }
      cur_elem++;
    }
  }
  delete [] block_ids;
  delete [] connectivity_table;
  
    // Finally, mark boundary nodes
  int num_fixed_nodes=0;
  int num_dist_in_set=0;
  if(ns_count>0){
    exo_err=ex_get_node_set_param(file_id,111,&num_fixed_nodes,
                                  &num_dist_in_set);
    if(exo_err<0){
      Message::print_warning("\nError opening nodeset 111, no boundary nodes marked.");
      num_fixed_nodes=0;
    }
  }
  int *fixed_nodes =NULL;
  if(num_fixed_nodes>0)
    fixed_nodes = new int[num_fixed_nodes];
  exo_err = ex_get_node_set(file_id, 111, fixed_nodes);
  if(exo_err<0){
    err.set_msg("Error retrieving fixed nodes.");
  }
  
    // See if this vertex is marked as a boundary vertex
  for (i=0; i < num_fixed_nodes; ++i)
  {
    fixed_nodes[i]-=1;
    unsigned char bit_flag = 1 << (fixed_nodes[i]%8);
    onBoundaryBits[fixed_nodes[i] / 8] |= bit_flag;  
  }
  if(fixed_nodes!=NULL)
    delete [] fixed_nodes;
  exo_err=ex_close(file_id);
  if(exo_err<0)
    err.set_msg("Error closing Exodus file.");
#endif
}
//!Writes an exodus file of the mesh.
void Mesquite::MeshImpl::write_exodus(const char*
#ifdef MSQ_USING_EXODUS
                                       out_filename
#endif                                     
                                       , Mesquite::MsqError &err)
{
    //just return an error if we don't have access to exodus
#ifndef MSQ_USING_EXODUS
  err.set_msg("Exodus not enabled in this build of Mesquite");
  return;
#else
  size_t i, j, counter;
  if (vertexArray == NULL)
  {
    err.set_msg("No vertices in MeshImpl.  Nothing written to file.");
    return;
  }
    //get some element info
    //We need to know how many of each element type we have.  We
    //are going to create an element block for each element type
    //that exists the mesh.  Block 1 will be tri3; block 2 will be
    //shell; block 3 will be tetra, and block 4 will be hex.
  int num_tri=0;
  int num_quad=0;
  int num_tet=0;
  int num_hex=0;
  int block_count=0;//one block for each element type
    //count each element 
  for(i=0;i<elementCount;++i){
    if(elementArray[i].mType==Mesquite::TRIANGLE){
      ++num_tri;
    }
    else if (elementArray[i].mType==Mesquite::QUADRILATERAL){
      ++num_quad;
    }
    else if(elementArray[i].mType==Mesquite::TETRAHEDRON){
      ++num_tet;
    }
    else if (elementArray[i].mType==Mesquite::HEXAHEDRON){
      ++num_hex;
    }
    else
      err.set_msg("Unrecognized element type");
  }
    //if an element of the specific type exists, we need a block for it
  if(num_tri>0)
    ++block_count;
  if(num_quad>0)
    ++block_count;
  if(num_tet>0)
    ++block_count;
  if(num_hex>0)
    ++block_count;
  if(block_count>4){
    err.set_msg("Too many element types in file");
  }
  if(block_count<1){
    err.set_msg("Too few element types in file");
  }
    //figure out if we have fixed nodes, if so, we need a nodeset
  int num_fixed_nodes=0;
  for (i = 0; i < vertexCount; i++)
  {
    if (onBoundaryBits[i/8] & ((unsigned char)(1) << i % 8))
      ++num_fixed_nodes; 
  }
    //write doubles instead of floats
  int app_float_size = sizeof(double);
  int file_float_size = sizeof(double);
  int exo_err = 0;
  
    // Create the file.  If it exists, clobber it.  This could be dangerous.
  int file_id = ex_create(out_filename, EX_CLOBBER, &app_float_size,
                          &file_float_size);

    // Make sure we opened the file correctly
  if (file_id < 0)
  {
    err.set_msg("Unable to create file");
    return;
  }
  
  char title[MAX_LINE_LENGTH]="Mesquite Generated Exodus File";
  int dim=3;
  int vert_count=get_total_vertex_count(err);
  int elem_count=get_total_element_count(err);
  
  int ns_count=0;
  if(num_fixed_nodes>0)
    ns_count=1;
  int ss_count=0;
  
    // put the initial info about the file
  exo_err = ex_put_init(file_id, title, dim, vert_count,
                        elem_count, block_count, ns_count, ss_count);
  if (exo_err < 0)
  {
    err.set_msg("Unable to initialize file data.");
    return;
  }
    //array of nodal coords.
  double* temp_doubles = new double[vertexCount*3];
 
  counter=0;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < vertexCount; j++)
    {
      temp_doubles[counter]=vertexArray[j].coords[i];
      ++counter;
    }
  }
  if(counter!= (3*vertexCount))
    err.set_msg("Counter at incorrect number.");
  
    //put the coords
  exo_err = ex_put_coord(file_id,
                         reinterpret_cast<void*>(temp_doubles),
                         reinterpret_cast<void*>(temp_doubles + vertexCount),
                         reinterpret_cast<void*>(temp_doubles +
                                                 vertexCount + vertexCount));
  
    // Make sure it worked
  if (exo_err < 0)
  {
    err.set_msg("Unable to put vertex coordinates in file.");
    return;
  } 
  delete [] temp_doubles;
    //put the names of the coordinates
  char *coord_names[3];
  char x_co[2]="x";
  char y_co[2]="y";
  char z_co[2]="z";
  coord_names[0]=x_co;
  coord_names[1]=y_co;
  coord_names[2]=z_co;
  exo_err = ex_put_coord_names(file_id, coord_names);
  
    // Process elements in a block
  int block_ids[4];
  block_ids[0]=1;//tri
  block_ids[1]=2;//quad
  block_ids[2]=3;//tet
  block_ids[3]=4;//hex
  
    //double check things
  if ((num_tri + num_quad + num_tet + num_hex) != elementCount){
    err.set_msg("Error in determining the element types.");
  }
  int num_atts=0;

    //for each element type that exists, set up the block
  if(num_tri>0)
    exo_err = ex_put_elem_block(file_id, block_ids[0], "TRI3",
                                num_tri, 3, num_atts);
  if(exo_err<0)
    err.set_msg("Error creating the tri block.");
  if(num_quad>0)
    exo_err = ex_put_elem_block(file_id, block_ids[1], "SHELL",
                                num_quad, 4, num_atts);
  if(exo_err<0)
    err.set_msg("Error creating the quad block.");
  if(num_tet>0)
    exo_err = ex_put_elem_block(file_id, block_ids[2], "TETRA",
                                num_tet, 4, num_atts);
  if(exo_err<0)
    err.set_msg("Error creating the tet block.");
  if(num_hex>0)
    exo_err = ex_put_elem_block(file_id, block_ids[3], "HEX",
                                num_hex, 8, num_atts);
  if(exo_err<0)
    err.set_msg("Error creating the hex block.");
    //alloc space for the connectivity arrays
  int *tri_connectivity = new int[3*num_tri];
  int *quad_connectivity = new int[4*num_quad];
  int *tet_connectivity = new int[4*num_tet];
  int *hex_connectivity = new int[8*num_hex];
    //counters for the different element types
  int tri_counter=0;
  int quad_counter=0;
  int tet_counter=0;
  int hex_counter=0;
    //put each element in the appropriate connectivity array
  for(i=0;i<elementCount;++i){
    if(elementArray[i].mType==Mesquite::TRIANGLE){
      tri_connectivity[3*tri_counter]=elementArray[i].vertexIndices[0]+1;
      tri_connectivity[1+(3*tri_counter)]=elementArray[i].vertexIndices[1]+1;
      tri_connectivity[2+(3*tri_counter)]=elementArray[i].vertexIndices[2]+1;
      ++tri_counter;
    }
    else if (elementArray[i].mType==Mesquite::QUADRILATERAL){
      quad_connectivity[4*quad_counter]=elementArray[i].vertexIndices[0]+1;
      quad_connectivity[1+(4*quad_counter)]=elementArray[i].vertexIndices[1]+1;
      quad_connectivity[2+(4*quad_counter)]=elementArray[i].vertexIndices[2]+1;
      quad_connectivity[3+(4*quad_counter)]=elementArray[i].vertexIndices[3]+1;
      ++quad_counter;
    }
    else if(elementArray[i].mType==Mesquite::TETRAHEDRON){
      tet_connectivity[4*tet_counter]=elementArray[i].vertexIndices[0]+1;
      tet_connectivity[1+(4*tet_counter)]=elementArray[i].vertexIndices[1]+1;
      tet_connectivity[2+(4*tet_counter)]=elementArray[i].vertexIndices[2]+1;
      tet_connectivity[3+(4*tet_counter)]=elementArray[i].vertexIndices[3]+1;
      ++tet_counter;
    }
    else if (elementArray[i].mType==Mesquite::HEXAHEDRON){
      hex_connectivity[8*hex_counter]=elementArray[i].vertexIndices[0]+1;
      hex_connectivity[1+(8*hex_counter)]=elementArray[i].vertexIndices[1]+1;
      hex_connectivity[2+(8*hex_counter)]=elementArray[i].vertexIndices[2]+1;
      hex_connectivity[3+(8*hex_counter)]=elementArray[i].vertexIndices[3]+1;
      hex_connectivity[4+(8*hex_counter)]=elementArray[i].vertexIndices[4]+1;
      hex_connectivity[5+(8*hex_counter)]=elementArray[i].vertexIndices[5]+1;
      hex_connectivity[6+(8*hex_counter)]=elementArray[i].vertexIndices[6]+1;
      hex_connectivity[7+(8*hex_counter)]=elementArray[i].vertexIndices[7]+1;
      ++hex_counter;
    }
    else
      err.set_msg("Unrecognized element type again");
  }
    //double check number of each element type
  if(tri_counter != num_tri)
    err.set_msg("Tri numbers not consistent");
  if(quad_counter != num_quad)
    err.set_msg("Quad numbers not consistent");
  if(tet_counter != num_tet)
    err.set_msg("Tet numbers not consistent");
  if(hex_counter != num_hex)
    err.set_msg("Hex numbers not consistent");
  
  //for each element type that exists, put the connectivity with the block
  if(num_tri>0)
    exo_err = ex_put_elem_conn(file_id, block_ids[0], tri_connectivity);
  if(exo_err<0)
    err.set_msg("Error in putting tri block");
  if(num_quad>0)
    exo_err = ex_put_elem_conn(file_id, block_ids[1], quad_connectivity);
  if(exo_err<0)
    err.set_msg("Error in putting quad block");
  if(num_tet>0)
    exo_err = ex_put_elem_conn(file_id, block_ids[2], tet_connectivity);
  if(exo_err<0)
    err.set_msg("Error in putting tet block");
  if(num_hex>0)
    exo_err = ex_put_elem_conn(file_id, block_ids[3], hex_connectivity);
  if(exo_err<0)
    err.set_msg("Error in putting hex block");
    //delete connectivity arrays.
  delete [] tri_connectivity;
  delete [] quad_connectivity;
  delete [] tet_connectivity;
  delete [] hex_connectivity;
  
    // Finally, mark boundary nodes
  
  if(num_fixed_nodes>0){
    exo_err=ex_put_node_set_param(file_id, 111, num_fixed_nodes, 0);
    if(exo_err<0)
      err.set_msg("Error while initializing node set.");
    int *fixed_nodes= new int[num_fixed_nodes];
    int fixed_node_counter=0;
    for (i = 0; i < vertexCount; i++)
    {
      if (onBoundaryBits[i/8] & ((unsigned char)(1) << i % 8)){
        fixed_nodes[fixed_node_counter]=i+1;
        ++fixed_node_counter;
      }
      
    }
    exo_err=ex_put_node_set(file_id, 111, fixed_nodes);
    if(exo_err<0)
      err.set_msg("Error while writing node set.");
  }
  exo_err=ex_close(file_id);
  if(exo_err<0)
    err.set_msg("Error closing Exodus file.");
  
#endif
}   
// Returns whether this mesh lies in a 2D or 3D coordinate system.
int Mesquite::MeshImpl::get_geometric_dimension(MsqError &/*err*/) const
{
  return numCoords;
}
    
// Returns the number of entities of the indicated type.
size_t Mesquite::MeshImpl::get_total_vertex_count(MsqError &/*err*/) const
{
  return vertexCount;
}
size_t Mesquite::MeshImpl::get_total_element_count(MsqError &/*err*/) const
{
  return elementCount;
}
    
// Fills array with handles to all vertices in the mesh.
#undef __FUNC__
#define __FUNC__ "MeshImpl::get_all_vertices"
void Mesquite::MeshImpl::get_all_vertices(
  Mesquite::Mesh::VertexHandle *vert_array,
  size_t array_size, MsqError &err)
{
  if (array_size > vertexCount)
    array_size = vertexCount;
  if (array_size < vertexCount) {
    err.set_msg("Array of insufficient size. "
                "Returning incomplete vertex list");
  }

  for (size_t i = 0; i < array_size; i++)
    vert_array[i] = vertexArray + i;
}

// Fills array with handles to all elements in the mesh.
void Mesquite::MeshImpl::get_all_elements(
  Mesquite::Mesh::ElementHandle *elem_array,
  size_t array_size, MsqError &/*err*/)
{
  if (array_size > elementCount)
    array_size = elementCount;
  
  for (size_t i = 0; i < array_size; i++)
    elem_array[i] = elementArray + i;
}

// Returns a pointer to an iterator that iterates over the
// set of all vertices in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If vertices are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
Mesquite::VertexIterator* Mesquite::MeshImpl::vertex_iterator(MsqError &/*err*/)
{
  return new Mesquite::MeshImpl_EntityIterator<Mesquite::MeshImpl::Vertex>(vertexArray, vertexCount);
}
    
// Returns a pointer to an iterator that iterates over the
// set of all top-level elements in this mesh.  The calling code should
// delete the returned iterator when it is finished with it.
// If elements are added or removed from the Mesh after obtaining
// an iterator, the behavior of that iterator is undefined.
Mesquite::ElementIterator* Mesquite::MeshImpl::element_iterator(MsqError &/*err*/)
{
  return new Mesquite::MeshImpl_EntityIterator<Mesquite::MeshImpl::Element>(elementArray, elementCount);
}

//************ Vertex Properties ********************
// Returns true or false, indicating whether the vertex
// is allowed to be repositioned.  True indicates that the vertex
// is fixed and cannot be moved.  Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
bool Mesquite::MeshImpl::vertex_is_fixed(Mesquite::Mesh::VertexHandle /*vertex*/,
                                         MsqError &/*err*/)
{
  return false;
}

// Sets on_bnd[] to true or false, indicating whether the vertex
// is on the boundary.  Boundary nodes may be treated as
// a special case by some algorithms or culling methods.
// Note that this is a read-only
// property; this flag can't be modified by users of the
// Mesquite::Mesh interface.
void Mesquite::MeshImpl::vertices_are_on_boundary(
 Mesquite::Mesh::VertexHandle vert_array[], bool on_bnd[],
 size_t num_vtx, MsqError &/*err*/)
{
  for (size_t i=0; i<num_vtx; ++i)
  {
    size_t j = (Vertex*)vert_array[i] - vertexArray;
    if ( (onBoundaryBits[j / 8] & ((unsigned char)(1) << j % 8)) != 0 )
      on_bnd[i] = true;
    else
      on_bnd[i] = false;
  }
}

// Get/set location of a vertex
void Mesquite::MeshImpl::vertices_get_coordinates(
  Mesquite::Mesh::VertexHandle vert_array[],
  Mesquite::MsqVertex* const &coordinates,
  const size_t &num_vtx,
  MsqError &/*err*/)
{
  for (size_t i=0; i<num_vtx; ++i) {
    coordinates[i].set(reinterpret_cast<double*>(vert_array[i]));
  }
}

void Mesquite::MeshImpl::vertex_set_coordinates(
  Mesquite::Mesh::VertexHandle vertex,
  const Vector3D &coordinates,
  MsqError &/*err*/)
{
  coordinates.get_coordinates(reinterpret_cast<Vertex*>(vertex)->coords);
}

// Each vertex has a byte-sized flag that can be used to store
// flags.  This byte's value is neither set nor used by the mesh
// implementation.  It is intended to be used by Mesquite algorithms.
// Until a vertex's byte has been explicitly set, its value is 0.
void Mesquite::MeshImpl::vertex_set_byte (
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char byte,
  MsqError &/*err*/)
{
  size_t index = reinterpret_cast<Vertex*>(vertex) - vertexArray;
  vertexMesquiteByte[index] = byte;
}

void Mesquite::MeshImpl::vertices_set_byte (
  Mesquite::Mesh::VertexHandle *vert_array,
  unsigned char *byte_array,
  size_t array_size,
  MsqError &/*err*/)
{
  for (size_t i = 0; i < array_size; i++)
  {
    size_t index = reinterpret_cast<Vertex*>(vert_array[i]) - vertexArray;
    vertexMesquiteByte[index] = byte_array[i];
  }
}

// Retrieve the byte value for the specified vertex or vertices.
// The byte value is 0 if it has not yet been set via one of the
// *_set_byte() functions.
void Mesquite::MeshImpl::vertex_get_byte(
  Mesquite::Mesh::VertexHandle vertex,
  unsigned char *byte,
  MsqError &/*err*/)
{
  size_t index = reinterpret_cast<Vertex*>(vertex) - vertexArray;
  *byte = vertexMesquiteByte[index];
}

void Mesquite::MeshImpl::vertices_get_byte(
  Mesquite::Mesh::VertexHandle *vertex,
  unsigned char *byte_array,
  size_t array_size,
  MsqError &/*err*/)
{
  for (size_t i = 0; i < array_size; i++)
  {
    size_t index = reinterpret_cast<Vertex*>(vertex[i]) - vertexArray;
    byte_array[i] = vertexMesquiteByte[index];
  }
}

//**************** Vertex Topology *****************

void Mesquite::MeshImpl::create_vertex_to_element_data(MsqError &/*err*/)
{
  if (v2E)
    return;
  
  v2eOffset = new size_t[vertexCount + 1];
  v2E = new size_t[totalVertexUses];
  
    // Initialize each use count to zero
  memset(v2eOffset, 0, (vertexCount+1)*sizeof(size_t));
  
  size_t elem_num, vert_num;
  
    // Go through each element, keep track of how many times
    // each vertex is used.
  for (elem_num = elementCount; elem_num--; )
  {
    for (vert_num =
           Mesquite::vertices_in_topology(elementArray[elem_num].mType);
         vert_num--;
         )
    {
      v2eOffset[elementArray[elem_num].vertexIndices[vert_num]]++;
    }
  }
  
    // Convert the uses counts to array offsets.
  elem_num = 0;
  for (vert_num = 0; vert_num < vertexCount; vert_num++)
  {
    size_t temp = v2eOffset[vert_num];
    v2eOffset[vert_num] = elem_num;
    elem_num += temp;
  }
  v2eOffset[vertexCount] = totalVertexUses;

    // Use newVertIndices to store how many elements
    // have already been added to v2E for each vertex.
    // Should already be initialized to zero
  
    // Finally, store the v2E data
  for (elem_num = 0; elem_num < elementCount; elem_num++)
  {
    for (vert_num =
           Mesquite::vertices_in_topology(elementArray[elem_num].mType);
         vert_num--;
         )
    {
      size_t vert_index = elementArray[elem_num].vertexIndices[vert_num];
      v2E[v2eOffset[vert_index] + newVertIndices[vert_index]++] = elem_num;
    }
  }
  
    // Reset newVertIndices to zero
  memset(newVertIndices, 0, vertexCount*sizeof(size_t));
}

// Gets the number of elements attached to this vertex.
// Useful to determine how large the "elem_array" parameter
// of the vertex_get_attached_elements() function must be.
size_t Mesquite::MeshImpl::vertex_get_attached_element_count(
  Mesquite::Mesh::VertexHandle vertex,
  MsqError &err) const
{
  const_cast<Mesquite::MeshImpl*>(this)->create_vertex_to_element_data(err);
  MSQ_CHKERR(err);
  
  size_t i = reinterpret_cast<MeshImpl::Vertex*>(vertex) - vertexArray;
  return v2eOffset[i+1] - v2eOffset[i];
}

// Gets the elements attached to this vertex.
void Mesquite::MeshImpl::vertex_get_attached_elements(
  Mesquite::Mesh::VertexHandle vertex,
  Mesquite::Mesh::ElementHandle* elem_array,
  size_t sizeof_elem_array,
  MsqError &err)
{
  create_vertex_to_element_data(err); MSQ_CHKERR(err);
  
  size_t index = reinterpret_cast<MeshImpl::Vertex*>(vertex) - vertexArray;
  
  if (sizeof_elem_array > v2eOffset[index+1] - v2eOffset[index])
    sizeof_elem_array = v2eOffset[index+1] - v2eOffset[index];
  
  for ( ; sizeof_elem_array--; )
  {
    elem_array[sizeof_elem_array] =
      elementArray + v2E[v2eOffset[index] + sizeof_elem_array];
  }
}


// Gets the number of vertices in this element.
// This data can also be found by querying the
// element's topology and getting the number
// of vertices per element for that topology type.
size_t Mesquite::MeshImpl::element_get_attached_vertex_count(
  Mesquite::Mesh::ElementHandle elem,
  MsqError &/*err*/) const
{
  return Mesquite::vertices_in_topology(reinterpret_cast<MeshImpl::Element*>(elem)->mType);
}

// Returns the vertices that are part of the topological definition of each
// element in the "elem_handles" array.  When this function is called, the
// following must be true:
//   a) "elem_handles" points at an array of "num_elems" element handles.
//   b) "vert_handles" points at an array of size "sizeof_vert_handles"
//   c) "csr_data" points at an array of size "sizeof_csr_data"
//   d) "csr_offsets" points at an array of size "num_elems+1"
//      
// When this function returns, adjacency information will be stored
// in csr format:
//    a) "vert_handles" stores handles to all vertices found in one
//       or more of the elements.  Each vertex appears only
//       once in "vert_handles", even if it is in multiple elements.
//    b) "sizeof_vert_handles" is set to the number of vertex
//       handles placed into "vert_handles".
//    c) "sizeof_csr_data" is set to the total number of vertex uses (for
//       example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
//       the two triangles share some vertices).
//    c) "csr_offsets" is filled such that csr_offset[i] indicates the location
//       of entity i's first adjacency in "csr_data".  The number of vertices
//       in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
//       reason, csr_offsets[num_elems] is set to the new value of
//       "sizeof_csr_data".
//    d) "csr_data" stores integer offsets which give the location of
//       each adjacency in the "vert_handles" array.
//
// As an example of how to use this data, you can get the handle of the first
// vertex in element #3 like this:
//   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ]
//
// and the second vertex of element #3 like this:
//   VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ]
// 
void Mesquite::MeshImpl::elements_get_attached_vertices(
  Mesquite::Mesh::ElementHandle *elem_handles,
  size_t num_elems,
  Mesquite::Mesh::VertexHandle *vert_handles,
  size_t &sizeof_vert_handles,
  size_t *csr_data,
  size_t &sizeof_csr_data,
  size_t *csr_offsets,
  MsqError &err)
{
  if (num_elems == 0)
    return;
  
  size_t total_verts = 0;
  csr_offsets[0] = 0;
  
    // for each element
  for (size_t i = 0; i < num_elems; i++)
  {
    MeshImpl::Element* elem =
      reinterpret_cast<MeshImpl::Element*>(elem_handles[i]);
    size_t verts_in_elem = Mesquite::vertices_in_topology(elem->mType);
    csr_offsets[i+1] = csr_offsets[i] + verts_in_elem;

      // Make sure we've got enough room in csr_data
    if (sizeof_csr_data < csr_offsets[i+1]) {
        // Error!!!
      err.set_msg("Arg. sizeof_csr_data is too small.");
      return;
    }
    
      // for each vertex in this element
    for (size_t j = 0; j < verts_in_elem; j++)
    {
        // Get the index for this vertex
      size_t vert_index = elem->vertexIndices[j];
      size_t new_vert_index = newVertIndices[vert_index];
      if (new_vert_index == 0)
      {
        new_vert_index = newVertIndices[vert_index] = ++total_verts;
        if (total_verts > sizeof_vert_handles) // ERROR!!!
          return;
        vert_handles[total_verts-1] =
          reinterpret_cast<Mesquite::Mesh::VertexHandle>(vertexArray+vert_index);
      }
      
        // Place that vertex into the csr_data array
      csr_data[csr_offsets[i] + j] = new_vert_index-1;
    }
  }
  
    // Set the amount of data we are returning
  sizeof_csr_data = csr_offsets[num_elems];
  sizeof_vert_handles = total_verts;
  
    // Set newVertIndices back to 0
  for ( ; total_verts--; )
  {
    newVertIndices[reinterpret_cast<Mesquite::MeshImpl::Vertex*>(vert_handles[total_verts]) - vertexArray] = 0;
  }
}

// Identifies the vertices attached to the elements by returning
// each vertex's global index.  The vertex's global index indicates
// where that element can be found in the array returned by
// Mesh::get_all_vertices. The indexes can be repeated.
// \param elems Array of element handles.
// \param num_elems number of elements in the array elems
// \param index_array Array containing the indexes of the elements vertices.
//                    Indexes can be repeated.
// \param index_array_size indicates the size of index_array
// \param offsets An array of offsets into the index_array that indicates where
//                the indexes corresponding to an element start. First entry is 0.
// For example element, to get the vertex indices of elems[3], we look at the entries
// index_array[offsets[3]] to index_array[offsets[4]] . 
#undef __FUNC__
#define __FUNC__ "MeshImpl::elements_get_attached_vertex_indices"
void Mesquite::MeshImpl::elements_get_attached_vertex_indices(
  Mesquite::Mesh::ElementHandle elems[],
  size_t num_elems,
  size_t index_array[],
  size_t array_size,
  size_t* offsets,
  MsqError &err)
{
  offsets[0] = 0;
  for (size_t i=0; i<num_elems; ++i) {
    Mesquite::MeshImpl::Element* elem =
      reinterpret_cast<Mesquite::MeshImpl::Element*>(elems[i]);

    size_t num_verts = Mesquite::vertices_in_topology(elem->mType);
    offsets[i+1] = offsets[i] + num_verts;
    
    for ( ; num_verts--; )
    {
      index_array[offsets[i]+num_verts] = elem->vertexIndices[num_verts];
    }
  }
}

// Returns the topology of the given entity.
Mesquite::EntityTopology Mesquite::MeshImpl::element_get_topology(
  Mesquite::Mesh::ElementHandle entity_handle,
  MsqError &/*err*/) const 
{
  return reinterpret_cast<MeshImpl::Element*>(entity_handle)->mType;
}

// Returns the topologies of the given entities.  The "entity_topologies"
// array must be at least "num_elements" in size.
void Mesquite::MeshImpl::elements_get_topologies(
  Mesquite::Mesh::ElementHandle element_handle_array[],
  Mesquite::EntityTopology element_topologies[],
  size_t num_elements,
  MsqError &/*err*/)
{
  for (size_t i = 0; i < num_elements; i++)
  {
    element_topologies[i] =
      reinterpret_cast<MeshImpl::Element*>(element_handle_array[i])->mType;
  }
}


//*************** Dense Tags (i.e. tags set on all entities) ***********

#undef __FUNC__
#define __FUNC__ "MeshImpl::element_tag_create"
//! only dense tags are implemented in MeshImpl for now.
void Mesquite::MeshImpl::element_tag_create(const string tag_name, int tag_size,
                                        TagHandle& tag_handle,
                                        MsqError &err)
{
  tag new_tag;
  new_tag.pt = 0;
  new_tag.size = tag_size;
  denseTags[tag_name] = new_tag;
  tag_handle = tag_get_handle(tag_name, err); MSQ_CHKERR(err);
}
    
#undef __FUNC__
#define __FUNC__ "MeshImpl::tag_get_handle"
//! only dense tags are implemented in MeshImpl for now.
void* Mesquite::MeshImpl::tag_get_handle(const string tag_name, MsqError &err)
{
  return (void*) &(denseTags[tag_name]);
}
    
#undef __FUNC__
#define __FUNC__ "MeshImpl::elements_set_tag_data"
//! only dense tags are implemented in MeshImpl for now.
void Mesquite::MeshImpl::elements_set_tag_data(
                                   const size_t num_elements, 
                                   const TagHandle tag_handle,
                                   TagDataPt const tag_data_array,
                                   const int& tag_size,
                                   MsqError &err)
{
  if (num_elements != elementCount) {
    err.set_msg("Incorrect num_elements. Must be equal to the total number of elements "
                "since only dense tags are supported.");
    return;
  }
  else if (((tag*)tag_handle)->size != tag_size) {
    err.set_msg("tag_size does not correspong to actual tag size.");
    return;
  }
  else
   ((tag*)tag_handle)->pt = tag_data_array;
}


#undef __FUNC__
#define __FUNC__ "MeshImpl::elements_get_tag_data"
//! only dense tags are implemented in MeshImpl for now.
void Mesquite::MeshImpl::elements_get_tag_data(
                                       const size_t num_elements,
                                       const TagHandle tag_handle,
                                       TagDataPt &tag_data_array,
                                       int& tag_size,
                                       MsqError &err)
{
  if (num_elements != elementCount) {
    err.set_msg("Incorrect num_elements. Must be equal to the total number of elements "
                "since only dense tags are supported.");
    return;
  }
  else if (((tag*)tag_handle)->size != tag_size) {
    err.set_msg("tag_size does not correspong to actual tag size.");
    return;
  }
  else
   tag_data_array = ((tag*)tag_handle)->pt;
  
}
    



//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
void Mesquite::MeshImpl::release_entity_handles(
  Mesquite::Mesh::EntityHandle */*handle_array*/,
  size_t /*num_handles*/,
  MsqError &/*err*/)
{
    // Do nothing
}

// Instead of deleting a Mesh when you think you are done,
// call release().  In simple cases, the implementation could
// just call the destructor.  More sophisticated implementations
// may want to keep the Mesh object to live longer than Mesquite
// is using it.
void Mesquite::MeshImpl::release()
{
  delete this;
}
