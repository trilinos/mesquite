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
#include "FileTokenizer.hpp"
#include "Vector3D.hpp"
#include "MsqVertex.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <string.h>
#  include <vector.h>
#else
#  include <string>
#  include <vector>
   using std::string;
   using std::vector;
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <fstream.h>
#  include <iomanip.h>
#else
#  include <fstream>
#  include <iomanip>
   using std::ifstream;
   using std::ofstream;
   using std::endl;
#endif

#ifdef MSQ_USING_EXODUS
#include "exodusII.h"
#endif

#include "MsqDebug.hpp"
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
  clear();
}

void Mesquite::MeshImpl::clear()
{
  delete [] vertexArray;
  delete [] elementArray;
  delete [] onBoundaryBits;
  delete [] vertexMesquiteByte;
  delete [] newVertIndices;
  delete [] v2eOffset;
  delete [] v2E;
    
  vertexArray = 0;
  elementArray = 0;
  onBoundaryBits = 0;
  vertexMesquiteByte = 0;
  newVertIndices = 0;
  v2eOffset = 0;
  v2E = 0;
  
  vertexCount = 0;
  elementCount = 0;
  numCoords = 0;
  totalVertexUses = 0;
  
  for (unsigned i = 0; i < tagList.size(); ++i)
    delete tagList[i];
  tagList.clear();
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
    MSQ_SETERR(err)( MsqError::FILE_ACCESS );
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
    file << vertexArray[i].coords[0] << ' '
         << vertexArray[i].coords[1] << ' '
         << vertexArray[i].coords[2] << '\n';
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
      MSQ_SETERR(err)( "element type not implemented",MsqError::NOT_IMPLEMENTED );
      break;
    }
    file << (int)type_id << '\n';
  }
  
    // Write out which points are fixed.
  file << "POINT_DATA " << vertexCount
       << "\nSCALARS fixed bit\nLOOKUP_TABLE default\n";
  for (i = 0; i < vertexCount; i++)
  {
    if (onBoundaryBits[i/8] & ((unsigned char)(1) << i % 8))
      file << "1\n";
    else
      file << "0\n";
  }
  
    // Write vertex tag data to vtk attributes
  for (i = 0; i < tagList.size(); ++i)
    if (tagList[i] && tagList[i]->vertexData)
    { 
      vtk_write_attrib_data( file, 
                             tagList[i]->desc, 
                             tagList[i]->vertexData,
                             vertexCount,
                             err );
      MSQ_ERRRTN(err);
    }
  
    // If there are any element attributes, write them
  for (i = 0; i< tagList.size(); ++i)
    if (tagList[i] && tagList[i]->elementData)
    {
      file << "\nCELL_DATA " << elementCount << "\n";
      break;
    }
  for (i = 0; i < tagList.size(); ++i)
    if (tagList[i] && tagList[i]->elementData)
    { 
      vtk_write_attrib_data( file, 
                             tagList[i]->desc, 
                             tagList[i]->elementData,
                             elementCount,
                             err );
      MSQ_ERRRTN(err);
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
  MSQ_SETERR(err)( MsqError::NOT_IMPLEMENTED );
  return;
#else
  if (vertexArray != NULL)
  {
    MSQ_SETERR(err)( MsqError::INVALID_STATE );
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
    MSQ_SETERR(err)( MsqError::FILE_ACCESS );
    return;
  }
  
    // make sure the file is saved as doubles
  if (file_float_size != sizeof(double))
  {
    MSQ_SETERR(err)("File saved with float-sized reals.  Can only read files "
                    "saved with doubles.", MsqError::NOT_IMPLEMENTED );
    return;
  }

  char title[MAX_LINE_LENGTH];
  int dim, vert_count, elem_count, block_count, ns_count, ss_count;
  
    // get info about the file
  exo_err = ex_get_init(file_id, title, &dim, &vert_count,
                        &elem_count, &block_count, &ns_count, &ss_count);
  if (exo_err < 0)
  {
    MSQ_SETERR(err)("Unable to get entity counts from file.", 
                    MsqError::PARSE_ERROR);
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
    MSQ_SETERR(err)("Unable to retrieve vertex coordinates from file.",
                    MsqError::PARSE_ERROR);
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
    MSQ_SETERR(err)("Unable to read block IDs from file.", MsqError::PARSE_ERROR);
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
      MSQ_SETERR(err)("Unable to read parameters for block.",MsqError::PARSE_ERROR);
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
      MSQ_SETERR(err)("Unrecognized element type in block",
                      MsqError::PARSE_ERROR);
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
      MSQ_SETERR(err)("Unable to read element block connectivity.",
                      MsqError::PARSE_ERROR);
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
      MSQ_PRINT(1)("\nError opening nodeset 111, no boundary nodes marked.");
      num_fixed_nodes=0;
    }
  }
  int *fixed_nodes =NULL;
  if(num_fixed_nodes>0)
    fixed_nodes = new int[num_fixed_nodes];
  exo_err = ex_get_node_set(file_id, 111, fixed_nodes);
  if(exo_err<0){
    MSQ_SETERR(err)("Error retrieving fixed nodes.", MsqError::PARSE_ERROR);
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
    MSQ_SETERR(err)("Error closing Exodus file.", MsqError::IO_ERROR);
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
  MSQ_SETERR(err)("Exodus not enabled in this build of Mesquite",
                  MsqError::NOT_IMPLEMENTED);
  return;
#else
  size_t i, j, counter;
  if (vertexArray == NULL)
  {
    MSQ_SETERR(err)("No vertices in MeshImpl.  Nothing written to file.", 
                    MsqError::PARSE_ERROR);
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
    else {
      MSQ_SETERR(err)("Unrecognized element type", MsqError::NOT_IMPLEMENTED);
      return;
    }
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
    MSQ_SETERR(err)("Too many element types in file",MsqError::PARSE_ERROR);
    return;
  }
  if(block_count<1){
    MSQ_SETERR(err)("No elements in file.",MsqError::PARSE_ERROR);
    return;
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
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }
  
  char title[MAX_LINE_LENGTH]="Mesquite Generated Exodus File";
  int dim=3;
  
  size_t vert_count=0;
  size_t elem_count=0;
  size_t temp_var=0;
  get_all_sizes(vert_count, elem_count, temp_var, err);
  
  int ns_count=0;
  if(num_fixed_nodes>0)
    ns_count=1;
  int ss_count=0;
  
    // put the initial info about the file
  exo_err = ex_put_init(file_id, title, dim, vert_count,
                        elem_count, block_count, ns_count, ss_count);
  if (exo_err < 0)
  {
    MSQ_SETERR(err)("Unable to initialize file data.", MsqError::IO_ERROR);
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
  {
    MSQ_SETERR(err)("Counter at incorrect number.", MsqError::INTERNAL_ERROR);
    return;
  }
    //put the coords
  exo_err = ex_put_coord(file_id,
                         reinterpret_cast<void*>(temp_doubles),
                         reinterpret_cast<void*>(temp_doubles + vertexCount),
                         reinterpret_cast<void*>(temp_doubles +
                                                 vertexCount + vertexCount));
  
    // Make sure it worked
  if (exo_err < 0)
  {
    MSQ_SETERR(err)("Unable to put vertex coordinates in file.",MsqError::IO_ERROR);
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
    MSQ_SETERR(err)("Error in determining the element types.",MsqError::INTERNAL_ERROR);
    return;
  }
  int num_atts=0;

    //for each element type that exists, set up the block
  if(num_tri>0)
    exo_err = ex_put_elem_block(file_id, block_ids[0], "TRI3",
                                num_tri, 3, num_atts);
  if(exo_err<0)
  {
    MSQ_SETERR(err)("Error creating the tri block.", MsqError::IO_ERROR);
    return;
  }
  if(num_quad>0)
    exo_err = ex_put_elem_block(file_id, block_ids[1], "SHELL",
                                num_quad, 4, num_atts);
  if(exo_err<0)
  {
    MSQ_SETERR(err)("Error creating the quad block.",MsqError::IO_ERROR);
    return;
  }
  if(num_tet>0)
    exo_err = ex_put_elem_block(file_id, block_ids[2], "TETRA",
                                num_tet, 4, num_atts);
  if(exo_err<0)
  {
    MSQ_SETERR(err)("Error creating the tet block.",MsqError::IO_ERROR);
    return;
  }
  if(num_hex>0)
    exo_err = ex_put_elem_block(file_id, block_ids[3], "HEX",
                                num_hex, 8, num_atts);
  if(exo_err<0)
  {
    MSQ_SETERR(err)("Error creating the hex block.",MsqError::IO_ERROR);
    return;
  }
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
    {
      MSQ_SETERR(err)("Unrecognized element type again", MsqError::NOT_IMPLEMENTED);
      return;
    }
  }
    //double check number of each element type
  if(tri_counter != num_tri)
  {
    MSQ_SETERR(err)("Tri numbers not consistent", MsqError::INTERNAL_ERROR);
    return;
  }
  if(quad_counter != num_quad)
  {
    MSQ_SETERR(err)("Quad numbers not consistent", MsqError::INTERNAL_ERROR);
    return;
  }
  if(tet_counter != num_tet)
  {
    MSQ_SETERR(err)("Tet numbers not consistent", MsqError::INTERNAL_ERROR);
    return;
  }
  if(hex_counter != num_hex)
  {
    MSQ_SETERR(err)("Hex numbers not consistent", MsqError::INTERNAL_ERROR);
    return;
  }
  
  //for each element type that exists, put the connectivity with the block
  if(num_tri>0)
    exo_err = ex_put_elem_conn(file_id, block_ids[0], tri_connectivity);
  if(exo_err<0) {
    MSQ_SETERR(err)("Error in putting tri block", MsqError::IO_ERROR);
    return;
  }
  if(num_quad>0)
    exo_err = ex_put_elem_conn(file_id, block_ids[1], quad_connectivity);
  if(exo_err<0) {
    MSQ_SETERR(err)("Error in putting quad block", MsqError::IO_ERROR);
    return;
  }
  if(num_tet>0)
    exo_err = ex_put_elem_conn(file_id, block_ids[2], tet_connectivity);
  if(exo_err<0) {
    MSQ_SETERR(err)("Error in putting tet block", MsqError::IO_ERROR);
    return;
  }
  if(num_hex>0)
    exo_err = ex_put_elem_conn(file_id, block_ids[3], hex_connectivity);
  if(exo_err<0) {
    MSQ_SETERR(err)("Error in putting hex block", MsqError::IO_ERROR);
    return;
  }
    //delete connectivity arrays.
  delete [] tri_connectivity;
  delete [] quad_connectivity;
  delete [] tet_connectivity;
  delete [] hex_connectivity;
  
    // Finally, mark boundary nodes
  
  if(num_fixed_nodes>0){
    exo_err=ex_put_node_set_param(file_id, 111, num_fixed_nodes, 0);
    if(exo_err<0) {
      MSQ_SETERR(err)("Error while initializing node set.", MsqError::IO_ERROR);
      return;
    }
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
    if(exo_err<0) {
      MSQ_SETERR(err)("Error while writing node set.", MsqError::IO_ERROR);
      return;
    }
  }
  exo_err=ex_close(file_id);
  if(exo_err<0)
    MSQ_SETERR(err)("Error closing Exodus file.", MsqError::IO_ERROR);
  
#endif
}   
// Returns whether this mesh lies in a 2D or 3D coordinate system.
int Mesquite::MeshImpl::get_geometric_dimension(MsqError &/*err*/)
{
  return numCoords;
}


void Mesquite::MeshImpl::get_all_sizes( size_t& vertex_count,
                                        size_t& element_count,
                                        size_t& vertex_use_count,
                                        MsqError& )
{
  vertex_count = vertexCount;
  element_count = elementCount;
  vertex_use_count = totalVertexUses;
}

void Mesquite::MeshImpl::get_all_mesh( 
                               VertexHandle*  vert_array, size_t vert_len,
                               ElementHandle* elem_array, size_t elem_len,
                               size_t* elem_conn_offsets, size_t offset_len,
                               size_t* elem_conn_indices, size_t index_len,
                               MsqError& err )
{
  if (vert_len < vertexCount ||
      elem_len < elementCount ||
      offset_len < elementCount+1 ||
      index_len < totalVertexUses) {
    MSQ_SETERR(err)("Insufficient space in array", MsqError::INVALID_ARG );
    return;
  }
  
  for( size_t i = 0; i < vertexCount; ++i)
    vert_array[i] = vertexArray + i;
  
  size_t index = 0;
  ElementHandle* handle_iter = elem_array;
  size_t* conn_iter = elem_conn_indices;
  size_t* offset_iter = elem_conn_offsets;
  size_t* const conn_end = conn_iter + index_len;
  Element* const end = elementArray + elementCount;
  for(Element* iter = elementArray; iter != end; ++iter)
  {
    *handle_iter = iter; ++handle_iter;
    *offset_iter = index; ++offset_iter;
    size_t len = Mesquite::vertices_in_topology(iter->mType);
    index += len;
    if (conn_iter + len > conn_end) {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
      return;
    }
    memcpy( conn_iter, iter->vertexIndices, len*sizeof(size_t));
    conn_iter += len;
  }
  *offset_iter = index;
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
  MsqError &err) 
{
  const_cast<Mesquite::MeshImpl*>(this)->create_vertex_to_element_data(err);
  MSQ_ERRZERO(err);
  
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
  create_vertex_to_element_data(err); MSQ_ERRRTN(err);
  
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
  MsqError &/*err*/)
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
      MSQ_SETERR(err)("Arg. sizeof_csr_data is too small.",MsqError::INVALID_ARG);
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
        if (total_verts > sizeof_vert_handles) {// ERROR!!!
          MSQ_SETERR(err)("Insufficient space for vertex handles.", MsqError::INVALID_ARG);
          return;
        }
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


// Returns the topology of the given entity.
Mesquite::EntityTopology Mesquite::MeshImpl::element_get_topology(
  Mesquite::Mesh::ElementHandle entity_handle,
  MsqError &/*err*/) 
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





//**************** Memory Management ****************
// Tells the mesh that the client is finished with a given
// entity handle.  
void Mesquite::MeshImpl::release_entity_handles(
  Mesquite::Mesh::EntityHandle* /*handle_array*/,
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



namespace Mesquite {


const char* const vtk_type_names[] = { "bit",
                                       "char",
                                       "unsigned_char",
                                       "short",
                                       "unsigned_short",
                                       "int",
                                       "unsigned_int",
                                       "long",
                                       "unsigned_long",
                                       "float",
                                       "double",
                                       0 };

void MeshImpl::read_vtk( const char* filename, Mesquite::MsqError &err )
{
  int major, minor;
  char vendor_string[257];
  size_t i;
  
  FILE* file = fopen( filename, "r" );
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }
  
    // Read file header
    
  if (!fgets( vendor_string, sizeof(vendor_string), file ))
  {
    MSQ_SETERR(err)( MsqError::IO_ERROR );
    fclose( file );
    return;
  }
  
  if (!strchr( vendor_string, '\n' ) ||
      2 != sscanf( vendor_string, "# vtk DataFile Version %d.%d", &major, &minor ))
  {
    MSQ_SETERR(err)( MsqError::FILE_FORMAT );
    fclose( file );
    return;
  }
  
  if (!fgets( vendor_string, sizeof(vendor_string), file )) 
  {
    MSQ_SETERR(err)( MsqError::IO_ERROR );
    fclose( file );
    return;
  }
  
    // VTK spec says this should not exceed 256 chars.
  if (!strchr( vendor_string, '\n' ))
  {
    MSQ_SETERR(err)( "Vendor string (line 2) exceeds 256 characters.",
                      MsqError::PARSE_ERROR);
    fclose( file );
    return;
  }
  
  
    // Check file type
  
  FileTokenizer tokens( file );
  const char* const file_type_names[] = { "ASCII", "BINARY", 0 };
  int filetype = tokens.match_token( file_type_names, err ); MSQ_ERRRTN(err);
  if (2 == filetype) {
    MSQ_SETERR(err)( "Cannot read BINARY VTK files -- use ASCII.",
                     MsqError::NOT_IMPLEMENTED );
    return;
  }

    // Clear any existing data
  this->clear();

    // Read the mesh
  tokens.match_token( "DATASET", err ); MSQ_ERRRTN(err);
  vtk_read_dataset( tokens, err ); MSQ_ERRRTN(err);
  
    // Make sure file actually contained some mesh
  if (elementCount == 0)
  {
    MSQ_SETERR(err)("File contained no mesh.", MsqError::PARSE_ERROR);
    return;
  }
  
    // Read attribute data until end of file.
  const char* const block_type_names[] = { "POINT_DATA", "CELL_DATA", 0 };
  int blocktype = 0;
  while (!tokens.eof())
  {
      // get POINT_DATA or CELL_DATA
    int new_block_type = tokens.match_token( block_type_names, err );
    if (tokens.eof())
    {
      err.clear();
      break;
    }
    if (err)
    {
        // If next token was neither POINT_DATA nor CELL_DATA,
        // then there's another attribute under the current one.
      if (blocktype)
      {
        tokens.unget_token();
        err.clear();
      }
      else
      {
        MSQ_ERRRTN(err);
      }
    }
    else
    {
      blocktype = new_block_type;
      long count;
      tokens.get_long_ints( 1, &count, err); MSQ_ERRRTN(err);
      
      if (blocktype == 1 && (unsigned long)count != vertexCount)
      {
        MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                         "Count inconsistent with number of vertices" 
                         " at line %d.", tokens.line_number());
        return;
      }
      else if (blocktype == 2 && (unsigned long)count != elementCount)
      {
         MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                         "Count inconsistent with number of elements" 
                         " at line %d.", tokens.line_number());
        return;
      }
    }
      
   
    if (blocktype == 1)
      vtk_read_point_data( tokens, err );
    else
      vtk_read_cell_data ( tokens, err );
    MSQ_ERRRTN(err);
  }
  
    // It seems that this variable needs to be initialized here
  totalVertexUses = 0;
  for (i = 0; i < elementCount; i++)
    switch( elementArray[i].mType ) {
      case TRIANGLE     : totalVertexUses += 3; break;
      case QUADRILATERAL: totalVertexUses += 4; break;
      case TETRAHEDRON  : totalVertexUses += 4; break;
      case HEXAHEDRON   : totalVertexUses += 8; break;
      case PRISM        : totalVertexUses += 6; break;
      case PYRAMID      : totalVertexUses += 5; break;
      case SEPTAHEDRON  : totalVertexUses += 7; break;
      default:
        MSQ_SETERR(err)( "Unimplemented element type encountered during vertex-use calculation", 
                         MsqError::INTERNAL_ERROR );
        return;
    }
  
    // Apparently MeshImpl needs all these arrays created.
  size_t num_bitset_bytes = (vertexCount / 8);
  num_bitset_bytes += (vertexCount % 8) ? 1 : 0;
  onBoundaryBits = new unsigned char[num_bitset_bytes];
  memset( onBoundaryBits, 0, num_bitset_bytes );
  vertexMesquiteByte = new unsigned char[vertexCount];
  memset( vertexMesquiteByte, 0, vertexCount );
  newVertIndices = new size_t[vertexCount];
  memset( newVertIndices, 0, vertexCount * sizeof(size_t) );
  numCoords = 3;

    // Convert tag data for fixed nodes to internal bitmap
  MsqError tmperr;
  TagHandle handle = tag_get( "fixed", tmperr );
  if (tmperr) return;
  
  TagData* tag = tag_from_handle( handle, err ); MSQ_ERRRTN(err);
  if (tag->desc.size / size_from_tag_type(tag->desc.type) != 1)
  {
    MSQ_SETERR(err)("'fixed' attribute is not a scalar value", MsqError::FILE_FORMAT);
    return;
  } 
  
  if (!tag->vertexData || tag->vertIsDefault)
  {
    MSQ_SETERR(err)("'fixed' attribute on elements, not vertices", MsqError::FILE_FORMAT);
    return;
  }

  void* tagdata = tag->vertexData;
  switch( tag->desc.type )
  {
    case BYTE:
      for (i = 0; i < vertexCount; ++i)
        if (((char*)tagdata)[i])
          onBoundaryBits[i / 8] |= ((unsigned char)1 << (i%8));
      break;
    case BOOL:  
      for (i = 0; i < vertexCount; ++i)
        if (((bool*)tagdata)[i])
          onBoundaryBits[i / 8] |= ((unsigned char)1 << (i%8));
      break;
    case INT :  
      for (i = 0; i < vertexCount; ++i)
        if (((int*)tagdata)[i])
          onBoundaryBits[i / 8] |= ((unsigned char)1 << (i%8));
      break;
    case DOUBLE:
      for (i = 0; i < vertexCount; ++i)
        if (((double*)tagdata)[i])
          onBoundaryBits[i / 8] |= ((unsigned char)1 << (i%8));
      break;
    default:
      MSQ_SETERR(err)("'fixed' attribute has invalid type", MsqError::PARSE_ERROR);
      return;
  }
  
  tag_delete( handle, err );
}

void MeshImpl::vtk_read_dataset( FileTokenizer& tokens, MsqError& err )
{
  const char* const data_type_names[] = { "STRUCTURED_POINTS",
                                          "STRUCTURED_GRID",
                                          "UNSTRUCTURED_GRID",
                                          "POLYDATA",
                                          "RECTILINEAR_GRID",
                                          "FIELD",
                                          0 };
  int datatype = tokens.match_token( data_type_names, err ); MSQ_ERRRTN(err);
  switch( datatype )
  {
    case 1: vtk_read_structured_points( tokens, err ); break;
    case 2: vtk_read_structured_grid  ( tokens, err ); break;
    case 3: vtk_read_unstructured_grid( tokens, err ); break;
    case 4: vtk_read_polydata         ( tokens, err ); break;
    case 5: vtk_read_rectilinear_grid ( tokens, err ); break;
    case 6: vtk_read_field            ( tokens, err ); break;
  }
}


void MeshImpl::vtk_read_structured_points( FileTokenizer& tokens, MsqError& err )
{
  long i, j, k;
  long dims[3];
  double origin[3], space[3];
 
  tokens.match_token( "DIMENSIONS", err ); MSQ_ERRRTN(err);
  tokens.get_long_ints( 3, dims, err );    MSQ_ERRRTN(err);
  tokens.get_newline( err );               MSQ_ERRRTN(err);
  
  if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1)
  {
    MSQ_SETERR(err)(MsqError::PARSE_ERROR,
                   "Invalid dimension at line %d", 
                   tokens.line_number());
    return;
  }
  
  tokens.match_token( "ORIGIN", err );     MSQ_ERRRTN(err);
  tokens.get_doubles( 3, origin, err );    MSQ_ERRRTN(err);
  tokens.get_newline( err );               MSQ_ERRRTN(err);
  
  const char* const spacing_names[] = { "SPACING", "ASPECT_RATIO", 0 };
  tokens.match_token( spacing_names, err );MSQ_ERRRTN(err);
  tokens.get_doubles( 3, space, err );     MSQ_ERRRTN(err);
  tokens.get_newline( err );               MSQ_ERRRTN(err);
  
  vertexCount = dims[0] * dims[1] * dims[2];
  Vertex* vtx_ptr = vertexArray = new Vertex[vertexCount];

  for (k = 0; k < dims[2]; ++k)
    for (j = 0; j < dims[1]; ++j)
      for (i = 0; i < dims[0]; ++i)
      {
        vtx_ptr->coords[0] = origin[0] + i * space[0];
        vtx_ptr->coords[1] = origin[1] + j * space[1];
        vtx_ptr->coords[2] = origin[2] + k * space[2];
        ++vtx_ptr;
      }
  
  vtk_create_structured_elems( dims, err ); MSQ_ERRRTN(err);
}

void MeshImpl::vtk_read_structured_grid( FileTokenizer& tokens, MsqError& err )
{
  long num_verts, dims[3];
 
  tokens.match_token( "DIMENSIONS", err );    MSQ_ERRRTN(err);
  tokens.get_long_ints( 3, dims, err );      MSQ_ERRRTN(err);
  tokens.get_newline( err );                 MSQ_ERRRTN(err); 
  
  if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1)
  {
    MSQ_SETERR(err)(MsqError::PARSE_ERROR,
                   "Invalid dimension at line %d", 
                   tokens.line_number());
    return;
  }
  
  tokens.match_token( "POINTS", err );        MSQ_ERRRTN(err); 
  tokens.get_long_ints( 1, &num_verts, err ); MSQ_ERRRTN(err); 
  tokens.match_token( vtk_type_names, err );  MSQ_ERRRTN(err); 
  tokens.get_newline( err );                  MSQ_ERRRTN(err); 
  
  if (num_verts != (dims[0] * dims[1] * dims[2]))
  {
    MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                    "Point count not consistent with dimensions "
                    "at line %d", tokens.line_number() );
    return;
  }
  
  vertexCount = dims[0] * dims[1] * dims[2];
  vertexArray = new Vertex[num_verts];
  
  tokens.get_doubles( 3 * vertexCount, (double*)vertexArray, err ); MSQ_ERRRTN(err);
  
  vtk_create_structured_elems( dims, err ); MSQ_ERRRTN(err);
}

void MeshImpl::vtk_read_rectilinear_grid( FileTokenizer& tokens, MsqError& err )
{
  int i, j, k;
  long dims[3];
  const char* labels[] = { "X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES" };
  vector<double> coordinates[3];
  
  tokens.match_token( "DIMENSIONS", err );                   MSQ_ERRRTN(err);
  tokens.get_long_ints( 3, dims, err );                     MSQ_ERRRTN(err);
  tokens.get_newline( err );                                MSQ_ERRRTN(err);
  
  if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1)
  {
     MSQ_SETERR(err)(MsqError::PARSE_ERROR,
                   "Invalid dimension at line %d", 
                   tokens.line_number());
    return;
  }

  for (i = 0; i < 3; i++)
  {
    long count;
    tokens.match_token( labels[i], err );                   MSQ_ERRRTN(err);
    tokens.get_long_ints( 1, &count, err );                 MSQ_ERRRTN(err);
    tokens.match_token( vtk_type_names, err );              MSQ_ERRRTN(err);
    tokens.get_newline( err );                              MSQ_ERRRTN(err);
    
    if (count != dims[i])
    {
      MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                       "Coordinate count inconsistent with dimensions"
                       " at line %d", tokens.line_number());
      return;
    }
    
    coordinates[i].resize(count);
    tokens.get_doubles( count, &coordinates[i][0], err );   MSQ_ERRRTN(err);
  }
  
  vertexCount = dims[0] * dims[1] * dims[2];
  Vertex* vtx_ptr = vertexArray = new Vertex[vertexCount];

  for (k = 0; k < dims[2]; ++k)
    for (j = 0; j < dims[1]; ++j)
      for (i = 0; i < dims[0]; ++i)
      {
        vtx_ptr->coords[0] = coordinates[0][i];
        vtx_ptr->coords[1] = coordinates[1][j];
        vtx_ptr->coords[2] = coordinates[2][k];
        ++vtx_ptr;
      }
  
  
  vtk_create_structured_elems( dims, err );                 MSQ_ERRRTN(err);
}

void MeshImpl::vtk_read_polydata( FileTokenizer& tokens, MsqError& err )
{
  long num_verts;
  vector<int> connectivity;
  const char* const poly_data_names[] = { "VERTICES",
                                          "LINES",
                                          "POLYGONS",
                                          "TRIANGLE_STRIPS", 
                                          0 };
  
  tokens.match_token( "POINTS", err );                      MSQ_ERRRTN(err);
  tokens.get_long_ints( 1, &num_verts, err );               MSQ_ERRRTN(err);
  tokens.match_token( vtk_type_names, err );                MSQ_ERRRTN(err);
  tokens.get_newline( err );                                MSQ_ERRRTN(err);
  
  if (num_verts < 1)
  {
    MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                     "Invalid point count at line %d", tokens.line_number());
    return;
  }
  
  Vertex* vtx_ptr = vertexArray = new Vertex[num_verts];
  tokens.get_doubles( 3*num_verts, (double*)vtx_ptr, err ); MSQ_ERRRTN(err);
  vertexCount = num_verts;

  int poly_type = tokens.match_token( poly_data_names, err );MSQ_ERRRTN(err);
  switch (poly_type)
  {
    case 3:
      vtk_read_polygons( tokens, err );                     MSQ_ERRRTN(err);
      break;
    case 4:
      MSQ_SETERR(err)( MsqError::NOT_IMPLEMENTED, 
                       "Unsupported type: triangle strips at line %d",
                       tokens.line_number() );
      return;
    case 1:
    case 2:
      MSQ_SETERR(err)( MsqError::NOT_IMPLEMENTED, 
                       "Entities of dimension < 2 at line %d",
                       tokens.line_number() );
      return;
  }
}

void MeshImpl::vtk_read_polygons( FileTokenizer& tokens, MsqError& err )
{
  long size[2];
  
  tokens.get_long_ints( 2, size, err );                     MSQ_ERRRTN(err);
  tokens.get_newline( err );                                MSQ_ERRRTN(err);
  
  elementCount = size[0];
  Element* elem_ptr = elementArray = new Element[elementCount];

  for (int i = 0; i < size[0]; ++i)
  {
    long count;
    tokens.get_long_ints( 1, &count, err );                 MSQ_ERRRTN(err);
    switch (count)
    {
      case 3: elem_ptr->mType = TRIANGLE;      break;
      case 4: elem_ptr->mType = QUADRILATERAL; break;
      default:
        MSQ_SETERR(err)( MsqError::NOT_IMPLEMENTED,
                         "Unsupported polygon, not 3- or 4-sided at line %d",
                         tokens.line_number() );
        return;
    }
    
    size_t* vtx_idx_itr = elem_ptr->vertexIndices;
    for (int j = 0; j < count; ++j, ++vtx_idx_itr)
    {
      long val;
      tokens.get_long_ints( 1, &val, err);                  MSQ_ERRRTN(err);
      if (val < 0 || (unsigned long)val >= vertexCount)
      {
        MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                         "Vertex index out of bounds at line %d",
                         tokens.line_number() );
        return;
      }
      *vtx_idx_itr = val;
    }
  }
}



void MeshImpl::vtk_read_unstructured_grid( FileTokenizer& tokens, MsqError& err )
{
  long i, num_verts, num_elems[2];
  
  tokens.match_token( "POINTS", err );                      MSQ_ERRRTN(err);
  tokens.get_long_ints( 1, &num_verts, err );               MSQ_ERRRTN(err);
  tokens.match_token( vtk_type_names, err );                MSQ_ERRRTN(err);
  tokens.get_newline( err );                                MSQ_ERRRTN(err);
  
  if (num_verts < 1)
  {
    MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                     "Invalid point count at line %d", tokens.line_number());
    return;
  }
  
  
  vertexArray = new Vertex[num_verts];
  tokens.get_doubles( num_verts * 3, (double*)vertexArray, err );MSQ_ERRRTN(err);
  vertexCount = num_verts;
  
  tokens.match_token( "CELLS", err );                       MSQ_ERRRTN(err);
  tokens.get_long_ints( 2, num_elems, err );                MSQ_ERRRTN(err);
  tokens.get_newline( err );                                MSQ_ERRRTN(err);

  elementCount = num_elems[0];
  Element* elem_ptr = elementArray = new Element[elementCount];
  for (i = 0; i < num_elems[0]; ++i, ++elem_ptr)
  {
    long count;
    tokens.get_long_ints( 1, &count, err);                  MSQ_ERRRTN(err);
    
    if (count > MSQ_MAX_NUM_VERT_PER_ENT)
    {
      MSQ_SETERR(err)( MsqError::NOT_IMPLEMENTED,
                       "Too many vertices in element at line %d",
                       tokens.line_number() );
      return;
    }
    
    size_t* vtx_idx_itr = elem_ptr->vertexIndices;
    for (int j = 0; j < count; ++j, ++vtx_idx_itr)
    {
      long val;
      tokens.get_long_ints( 1,&val, err);                   MSQ_ERRRTN(err);
      
      if (val < 0 || (unsigned long)val >= vertexCount)
      {
        MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                         "Vertex index out of bounds at line %d",
                         tokens.line_number() );
        return;
      }
      *vtx_idx_itr = val;
    }
  }
  
  tokens.match_token( "CELL_TYPES", err );                  MSQ_ERRRTN(err);
  tokens.get_long_ints( 1, &num_elems[1], err );            MSQ_ERRRTN(err);
  tokens.get_newline( err );                                MSQ_ERRRTN(err);
  
  if (num_elems[0] != num_elems[1])
  {
    MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                    "Number of element types does not match number of elements"
                    "at line %d", tokens.line_number() );
    return;
  }
  
  elem_ptr = elementArray ;
  for (i = 0; i < num_elems[0]; ++i, ++elem_ptr)
  {
    long type;
    tokens.get_long_ints( 1, &type, err );                  MSQ_ERRRTN(err);
    switch (type)
    {
      case 5: 
        elem_ptr->mType = TRIANGLE; 
        break;
      case 8:
        std::swap( elem_ptr->vertexIndices[2], elem_ptr->vertexIndices[3] );
      case 9:
        elem_ptr->mType = QUADRILATERAL;
        break;
      case 10:
        elem_ptr->mType = TETRAHEDRON;
        break;
      case 11:
        std::swap( elem_ptr->vertexIndices[2], elem_ptr->vertexIndices[3] );
        std::swap( elem_ptr->vertexIndices[6], elem_ptr->vertexIndices[7] );
      case 12:
        elem_ptr->mType = HEXAHEDRON;
        break;
      case 13:
        std::swap( elem_ptr->vertexIndices[1], elem_ptr->vertexIndices[2] );
        std::swap( elem_ptr->vertexIndices[4], elem_ptr->vertexIndices[5] );
        elem_ptr->mType = PRISM;
        break;
      case 14:
        elem_ptr->mType = PYRAMID;
        break;
      default:
        MSQ_SETERR(err)( MsqError::NOT_IMPLEMENTED,
                         "Unsupported cell type (%ld) at line %d.",
                         type, tokens.line_number() );
        return;
    }
  }
}

void MeshImpl::vtk_create_structured_elems( const long* dims, 
                                            MsqError& err )
{
    //NOTE: this should be work fine for edges also if 
    //      Mesquite ever supports them.  Just add the
    //      type for dimension 1 to the switch statement.
    
  int non_zero[3] = {0,0,0};  // True if dim > 0 for x, y, z respectively
  long elem_dim = 0;          // Element dimension (2->quad, 3->hex)
  long num_elems = 1;         // Total number of elements
  long vert_per_elem;         // Element connectivity length
  long edims[3] = { 1, 1, 1 };// Number of elements in each grid direction
  
    // Populate above data
  for (int d = 0; d < 3; d++) 
    if (dims[d] > 1)
    {
      non_zero[elem_dim++] = d;
      edims[d] = dims[d] - 1;
      num_elems *= edims[d];
    }
  vert_per_elem = 1 << elem_dim;
  
    // Get element type from element dimension
  EntityTopology type;
  switch( elem_dim )
  {
  //case 1: type = EDGE;          break;
    case 2: type = QUADRILATERAL; break;
    case 3: type = HEXAHEDRON;    break;
    default:
      MSQ_SETERR(err)( "Cannot create structured mesh with elements "
                       "of dimension < 2 or > 3.",
                       MsqError::NOT_IMPLEMENTED );
      return;
  }

    // Allocate storage for elements
  Element* elem_ptr = elementArray = new Element[num_elems];
  
    // Offsets of element vertices in grid relative to corner closest to origin 
  long k = dims[0]*dims[1];
  const long corners[8] = { 0, 1, 1+dims[0], dims[0], k, k+1, k+1+dims[0], k+dims[0] };
                             
    // Populate element list
  for (long z = 0; z < edims[2]; ++z)
    for (long y = 0; y < edims[1]; ++y)
      for (long x = 0; x < edims[0]; ++x)
      {
        const long index = x + y*dims[0] + z*(dims[0]*dims[1]);
        elem_ptr->mType = type;
        for (long j = 0; j < vert_per_elem; ++j)
          elem_ptr->vertexIndices[j] = index + corners[j];
        ++elem_ptr;
      }
  
  elementCount = num_elems;
}

void MeshImpl::vtk_read_field( FileTokenizer& tokens, MsqError& err )
{
    // This is not supported yet.
    // Parse the data but throw it away because
    // Mesquite has no internal representation for it.
  
    // Could save this in tags, but the only useful thing that
    // could be done with the data is to write it back out
    // with the modified mesh.  As there's no way to save the
    // type of a tag in Mesquite, it cannot be written back
    // out correctly either.
    // FIXME: Don't know what to do with this data.
    // For now, read it and throw it out.

  long num_arrays;
  tokens.get_long_ints( 1, &num_arrays, err );              MSQ_ERRRTN(err);
  
  for (long i = 0; i < num_arrays; ++i)
  {
    /*const char* name =*/ tokens.get_string( err );        MSQ_ERRRTN(err);
    
    long dims[2];
    tokens.get_long_ints( 2, dims, err );                   MSQ_ERRRTN(err);
    tokens.match_token( vtk_type_names, err );              MSQ_ERRRTN(err);
    
    long num_vals = dims[0] * dims[1];
    
    for (long j = 0; j < num_vals; j++)
    {
      double junk;
      tokens.get_doubles( 1, &junk, err );                  MSQ_ERRRTN(err);
    }
  }
}

void* MeshImpl::vtk_read_attrib_data( FileTokenizer& tokens, 
                                      long count,
                                      TagDescription& tag,
                                      MsqError& err )
{
  const char* const type_names[] = { "SCALARS",
                                     "COLOR_SCALARS",
                                     "VECTORS",
                                     "NORMALS",
                                     "TEXTURE_COORDINATES",
                                     "TENSORS",
                                     "FIELD",
                                     0 };

  int type = tokens.match_token( type_names, err );
  const char* name = tokens.get_string( err ); MSQ_ERRZERO(err);
  tag.name = name;
  
  void* data = 0;
  switch( type )
  {
    case 1: data = vtk_read_scalar_attrib ( tokens, count, tag, err ); 
            tag.vtkType = TagDescription::SCALAR; 
            break;
    case 2: data = vtk_read_color_attrib  ( tokens, count, tag, err ); 
            tag.vtkType = TagDescription::COLOR;
            break;
    case 3: data = vtk_read_vector_attrib ( tokens, count, tag, err );
            tag.vtkType = TagDescription::VECTOR;
            break;
    case 4: data = vtk_read_vector_attrib ( tokens, count, tag, err ); 
            tag.vtkType = TagDescription::NORMAL;
            break;
    case 5: data = vtk_read_texture_attrib( tokens, count, tag, err ); 
            tag.vtkType = TagDescription::TEXTURE;
            break;
    case 6: data = vtk_read_tensor_attrib ( tokens, count, tag, err ); 
            tag.vtkType = TagDescription::TENSOR;
            break;
    case 7: // Can't handle field data yet.
      MSQ_SETERR(err)( MsqError::NOT_IMPLEMENTED,
                       "Cannot read field data (line %d).",
                       tokens.line_number());
  }

  return data;
}

void MeshImpl::vtk_read_point_data( FileTokenizer& tokens, 
                                    MsqError& err )
{
  TagDescription tag;
  void* data = vtk_read_attrib_data( tokens, vertexCount, tag, err );
  MSQ_ERRRTN(err);
  
  TagHandle handle = tag_get( tag.name, err );
  TagData* tag_ptr;
  if (!err)
  {
    tag_ptr = tag_from_handle( handle, err );  MSQ_ERRRTN(err);
    if (tag_ptr->desc != tag)
    {
      MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                       "Inconsistent types between element "
                       "and vertex attributes of same name "
                       "at line %d", tokens.line_number() );
      free( data );
      return;
    }
    
    if (tag_ptr->vertIsDefault)
    {
      free(tag_ptr->vertexData);
      tag_ptr->vertexData = 0;
      tag_ptr->vertIsDefault = false;
    }
    else if (tag_ptr->vertexData)
    {
      MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                       "Duplicate vertex attribute name"
                       " at line %d", tokens.line_number() );
      free( data );
      return;
    }
  }
  else
  {
    err.clear();
    tag_ptr = new TagData( tag );
    tagList.push_back(tag_ptr);
  }
  tag_ptr->vertexData = data;
}


void MeshImpl::vtk_read_cell_data( FileTokenizer& tokens, 
                                   MsqError& err )
{
  TagDescription tag;
  void* data = vtk_read_attrib_data( tokens, elementCount, tag, err );
  MSQ_ERRRTN(err);
  
  TagHandle handle = tag_get( tag.name, err );
  TagData* tag_ptr;
  if (!err)
  {
    tag_ptr = tag_from_handle( handle, err );  MSQ_ERRRTN(err);
    if (tag_ptr->desc != tag)
    {
      MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                       "Inconsistent types between element "
                       "and vertex attributes of same name "
                       "at line %d", tokens.line_number() );
      free( data );
      return;
    }
    
    if (tag_ptr->elemIsDefault)
    {
      free(tag_ptr->elementData);
      tag_ptr->elementData = 0;
      tag_ptr->elemIsDefault = false;
    }
    else if (tag_ptr->elementData)
    {
      MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                       "Duplicate element attribute name"
                       " at line %d", tokens.line_number() );
      free( data );
      return;
    }
  }
  else
  {
    err.clear();
    tag_ptr = new TagData( tag );
    tagList.push_back(tag_ptr);
  }
  tag_ptr->elementData = data;
}

void* MeshImpl::vtk_read_typed_data( FileTokenizer& tokens, 
                                     int type, 
                                     size_t per_elem, 
                                     size_t num_elem,
                                     TagDescription& tag,
                                     MsqError &err )
{
  void* data_ptr;
  size_t count = per_elem * num_elem;
  switch ( type )
  {
    case 1:
      tag.size = per_elem*sizeof(bool);
      tag.type = BOOL;
      data_ptr = malloc( num_elem*tag.size );
      tokens.get_booleans( count, (bool*)data_ptr, err );
      break;
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
      tag.size = per_elem*sizeof(int);
      tag.type = INT;
      data_ptr = malloc( num_elem*tag.size );
      tokens.get_integers( count, (int*)data_ptr, err );
      break;
    case 10:
    case 11:
      tag.size = per_elem*sizeof(double);
      tag.type = DOUBLE;
      data_ptr = malloc( num_elem*tag.size );
      tokens.get_doubles( count, (double*)data_ptr, err );
      break;
    default:
      MSQ_SETERR(err)( "Invalid data type", MsqError::INVALID_ARG );
      return 0;
  }
  
  if (MSQ_CHKERR(err))
  {
    free( data_ptr );
    return 0;
  }
  
  return data_ptr;
}
  
      
      

void* MeshImpl::vtk_read_scalar_attrib( FileTokenizer& tokens,
                                        long count,
                                        TagDescription& desc,
                                        MsqError& err )
{
  int type = tokens.match_token( vtk_type_names, err );      MSQ_ERRZERO(err);
    
  long size;
  const char* tok = tokens.get_string( err );                MSQ_ERRZERO(err);
  const char* end = 0;
  size = strtol( tok, (char**)&end, 0 );
  if (*end)
  {
    size = 1;
    tokens.unget_token();
  }
  
    // VTK spec says cannot be greater than 4--do we care?
  if (size < 1 || size > 4)
  {
    MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                    "Scalar count out of range [1,4]" 
                    " at line %d", tokens.line_number());
    return 0;
  }
  
  tokens.match_token("LOOKUP_TABLE",err);                     MSQ_ERRZERO(err);
  tok = tokens.get_string(err);                               MSQ_ERRZERO(err);
  
    // If no lookup table, just read and return the data
  if (!strcmp( tok, "default" ))
  {
    void* ptr = vtk_read_typed_data( tokens, type, size, count, desc, err );
    MSQ_ERRZERO(err);
    return ptr;
  }
  
    // If we got this far, then the data has a lookup
    // table.  First read the lookup table and convert
    // to integers.
  string name = tok;
  vector<long> table( size*count );
  if (type > 0 && type < 10)  // Is an integer type
  {
    tokens.get_long_ints( table.size(), &table[0], err );   
    MSQ_ERRZERO(err);
  }
  else // Is a real-number type
  {
    for (msq_std::vector<long>::iterator iter = table.begin(); iter != table.end(); ++iter)
    {
      double data;
      tokens.get_doubles( 1, &data, err );
      MSQ_ERRZERO(err);

      *iter = (long)data;
      if ((double)*iter != data)
      {
        MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                         "Invalid lookup index (%.0f) at line %d",
                         data, tokens.line_number() );
        return 0;
      }
    }
  }
  
    // Now read the data - must be float RGBA color triples
  
  long table_count;
  tokens.match_token( "LOOKUP_TABLE", err );                  MSQ_ERRZERO(err);
  tokens.match_token( name.c_str(), err );                    MSQ_ERRZERO(err);
  tokens.get_long_ints( 1, &table_count, err );               MSQ_ERRZERO(err);
 
  vector<float> table_data(table_count*4);
  tokens.get_floats( table_data.size(), &table_data[0], err );MSQ_ERRZERO(err);
  
    // Create list from indexed data
  
  float* data = (float*)malloc( sizeof(float)*count*size*4 );
  float* data_iter = data;
  for (std::vector<long>::iterator idx = table.begin(); idx != table.end(); ++idx)
  {
    if (*idx < 0 || *idx >= table_count)
    {
      MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                       "LOOKUP_TABLE index %ld out of range.",
                       *idx );
      free( data );
      return  0;
    }
    
    for (int i = 0; i < 4; i++)
    {
      *data_iter = table_data[4 * *idx + i];
      ++data_iter;
    }
  }
  
  desc.size = size * 4 * sizeof(float);
  desc.type = DOUBLE;
  return data;
}

void* MeshImpl::vtk_read_color_attrib( FileTokenizer& tokens, 
                                       long count, 
                                       TagDescription& tag,
                                       MsqError& err )
{
  long size;
  tokens.get_long_ints( 1, &size, err );                    MSQ_ERRZERO(err);
  
  if (size < 1)
  {
    MSQ_SETERR(err)( MsqError::PARSE_ERROR, 
                     "Invalid size (%ld) at line %d",
                     size, tokens.line_number() );
    return 0;
  }
  
  float* data = (float*)malloc( sizeof(float)*count*size );
  tokens.get_floats( count*size, data, err );
  if (MSQ_CHKERR(err))
  {
    free( data );
    return 0;
  }
  
  tag.size = size*sizeof(float);
  tag.type = DOUBLE;
  return data;
}

void* MeshImpl::vtk_read_vector_attrib( FileTokenizer& tokens, 
                                        long count, 
                                        TagDescription& tag,
                                        MsqError& err )
{
  int type = tokens.match_token( vtk_type_names, err );
  MSQ_ERRZERO(err);
    
  void* result = vtk_read_typed_data( tokens, type, 3, count, tag, err );
  MSQ_ERRZERO(err);
  return result;
}

void* MeshImpl::vtk_read_texture_attrib( FileTokenizer& tokens,
                                         long count,
                                         TagDescription& tag,
                                         MsqError& err )
{
  int type, dim;
  tokens.get_integers( 1, &dim, err );
  MSQ_ERRZERO(err);
  type = tokens.match_token( vtk_type_names, err );
  MSQ_ERRZERO(err);
    
  if (dim < 1 || dim > 3)
  {
    MSQ_SETERR(err)( MsqError::PARSE_ERROR,
                     "Invalid dimension (%d) at line %d.",
                     dim, tokens.line_number() );
    return 0;
  }
  
  void* result = vtk_read_typed_data( tokens, type, dim, count, tag, err );
  MSQ_ERRZERO(err);
  return result;
}

void* MeshImpl::vtk_read_tensor_attrib( FileTokenizer& tokens,
                                        long count, 
                                        TagDescription& tag,
                                        MsqError& err )
{
  int type = tokens.match_token( vtk_type_names, err );
  MSQ_ERRZERO(err);
    
  void* result = vtk_read_typed_data( tokens, type, 9, count, tag, err );
  MSQ_ERRZERO(err);
  return result;
}  

void MeshImpl::vtk_write_attrib_data( msq_stdio::ostream& file,
                                      const TagDescription& desc,
                                      const void* data, size_t count,
                                      MsqError& err ) const
{
  if (desc.type == HANDLE)
  {
    MSQ_SETERR(err)("Cannot write HANDLE tag data to VTK file.",
                    MsqError::FILE_FORMAT);
    return;
  }
    
  
  TagDescription::VtkType vtk_type = desc.vtkType;
  unsigned vlen = desc.size / size_from_tag_type(desc.type);
    // guess one from data length if not set
  if (vtk_type == TagDescription::NONE)
  {
    switch ( vlen )
    {
      case 1: 
      case 2: 
      case 4: vtk_type = TagDescription::SCALAR; break;
      case 3: vtk_type = TagDescription::VECTOR; break;
      case 9: vtk_type = TagDescription::TENSOR; break;
      default:
        MSQ_SETERR(err)(MsqError::FILE_FORMAT,
          "Cannot map vector tag \"%s\" with length %u to a VTK attribute type.",
          desc.name.c_str(), vlen );
        return;
    }
  }
  
  const char* const typenames[] = { "unsigned_char", "bit", "int", "double" };
  
  int num_per_line;
  switch (vtk_type)
  {
    case TagDescription::SCALAR: 
      num_per_line = vlen; 
      file << "SCALARS " << desc.name << " " << typenames[desc.type] << " " << vlen << "\n";
      break;
    case TagDescription::COLOR : 
      num_per_line = vlen;
      file << "COLOR_SCALARS " << desc.name << " " << vlen << "\n";
      break;
    case TagDescription::VECTOR:
      num_per_line = 3;
      if (vlen != 3)
      {
        MSQ_SETERR(err)(MsqError::INTERNAL_ERROR,
         "Tag \"%s\" is labeled as a VTK vector attribute but has %u values.",
         desc.name.c_str(), vlen);
        return;
      }
      file << "VECTORS " << desc.name << " " << typenames[desc.type] << "\n";
      break;
    case TagDescription::NORMAL:
      num_per_line = 3;
      if (vlen != 3)
      {
        MSQ_SETERR(err)(MsqError::INTERNAL_ERROR,
         "Tag \"%s\" is labeled as a VTK normal attribute but has %u values.",
         desc.name.c_str(), vlen);
        return;
      }
      file << "NORMALS " << desc.name << " " << typenames[desc.type] << "\n";
      break;
    case TagDescription::TEXTURE:
      num_per_line = vlen;
      file << "TEXTURE_COORDINATES " << desc.name << " " << typenames[desc.type] << " " << vlen << "\n";
      break;
    case TagDescription::TENSOR:
      num_per_line = 3;
      if (vlen != 9)
      {
        MSQ_SETERR(err)(MsqError::INTERNAL_ERROR,
         "Tag \"%s\" is labeled as a VTK tensor attribute but has %u values.",
         desc.name.c_str(), vlen);
        return;
      }
      file << "TENSORS " << desc.name << " " << typenames[desc.type] << "\n";
      break;
    default:
      MSQ_SETERR(err)("Unknown VTK attribute type for tag.", MsqError::INTERNAL_ERROR );
      return;
  }
  
  size_t i = 0, total = count*vlen;
  char* space = new char[num_per_line];
  memset( space, ' ', num_per_line );
  space[0] = '\n';
  const unsigned char* odata = (const unsigned char*)data;
  const bool* bdata = (const bool*)data;
  const int* idata = (const int*)data;
  const double* ddata = (const double*)data;
  switch ( desc.type )
  {
    case BYTE:
      while (i < total)
        file << (unsigned int)odata[i++] << space[i%num_per_line];
      break;
    case BOOL:
      while (i < total)
        file << (bdata[i++] ? '1' : '0') << space[i%num_per_line];
      break;
    case INT:
      while (i < total)
        file << idata[i++] << space[i%num_per_line];
      break;
    case DOUBLE:
      while (i < total)
        file << ddata[i++] << space[i%num_per_line];
      break;
    default:
      MSQ_SETERR(err)("Unknown tag type.", MsqError::INTERNAL_ERROR);
  }
  delete [] space;
}
      
        
          
  


/**************************************************************************
 *                               TAGS
 **************************************************************************/


MeshImpl::TagData::~TagData() 
{
  if (elementData) 
    free(elementData);
  if (vertexData)
    free(vertexData);
}


MeshImpl::TagData* MeshImpl::tag_from_handle( TagHandle handle, MsqError& err )
{
  unsigned index = (unsigned) handle;
  if (index >= tagList.size() || !tagList[index])
  {
    MSQ_SETERR(err)("Invalid tag handle.",MsqError::INVALID_ARG);
    return 0;
  }
  
  return tagList[index];
}

size_t MeshImpl::size_from_tag_type( TagType type )
{
  switch( type ) {
    case BYTE:   return 1;
    case BOOL:   return sizeof(bool);
    case DOUBLE: return sizeof(double);
    case INT:    return sizeof(int);
    case HANDLE: return sizeof(EntityHandle);
    default: assert(0); return 0;
  }
}

TagHandle MeshImpl::tag_create( const string& name,
                                TagType type,
                                unsigned length,
                                const void* defval,
                                MsqError& err )
{
  TagHandle handle = tag_get( name, err );
  if (!err)
  {
    err.clear();
    MSQ_SETERR(err)(name, MsqError::TAG_ALREADY_EXISTS);
    return 0;
  }
  
  if (length == 0 || size_from_tag_type(type) == 0)
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    return 0;
  }
  
  TagData* tag = new TagData( name, type, length );
  handle = (TagHandle)(tagList.size());
  tagList.push_back(tag);
  
  if (defval)
  {
    tag->vertIsDefault = tag->elemIsDefault = true;
    tag->elementData = malloc( tag->desc.size );
    tag->vertexData = malloc( tag->desc.size );
    memcpy( tag->elementData, defval, tag->desc.size );
    memcpy( tag->vertexData, defval, tag->desc.size );
  }
  
  return handle;
}

void MeshImpl::tag_delete( TagHandle handle, MsqError& err )
{
  unsigned index = (unsigned)handle;
  if (index >= tagList.size() || 0 == tagList[index])
  {
    MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
    return ;
  }
  
  delete tagList[index];
  tagList[index] = 0;
}

TagHandle MeshImpl::tag_get( const msq_std::string& name, MsqError& err )
{
  for (unsigned i = 0; i < tagList.size(); ++i)
    if (tagList[i] && tagList[i]->desc.name == name)
      return (TagHandle)i;
      
  MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
  return 0;
}

void MeshImpl::tag_properties( TagHandle handle,
                               msq_std::string& name,
                               TagType& type,
                               unsigned& length,
                               MsqError& err )
{
  TagData* data = tag_from_handle(handle,err); MSQ_ERRRTN(err);
  name = data->desc.name;
  type = data->desc.type;
  length = data->desc.size / size_from_tag_type(data->desc.type);
}


void MeshImpl::tag_set_element_data( TagHandle handle,
                                     size_t num_elems,
                                     const ElementHandle* elem_array,
                                     const void* values,
                                     MsqError& err )
{
  size_t i;
  TagData* tag = tag_from_handle(handle,err); MSQ_ERRRTN(err);
  
    // Allocate space for data, if it has not yet been allocated
  if (tag->elemIsDefault)
  {
    tag->elemIsDefault = false;
    void* def = tag->elementData;
    tag->elementData = malloc( tag->desc.size * elementCount );
    for (i = 0; i < elementCount; ++i)
      memcpy( ((char*)tag->elementData) + i*tag->desc.size, def, tag->desc.size );
    free( def );
  }
  else if (!tag->elementData)
  {
    tag->elementData = malloc( tag->desc.size * elementCount );
    memset( tag->elementData, 0, tag->desc.size * elementCount );
  }
  
    // Store passed tag values
  char* data = (char*)tag->elementData;
  const char* iter = (const char*)values;
  for (size_t i = 0; i < num_elems; ++i)
  {
    size_t index = ((Element*)elem_array[i]) - elementArray;
    if (index > elementCount)
    {
      MSQ_SETERR(err)(MsqError::INVALID_ARG);
      return;
    }
    
    memcpy( data + index*tag->desc.size, iter, tag->desc.size );
    iter += tag->desc.size;
  }
}

void MeshImpl::tag_get_element_data( TagHandle handle,
                                     size_t num_elems,
                                     const ElementHandle* elem_array,
                                     void* values,
                                     MsqError& err )
{
  TagData* tag = tag_from_handle(handle,err); MSQ_ERRRTN(err);
  char* iter = (char*)values;
  const char* data = (const char*)tag->elementData;
  
    // If no data (not even a default value) is set, error
  if (!data)
  {
    MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
    return;
  }
    // If element data has default value, return that for all elements
  else if (tag->elemIsDefault)
  {
    for (size_t i = 0; i < num_elems; ++i)
      memcpy( iter, data, tag->desc.size );
  }
    // Otherwise return specific value for each requested entity
  else
  {
    for (size_t i = 0; i < num_elems; ++i)
    {
      size_t index = ((Element*)elem_array[i]) - elementArray;
      if (index > elementCount)
      {
        MSQ_SETERR(err)(MsqError::INVALID_ARG);
        return;
      }

      memcpy( iter, data + index*tag->desc.size, tag->desc.size );
      iter += tag->desc.size;
    }
  }
}

void MeshImpl::tag_set_vertex_data(  TagHandle handle,
                                     size_t num_elems,
                                     const VertexHandle* elem_array,
                                     const void* values,
                                     MsqError& err )
{
  size_t i;
  TagData* tag = tag_from_handle(handle,err); MSQ_ERRRTN(err);
  
  if (tag->vertIsDefault)
  {
    tag->vertIsDefault = false;
    void* def = tag->vertexData;
    tag->vertexData = malloc( tag->desc.size * vertexCount );
    for (i = 0; i < vertexCount; ++i)
      memcpy( ((char*)tag->vertexData) + i*tag->desc.size, def, tag->desc.size );
    free( def );
  }
  else if (!tag->vertexData)
  {
    tag->vertexData = malloc( tag->desc.size * vertexCount );
    memset( tag->vertexData, 0, tag->desc.size * vertexCount );
  }
  
  char* data = (char*)tag->vertexData;
  const char* iter = (const char*)values;
  for (size_t i = 0; i < num_elems; ++i)
  {
    size_t index = ((Vertex*)elem_array[i]) - vertexArray;
    if (index > vertexCount)
    {
      MSQ_SETERR(err)(MsqError::INVALID_ARG);
      return;
    }
    
    memcpy( data + index*tag->desc.size, iter, tag->desc.size );
    iter += tag->desc.size;
  }
}

void MeshImpl::tag_get_vertex_data(  TagHandle handle,
                                     size_t num_elems,
                                     const VertexHandle* elem_array,
                                     void* values,
                                     MsqError& err )
{
  TagData* tag = tag_from_handle(handle,err); MSQ_ERRRTN(err);
  char* iter = (char*)values;
  const char* data = (const char*)tag->vertexData;
  
    // If no data (not even a default value) is set, error
  if (!data)
  {
    MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND);
    return;
  }
    // If vertex data has default value, return that for all vertices
  else if (tag->vertIsDefault)
  {
    for (size_t i = 0; i < num_elems; ++i)
      memcpy( iter, data, tag->desc.size );
  }
    // Otherwise return specific value for each requested entity
  else
  {
    for (size_t i = 0; i < num_elems; ++i)
    {
      size_t index = ((Vertex*)elem_array[i]) - vertexArray;
      if (index > vertexCount)
      {
        MSQ_SETERR(err)(MsqError::INVALID_ARG);
        return;
      }

      memcpy( iter, data + index*tag->desc.size, tag->desc.size );
      iter += tag->desc.size;
    }
  }
}

bool MeshImpl::tag_has_vertex_data( TagHandle handle, MsqError& err ) 
{
  TagData* tag = tag_from_handle(handle,err); MSQ_ERRZERO(err);
  return 0 != tag->vertexData;
}  

bool MeshImpl::tag_has_element_data( TagHandle handle, MsqError& err ) 
{
  TagData* tag = tag_from_handle(handle,err); MSQ_ERRZERO(err);
  return 0 != tag->elementData;
}  


} // namespace Mesquite


