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
/*!
  \file   ParallelHelper.hpp
  \brief  

  Implements ParallelHelper Class 

  \author Martin Isenburg
  \date   2008-03-04
*/

#ifndef Mesquite_ParallelHelper_hpp 
#define Mesquite_ParallelHelper_hpp

#include "ParallelHelperInterface.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#include <vector.h>
#include <map.h>
#else
#include <vector>
#include <map>
#endif

#include <mpi.h>

namespace MESQUITE_NS
{
  typedef struct VertexIdMapKey {
    int id;
    int proc_id;
  } VertexIdMapKey;
  
  struct VertexIdLessFunc {
    bool operator()( const VertexIdMapKey &that1, const VertexIdMapKey& that2 ) const
    {
        return ( (that1.proc_id < that2.proc_id) || ((that1.proc_id==that2.proc_id)&&(that1.id<that2.id)) );
    }
  };
  
    typedef std::map<VertexIdMapKey,int,VertexIdLessFunc> VertexIdMap;

  class ParallelHelperImpl : public ParallelHelper
  {
  public:

    enum CommunicationModel{ TrulyNonBlocking = 0, TrulyNonBlockingAvoidAllReduce,
			     NonBlocking, NonBlockingAvoidAllReduce,
			     Blocking, BlockingAvoidAllReduce };

    ParallelHelperImpl();
    ~ParallelHelperImpl();

    // function called by application during set-up
    bool set_parallel_mesh(ParallelMesh* mesh);
    bool set_communicator(size_t comm);
    bool set_communication_model(int model);
    bool set_generate_random_numbers(int grn);

  protected:
    friend class VertexMover;
    // functions called by VertexMover::loop_over_mesh()
    bool smoothing_init();
    void compute_first_independent_set(msq_std::vector<Mesh::VertexHandle>& fixed_vertices);
    void communicate_first_independent_set();
    bool compute_next_independent_set();
    bool get_next_partition_boundary_vertex(Mesquite::Mesh::VertexHandle& vertex_handle);
    void communicate_next_independent_set();
    bool smoothing_close();

  protected:
    friend class QualityAssessor;
    // functions called by QualityAssessor::loop_over_mesh()
    int get_rank() const;
    int get_nprocs() const;
    bool is_our_element(Mesquite::Mesh::ElementHandle element_handle) const;
    bool is_our_vertex(Mesquite::Mesh::VertexHandle vertex_handle) const;
    void communicate_min_max_to_all(double* minimum, double* maximum) const;
    void communicate_min_max_to_zero(double* minimum, double* maximum) const;
    void communicate_sums_to_zero(size_t* freeElementCount, int* invertedElementCount, size_t* elementCount, int* invertedSampleCount, size_t* sampleCount, long unsigned int* count, long unsigned int* invalid, double* sum, double *sqrSum) const;
    void communicate_power_sum_to_zero(double* pMean) const;
    void communicate_histogram_to_zero(msq_std::vector<int> &histogram) const;

  private:
    ParallelMesh* mesh;

    MPI_Comm communicator;
    int communication_model;

    int rank;
    int nprocs;

    // variables for VertexMover::loop_over_mesh()
    int generate_random_numbers;
    msq_std::vector<Mesquite::Mesh::VertexHandle> *vertices;
    int num_vertex;
    char* vtx_in_partition_boundary;
    int num_vtx_partition_boundary;
    int num_vtx_partition_boundary_local;
    int num_vtx_partition_boundary_remote;
    msq_std::vector<Mesquite::Mesh::VertexHandle> *part_vertices;
    int* part_proc_owner;
    int* part_gid;
    int* part_smoothed_flag;
    double* part_rand_number;
    int num_exportVtx;
    int* exportVtxGIDs;
    int* exportVtxLIDs;
    int* exportProc;
    bool* in_independent_set;
    VertexIdMap* vid_map;
    int total_num_vertices_to_smooth;
    int total_num_vertices_to_recv;
    int* neighbourProcSend;
    int* neighbourProcRecv;
    int* neighbourProcSendRemain;
    int* neighbourProcRecvRemain;
    int num_already_smoothed_vertices;
    int num_already_recv_vertices;
    int* vtx_off_proc_list_size;
    int** vtx_off_proc_list;
    int num_neighbourProc;
    int* neighbourProc;
    int iteration;
    int global_work_remains;
    int next_vtx_partition_boundary;
    /* for exchanging unused ghost node information */
    int unghost_num_vtx;
    msq_std::vector<Mesquite::Mesh::VertexHandle> *unghost_vertices;
    int unghost_num_procs;
    int* unghost_procs;
    int* unghost_procs_num_vtx;
    int *unghost_procs_offset;
    int update_num_vtx;
    int* update_gid;
    int update_num_procs;
    int* update_procs;
    int* update_procs_num_vtx;
    int* update_procs_offset;

    // functions for VertexMover::loop_over_mesh()
    void compute_independent_set();
    int comm_smoothed_vtx_b();
    int comm_smoothed_vtx_b_no_all();
    int comm_smoothed_vtx_nb();
    int comm_smoothed_vtx_nb_no_all();
    int comm_smoothed_vtx_tnb();
    int comm_smoothed_vtx_tnb_no_all();
  };
  
} // namespace
#endif // Mesquite_ParallelHelper_hpp
