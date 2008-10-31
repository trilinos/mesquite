#include "ParallelHelper.hpp"

#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>

#include "MsqVertex.hpp"
#include "MsqError.hpp"

#define VERTEX_HEADER 1
#define VERTEX_BLOCK 1000

namespace Mesquite {

static void vertex_hash_insert(VertexIdHash* hash, int id, int proc_id, int value)
{
  VertexIdHashKey vid;
  vid.id = id;
  vid.proc_id = proc_id;
  hash->insert(VertexIdHash::value_type(vid, value));
}

static int vertex_hash_find(VertexIdHash* hash, int id, int proc_id)
{
  VertexIdHashKey vid;
  vid.id = id;
  vid.proc_id = proc_id;
  VertexIdHash::iterator hash_element = hash->find(vid);
  if (hash_element == hash->end())
  {
    return 0;
  }
  else
  {
    return (*hash_element).second;
  }
}

static void my_quicksort(int* a, int i, int j)
{
  int in_i = i;
  int in_j = j;
  int w;
  int key = a[(i+j)/2];
  do
    {
      while ( a[i] < key ) i++;
      while ( a[j] > key ) j--;
      if (i<j)
	{
	  w = a[i];
	  a[i] = a[j];
	  a[j] = w;
	}
    } while (++i<=--j);
  if (i == j+3)
    {
      i--;
      j++;
    }
  if (j>in_i) my_quicksort(a, in_i, j);
  if (i<in_j) my_quicksort(a, i, in_j);
}

static int hash6432shift(unsigned long long key)
{
  key = (~key) + (key << 18); // key = (key << 18) - key - 1;
  key = key ^ (key >> 31);
  key = key * 21; // key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (int) key;
}

static unsigned long long hash64shift(unsigned long long key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

static double generate_random_number(int generate_random_numbers, int proc_id, int id)
{
  if (generate_random_numbers == 1)
  {
    // count number of on bits
    int on = 0;
    int mist = id;
    while (mist)
    {
      if (mist & 1) on++;
      mist = mist >> 1;
    }
    unsigned short xsubi[3];
    ((int*)&(xsubi[0]))[0] = on&1 ? id : -id;
    xsubi[2] = proc_id;
    return erand48(xsubi);
  }
  else if (generate_random_numbers == 2)
  {
    // count number of on bits
    int on = 0;
    int mist = id;
    while (mist)
    {
      if (mist & 1) on++;
      mist = mist >> 1;
    }
    unsigned short xsubi[3];
    ((int*)&(xsubi[0]))[0] = on&1 ? id : -id;
    xsubi[2] = proc_id ^ xsubi[1];
    return erand48(xsubi);
  }
  else if (generate_random_numbers == 3)
  {
    unsigned long long key = id;
    key = key << 32;
    key = key | proc_id;
    return (double)hash64shift(key);
  }
  else if (generate_random_numbers == 4)
  {
    unsigned long long key = id;
    key = key << 32;
    key = key | proc_id;
    key = hash64shift(key);
    unsigned short xsubi[3];
    xsubi[0] = (unsigned short)(key >> 48);
    xsubi[1] = (unsigned short)(key >> 32);
    xsubi[2] = (unsigned short)(key >> 16);
    return erand48(xsubi);
  }
  else
  {
    return -1;
  }
}

ParallelHelperImpl::ParallelHelperImpl()
{
  mesh = 0;
  communicator = MPI_COMM_WORLD;
  communication_model = TrulyNonBlockingAvoidAllReduce;

  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &nprocs);

  // variable defaults for VertexMover::loop_over_mesh()
  generate_random_numbers = 2;

  // memory chunks referenced for VertexMover::loop_over_mesh()
  vertices = 0;
  vtx_in_partition_boundary = 0;
  part_vertices = 0;
  unused_ghost_vertices = 0;
  part_proc_owner = 0;
  part_gid = 0;
  part_smoothed_flag = 0;
  part_rand_number = 0;
  exportVtxGIDs = 0;
  exportVtxLIDs = 0;
  exportProc = 0;
  in_independent_set = 0;
  vid_hash = 0;
  neighbourProcSend = 0;
  neighbourProcRecv = 0;
  neighbourProcSendRemain = 0;
  neighbourProcRecvRemain = 0;
  vtx_off_proc_list_size = 0;
  vtx_off_proc_list = 0;
  neighbourProc = 0;
}

ParallelHelperImpl::~ParallelHelperImpl()
{
  // just making sure
  smoothing_close();
}

bool ParallelHelperImpl::set_parallel_mesh(ParallelMesh* mesh) {
  this->mesh = mesh;
  return true;
}

bool ParallelHelperImpl::set_communicator(size_t comm) {
  communicator = (MPI_Comm)comm;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &nprocs);
  return true;
}

bool ParallelHelperImpl::set_communication_model(int model) {
  communication_model = model;
  return true;
}

bool ParallelHelperImpl::set_generate_random_numbers(int grn) {
  this->generate_random_numbers = grn;
  return true;
}

bool ParallelHelperImpl::smoothing_init()
{
  int i,j,k,l;
  Mesquite::MsqError err;

  if (!mesh) return false;

  /* get the vertices */
  vertices = new msq_std::vector<Mesquite::Mesh::VertexHandle>;
  mesh->get_all_vertices(*vertices, err);
  if (err) {std::cout << err << std::endl; return false;}
  num_vertex = vertices->size();

  /* allocate the data arrays we'll use for smoothing */
  int* gid = new int[num_vertex];
  bool* app_fixed = new bool[num_vertex];
  int* proc_owner = new int[num_vertex];

  /* get the data from the mesquite mesh */
  mesh->vertices_get_global_id(&(*vertices)[0],&gid[0],num_vertex,err);
  if (err) {std::cout << err << std::endl; return 1; }

  mesh->vertices_get_fixed_flag(&(*vertices)[0],&app_fixed[0],num_vertex,err);
  if (err) {std::cout << err << std::endl; return 1; }

  mesh->vertices_get_processor_id(&(*vertices)[0],&proc_owner[0],num_vertex,err);
  if (err) {std::cout << err << std::endl; return 1; }

  /* create temporary Tag for the local IDs */
  int* lid = new int[num_vertex];
  for (i=0; i < num_vertex; i++) lid[i] = i;
  const char LOCAL_ID_NAME[] = "LOCAL_ID";
  TagHandle lid_tag = mesh->tag_create( LOCAL_ID_NAME, Mesh::INT, 1, NULL, err );
  mesh->tag_set_vertex_data( lid_tag, num_vertex, &(*vertices)[0], lid, err );

  /* get the elements */
  msq_std::vector<Mesquite::Mesh::ElementHandle> *elements = new msq_std::vector<Mesquite::Mesh::ElementHandle>;
  mesh->get_all_elements(*elements, err);
  if (err) {std::cout << err << std::endl; return 1; }
  int num_elems = elements->size();

  /****************************************************************
  PARTITION BOUNDARY VERTEX DETERMINATION
  ***********************************************************************/

  /* initialize the vertex partition boundary array */
  vtx_in_partition_boundary = new char[num_vertex];
  for (i=0;i<num_vertex;i++) vtx_in_partition_boundary[i] = 0;
  int incident_vtx, vtx_off_proc, vtx_on_proc;

  /* get the array that contains the adjacent vertices for each mesh element */
  msq_std::vector<Mesquite::Mesh::VertexHandle> *adj_vertices = new msq_std::vector<Mesquite::Mesh::VertexHandle>;
  msq_std::vector<size_t> *vtx_offsets = new msq_std::vector<size_t>;
  mesh->elements_get_attached_vertices(&(*elements)[0],num_elems,*adj_vertices,*vtx_offsets,err);
  delete elements; elements = 0;
  int* adj_vertices_lid = new int[(*adj_vertices).size()];
  mesh->tag_get_vertex_data( lid_tag, (*adj_vertices).size(), &(*adj_vertices)[0], adj_vertices_lid, err );
  delete adj_vertices; adj_vertices = 0;

  /* determine which vertices are smoothed as part of the boundary */
  num_vtx_partition_boundary = 0;
  num_vtx_partition_boundary_local = 0;
  num_vtx_partition_boundary_remote = 0;
  for (i=0;i<num_elems;i++) {
    /* count how many vertices of the current element are on/off a different processor */
    vtx_off_proc = 0;
    vtx_on_proc = 0;
    for (j=(*vtx_offsets)[i];j<(*vtx_offsets)[i+1];j++) {
	incident_vtx = adj_vertices_lid[j];
	/* obviously the vertex only counts if it is not app_fixed */
	if (!app_fixed[incident_vtx]) {
	  if (proc_owner[incident_vtx]!=rank) {
	    vtx_off_proc++;
	  }
	  else {
	    vtx_on_proc++;
	  }
	}
    }

    /* if vertices are on different processors mark all local vertices as boundary1 and all remote vertices as boundary2 */
    if (vtx_off_proc > 0 && vtx_on_proc > 0) {
	/* collect stats */
	//	smooth_stats.num_part_bndy_elem++;
	/* mark the vertices */
	for (j=(*vtx_offsets)[i];j<(*vtx_offsets)[i+1];j++) {
	  incident_vtx = adj_vertices_lid[j];
	  /* obviously the vertex does not need to be marked if it was already marked or if it is app_fixed*/
	  if (vtx_in_partition_boundary[incident_vtx] <= 0 && app_fixed[incident_vtx] == 0) {
	    /* collect stats */
	    //	      smooth_stats.num_part_bndy_vtx++;
	    /* mark and count the vertex */
	    if (proc_owner[incident_vtx]!=rank) {
	      vtx_in_partition_boundary[incident_vtx] = 2;
	      num_vtx_partition_boundary_remote++;
	    }
	    else {
	      vtx_in_partition_boundary[incident_vtx] = 1;
	      num_vtx_partition_boundary_local++;
	    }
	  }
	}
    }
    else if (vtx_off_proc > 0) {
      /* mark the vertices as boundary-1 if the element has only off-processor vertices */
	for (j=(*vtx_offsets)[i];j<(*vtx_offsets)[i+1];j++) {
	  incident_vtx = adj_vertices_lid[j];
	  /* obviously the vertex is not marked if it was already marked or if it is app_fixed*/
	  if (vtx_in_partition_boundary[incident_vtx] == 0 && app_fixed[incident_vtx] == 0) {
	    vtx_in_partition_boundary[incident_vtx] = -1;
	  }
	}
    }
  }    

  if (0)
  {
    printf("[%d]i%d local %d remote %d ",rank,iteration,num_vtx_partition_boundary_local,num_vtx_partition_boundary_remote);
    printf("[%d]i%d pb1 ",rank,iteration);
    for (i=0;i<num_vertex;i++) if (vtx_in_partition_boundary[i] == 1) printf("%d,%d ",i,gid[i]);
    printf("\n");
    printf("[%d]i%d pb2 ",rank,iteration);
    for (i=0;i<num_vertex;i++) if (vtx_in_partition_boundary[i] == 2) printf("%d,%d ",i,gid[i]);
    printf("\n");
    fflush(NULL);
  }

  num_vtx_partition_boundary = num_vtx_partition_boundary_local + num_vtx_partition_boundary_remote;

  delete [] app_fixed; app_fixed = 0;
  delete [] adj_vertices_lid; adj_vertices_lid = 0;
  delete vtx_offsets; vtx_offsets = 0;

  /********************************************************************
 COLLECT THE PARTITION BOUNDARY VERTICES
  ********************************************************************/

  /* create the vectors to stores the partition boundary vertex data */
  part_vertices = new msq_std::vector<Mesquite::Mesh::VertexHandle>;
  part_vertices->resize(num_vtx_partition_boundary);
  part_proc_owner = new int[num_vtx_partition_boundary];
  part_gid = new int[num_vtx_partition_boundary];
  part_smoothed_flag = new int[num_vtx_partition_boundary];
  part_rand_number = new double[num_vtx_partition_boundary];

  /* create the partition boundary map and its inverse */
  int* vtx_partition_boundary_map_inverse = new int[num_vertex];

  j=0;
  /* first we map the partition boundary vertices that we will smooth on this processor */
  for (i=0;i<num_vertex;i++) {
    if (vtx_in_partition_boundary[i]==1) {
	(*part_vertices)[j] = (*vertices)[i];
	part_proc_owner[j] = rank; assert(proc_owner[i] == rank);
	part_gid[j] = gid[i];
	vtx_partition_boundary_map_inverse[i] = j;
	j++;
    }
  }

  vid_hash = new VertexIdHash;

  /* then we map the ghost vertices that will be smoothed on other processors */
  for (i=0;i<num_vertex;i++) {
    if (vtx_in_partition_boundary[i]==2) {
	(*part_vertices)[j] = (*vertices)[i];
	part_proc_owner[j] = proc_owner[i];  assert(proc_owner[i] != rank);
	part_gid[j] = gid[i];
	vtx_partition_boundary_map_inverse[i] = j;
	/* only insert those vertices in the hash that are smoothed on other processors */
	vertex_hash_insert(vid_hash, part_gid[j], part_proc_owner[j], j);
	j++;
    }
  }

  /* count the number of un-used ghost vertices */
  int num_unused_ghost_vtx = 0;
  for (i=0;i<num_vertex;i++) {
    if (vtx_in_partition_boundary[i]==-1) {
      num_unused_ghost_vtx++;
    }
  }

  if (0 && num_unused_ghost_vtx) {printf("[%d] found %d unused ghost vertices (total %d local %d remote %d) \n",rank, num_unused_ghost_vtx, num_vertex, num_vtx_partition_boundary_local,num_vtx_partition_boundary_remote); fflush(NULL);}

  if (num_unused_ghost_vtx) {
    unused_ghost_vertices = new msq_std::vector<Mesquite::Mesh::VertexHandle>;
    unused_ghost_vertices->resize(num_unused_ghost_vtx);
    j = 0;
    for (i=0;i<num_vertex;i++) {
      if (vtx_in_partition_boundary[i]==-1) {
	(*unused_ghost_vertices)[j] = (*vertices)[i];
	j++;
      }
    }
  }

  /* no longer needed */
  delete [] gid; gid = 0;
  delete [] proc_owner; proc_owner = 0;

  /* create our own 'very pseudo random' numbers */
  for (i=0;i<num_vtx_partition_boundary;i++)
    part_rand_number[i] = generate_random_number(generate_random_numbers, part_proc_owner[i], part_gid[i]);

  /***********************************************************************
 COMPUTE THE SET OF NEIGHBORS THAT OUR VERTICES HAVE ON OTHER PROCESSORS
               COMPUTE THE SET OF NEIGHBORS PROCESSORS
  ***********************************************************************/

  if (num_vtx_partition_boundary_local == 0)
  {
    /* this processor does not partake in the boundary smoothing */
    num_neighbourProc = 0;
    mesh->tag_delete( lid_tag, err );
    delete [] vtx_partition_boundary_map_inverse;
    return true;
  }

  /* init the neighbour processor list */
  num_neighbourProc = 0;
  neighbourProc = (int*)malloc(sizeof(int)*10);
  for (i=0;i<9;i++) {
    neighbourProc[i] = 0;
  }
  neighbourProc[i] = -1;


  /* init the neighbour lists */

  int list_size;
  vtx_off_proc_list_size = new int[num_vtx_partition_boundary_local];
  vtx_off_proc_list = new int*[num_vtx_partition_boundary_local];
  for (i=0;i<num_vtx_partition_boundary_local;i++) {
    vtx_off_proc_list_size[i] = 0;
    vtx_off_proc_list[i] = (int*)malloc(sizeof(int)*15);
    for (j=0;j<14;j++) {
	vtx_off_proc_list[i][j] = 0;
    }
    vtx_off_proc_list[i][j] = -1;
  }

  /* get the adjacency arrays that we need */
  msq_std::vector<Mesquite::Mesh::ElementHandle> *adj_elements = new msq_std::vector<Mesquite::Mesh::ElementHandle>;
  msq_std::vector<Mesquite::Mesh::VertexHandle> *adj_adj_vertices = new msq_std::vector<Mesquite::Mesh::VertexHandle>;
  msq_std::vector<size_t> *elem_offsets = new msq_std::vector<size_t>;
  msq_std::vector<size_t> *adj_vtx_offsets = new msq_std::vector<size_t>;
  mesh->vertices_get_attached_elements(&(*part_vertices)[0],num_vtx_partition_boundary_local,
					 *adj_elements,*elem_offsets,err);
  mesh->elements_get_attached_vertices(&(*adj_elements)[0],adj_elements->size(),
					 *adj_adj_vertices,*adj_vtx_offsets,err);
  delete adj_elements; adj_elements = 0;
  int* adj_adj_vertices_lid = new int[(*adj_adj_vertices).size()];
  mesh->tag_get_vertex_data( lid_tag, (*adj_adj_vertices).size(), &(*adj_adj_vertices)[0], adj_adj_vertices_lid, err );
  delete adj_adj_vertices; adj_adj_vertices = 0;
  mesh->tag_delete( lid_tag, err );

  for (i=0;i<num_vtx_partition_boundary_local;i++) {
    /* loop over the elements surrounding that vertex */
    for (j=(*elem_offsets)[i];j<(*elem_offsets)[i+1];j++) {
      /* loop over the neighbors of the considered vertex (i.e. the vertices of these element) */
      for (k=(*adj_vtx_offsets)[j];k<(*adj_vtx_offsets)[j+1];k++) {
	/* get the next neighbour */
	incident_vtx = adj_adj_vertices_lid[k];
	/* if this neighbour is a vertex that is smoothed on a different processor */
	if (vtx_in_partition_boundary[incident_vtx] == 2) {
	  /* then map it into our domain */
	  incident_vtx = vtx_partition_boundary_map_inverse[incident_vtx];
	  /* is this vertex already in our neighbour list ? */
	  list_size = vtx_off_proc_list_size[i];
	  /* check by scanning the list for this vertex */
	  for (l = 0; l < list_size; l++) {
	    if (vtx_off_proc_list[i][l] == incident_vtx) {
	      /* no need to add this vertex to the list ... it is already there */
	      incident_vtx = -1;
	      /* end the loop */
	      l = list_size;
	    }
	  }
	  if (incident_vtx != -1) {
	    /* if the vertex is not in the list yet ... add it */
	    if (vtx_off_proc_list[i][list_size] == -1) {
	      /* need to make the list longer */
	      vtx_off_proc_list[i] = (int*)realloc(vtx_off_proc_list[i],sizeof(int)*list_size*2);
	      for (l=list_size;l<list_size*2-1;l++) {
		vtx_off_proc_list[i][l] = 0;
	      }
	      vtx_off_proc_list[i][l] = -1;
	    }
	    vtx_off_proc_list[i][list_size] = incident_vtx;
	    vtx_off_proc_list_size[i]++;
	    /* is the processor of this vertex already in the processor list */
	    incident_vtx = part_proc_owner[incident_vtx];
	    /* check by scanning the list for this processor */
	    for (l = 0; l < num_neighbourProc; l++) {
	      if (neighbourProc[l] == incident_vtx) {
		/* no need to add this processor to the list ... it is already there */
		incident_vtx = -1;
		/* end the loop */
		l = num_neighbourProc;
	      }
	    }
	    if (incident_vtx != -1) {
	      /* the processor is not in the list yet ... add it */
	      if (neighbourProc[num_neighbourProc] == -1) {
		/* need to make the list longer */
		neighbourProc = (int*)realloc(neighbourProc,sizeof(int)*num_neighbourProc*2);
		for (l=num_neighbourProc;l<num_neighbourProc*2-1;l++) {
		  neighbourProc[l] = 0;
		}
		neighbourProc[l] = -1;
	      }
	      neighbourProc[num_neighbourProc] = incident_vtx;
	      num_neighbourProc++;
	    }
	  }
	}
      }
    }
  }

  /* sort the list of neighbour processors */

  my_quicksort(neighbourProc,0,num_neighbourProc-1);

  delete [] vtx_partition_boundary_map_inverse;
  delete [] adj_adj_vertices_lid; adj_adj_vertices_lid = 0;
  delete elem_offsets; elem_offsets = 0;
  delete adj_vtx_offsets; adj_vtx_offsets = 0;

  /***********************************************************************
    COMPUTE HOW MANY VERTICES WE NEED TO SEND/RECV FROM EACH PROCESSOR
  ***********************************************************************/
  neighbourProcSend = 0;
  neighbourProcRecv = 0;
  neighbourProcSendRemain = 0;
  neighbourProcRecvRemain = 0;

  if (communication_model & 1) // AVOID_ALL_REDUCE
  {
    total_num_vertices_to_smooth = num_vtx_partition_boundary_local;
    total_num_vertices_to_recv = num_vtx_partition_boundary_remote;
    neighbourProcSend = new int[num_neighbourProc];
    neighbourProcRecv = new int[num_neighbourProc];
    for (i=0;i<num_neighbourProc;i++) {
	neighbourProcSend[i] = 0;
	neighbourProcRecv[i] = 0;
    }

    /* for each vertex we smooth find the processors we need to send it too */
    for (i=0;i<num_vtx_partition_boundary_local;i++) {
	list_size = vtx_off_proc_list_size[i];
	/* loop over its adjacent off-processor vertices */
	for (j=0;j<list_size; j++) {
	  /* get the processor id of these vertices */
	  incident_vtx = part_proc_owner[vtx_off_proc_list[i][j]];
	  /* check if we got this processor id before */
	  for (k=0;k<j;k++) {
	    if (incident_vtx == part_proc_owner[vtx_off_proc_list[i][k]]) {
	      /* if we have has this procesor id already we do not need to count it again */
	      incident_vtx = -1;
	      break;
	    }
	  }
	  /* if this was a new processor id */
	  if (incident_vtx != -1) {
	    /* find the processor in the list and increment its counter */
	    for (l = 0; l < num_neighbourProc; l++) {
	      if (neighbourProc[l] == incident_vtx) {
		neighbourProcSend[l]++;
		break;
	      }
	    }
	  }
	}
    }
    for (i=num_vtx_partition_boundary_local;i<num_vtx_partition_boundary;i++) {
	incident_vtx = part_proc_owner[i];
	for (l = 0; l < num_neighbourProc; l++) {
	  if (neighbourProc[l] == incident_vtx) {
	    neighbourProcRecv[l]++;
	    break;
	  }
	}
    }
    neighbourProcSendRemain = new int[num_neighbourProc];
    neighbourProcRecvRemain = new int[num_neighbourProc];
  }

  exportVtxGIDs = new int[num_vtx_partition_boundary];
  exportVtxLIDs = new int[num_vtx_partition_boundary];
  exportProc = new int[num_vtx_partition_boundary];
  in_independent_set = new bool[num_vtx_partition_boundary_local];

  return true;
}


void ParallelHelperImpl::compute_first_independent_set(msq_std::vector<Mesh::VertexHandle>& fixed_vertices)
{
  int i;

  // to avoid all reduce we need to know how many vertices we send & receive
  if (communication_model & 1) // AVOID_ALL_REDUCE
  {
    for (i=0;i<num_neighbourProc;i++) {
      neighbourProcSendRemain[i] = neighbourProcSend[i];
      neighbourProcRecvRemain[i] = neighbourProcRecv[i];
    }
  }

  // this is iteration zero of the bounday smooting process
  iteration = 0;

  // mark all boundary partition vertices as not smoothed
  for (i=0;i<num_vtx_partition_boundary;i++) {
    part_smoothed_flag[i] = 0;
  }

  // populates the in_independent_set and the vertex export arrays
  compute_independent_set();

  // counts how many vertices are already smoothed
  num_already_smoothed_vertices = 0;

  fixed_vertices.clear();

  /* mark which local boundary partition vertices are already smoothed (i.e. they are in the 1st independent set) */
  for (i=0;i<num_vtx_partition_boundary_local;i++) {
    if (in_independent_set[i]) {
      part_smoothed_flag[i] = 1;
      num_already_smoothed_vertices++;
    }
    else {
      fixed_vertices.push_back((*part_vertices)[i]); // fix vertices *not* in the independent set
    }
  }

  if (0) {printf("[%d]i%d after first we smoothed %d of %d\n",rank,iteration,num_already_smoothed_vertices,num_vtx_partition_boundary_local); fflush(NULL);}

  // fix the ghost vertices that are smoothed on another processor
  for (i=num_vtx_partition_boundary_local;i<num_vtx_partition_boundary;i++) {
    fixed_vertices.push_back((*part_vertices)[i]);
  }

  // fix the ghost vertices that are unused
  if (unused_ghost_vertices) {
    for (i=0;i<unused_ghost_vertices->size();i++) {
      fixed_vertices.push_back((*unused_ghost_vertices)[i]);
    }
  }
}

void ParallelHelperImpl::communicate_first_independent_set()
{
  switch (communication_model)
  {
  case TrulyNonBlocking:
    num_already_recv_vertices = comm_smoothed_vtx_tnb();
    break;
  case TrulyNonBlockingAvoidAllReduce:
    num_already_recv_vertices = comm_smoothed_vtx_tnb_no_all();
    break;
  case NonBlocking:
    num_already_recv_vertices = comm_smoothed_vtx_nb();
    break;
  case NonBlockingAvoidAllReduce:
    num_already_recv_vertices = comm_smoothed_vtx_nb_no_all();
    break;
  case Blocking:
    num_already_recv_vertices = comm_smoothed_vtx_b();
    break;
  case BlockingAvoidAllReduce:
    num_already_recv_vertices = comm_smoothed_vtx_b_no_all();
    break;
  }
  global_work_remains = (num_neighbourProc ? 1 : 0);
}

bool ParallelHelperImpl::compute_next_independent_set()
{
  if (global_work_remains && (iteration<20))
  {
    iteration++;
    compute_independent_set();
    next_vtx_partition_boundary = 0;
    return true;
  }
  else
  {
    return false;
  }
}

bool ParallelHelperImpl::get_next_partition_boundary_vertex(Mesquite::Mesh::VertexHandle& vertex_handle)
{
  while (next_vtx_partition_boundary < num_vtx_partition_boundary_local)
  {
    if (in_independent_set[next_vtx_partition_boundary]) {
      vertex_handle = (*part_vertices)[next_vtx_partition_boundary];
      num_already_smoothed_vertices++;
      assert(part_smoothed_flag[next_vtx_partition_boundary] == 0);
      part_smoothed_flag[next_vtx_partition_boundary] = 1;
      next_vtx_partition_boundary++;
      return true;
    }
    next_vtx_partition_boundary++;
  }
  if (0) {printf("[%d]i%d after next we smoothed %d of %d\n",rank,iteration,num_already_smoothed_vertices,num_vtx_partition_boundary_local); fflush(NULL);}
  return false;
}

void ParallelHelperImpl::communicate_next_independent_set()
{
  switch (communication_model)
  {
  case TrulyNonBlocking:
    num_already_recv_vertices += comm_smoothed_vtx_tnb();
    break;
  case TrulyNonBlockingAvoidAllReduce:
    num_already_recv_vertices += comm_smoothed_vtx_tnb_no_all();
    break;
  case NonBlocking:
    num_already_recv_vertices += comm_smoothed_vtx_nb();
    break;
  case NonBlockingAvoidAllReduce:
    num_already_recv_vertices += comm_smoothed_vtx_nb_no_all();
    break;
  case Blocking:
    num_already_recv_vertices += comm_smoothed_vtx_b();
    break;
  case BlockingAvoidAllReduce:
    num_already_recv_vertices += comm_smoothed_vtx_b_no_all();
    break;
  }

  if (communication_model & 1) // AVOID_ALL_REDUCE
  {
    global_work_remains = (total_num_vertices_to_smooth - num_already_smoothed_vertices) + (total_num_vertices_to_recv - num_already_recv_vertices); 
    if (0) {printf("[%d]i%d %d - %d + %d  - %d = %d \n",rank,iteration,total_num_vertices_to_smooth,num_already_smoothed_vertices,total_num_vertices_to_recv,num_already_recv_vertices,global_work_remains); fflush(NULL);}
  }
  else
  {
    int i, work_remains = 0;
    for (i=0;i<num_vtx_partition_boundary_local;i++) {
      if (part_smoothed_flag[i] == 0) {
	work_remains++;
      }
    }
    MPI_Allreduce(&work_remains, &global_work_remains, 1, MPI_INT, MPI_SUM, communicator);
  }
}

bool ParallelHelperImpl::smoothing_close()
{
  int i;

  if (vertices) delete vertices; vertices = 0;
  if (vtx_in_partition_boundary) delete [] vtx_in_partition_boundary; vtx_in_partition_boundary = 0;
  if (part_vertices) delete part_vertices; part_vertices = 0;
  if (unused_ghost_vertices) delete unused_ghost_vertices; unused_ghost_vertices = 0;
  if (part_proc_owner) delete [] part_proc_owner; part_proc_owner = 0;
  if (part_gid) delete [] part_gid; part_gid = 0;
  if (part_smoothed_flag) delete [] part_smoothed_flag; part_smoothed_flag = 0; 
  if (part_rand_number) delete [] part_rand_number; part_rand_number = 0;
  if (exportVtxGIDs) delete [] exportVtxGIDs; exportVtxGIDs = 0;
  if (exportVtxLIDs) delete [] exportVtxLIDs; exportVtxLIDs = 0;
  if (exportProc) delete [] exportProc; exportProc = 0;
  if (in_independent_set) delete [] in_independent_set; in_independent_set = 0;
  if (vid_hash) delete vid_hash; vid_hash = 0;
  if (neighbourProcSend) delete [] neighbourProcSend; neighbourProcSend = 0;
  if (neighbourProcRecv) delete [] neighbourProcRecv; neighbourProcRecv = 0;
  if (neighbourProcSendRemain) delete [] neighbourProcSendRemain; neighbourProcSendRemain = 0;
  if (neighbourProcRecvRemain) delete [] neighbourProcRecvRemain; neighbourProcRecvRemain = 0;
  if (vtx_off_proc_list_size) delete [] vtx_off_proc_list_size; vtx_off_proc_list_size = 0;
  for (i = 0; i < num_vtx_partition_boundary_local; i++) free(vtx_off_proc_list[i]);
  if (vtx_off_proc_list) delete [] vtx_off_proc_list; vtx_off_proc_list = 0;
  if (neighbourProc) free(neighbourProc); neighbourProc = 0;

  return true;
}

typedef struct VertexPack {
  double x;
  double y;
  double z;
  double glob_id;
} VertexPack;

int ParallelHelperImpl::comm_smoothed_vtx_tnb()
{
  int i,j,k;
  Mesquite::MsqError msq_err;

  // printf("[%d]i%d truly non blocking\n",rank, iteration);fflush(NULL);

  /* compute how many vertices we send to each processor */

  int numVtxPerProcSend[num_neighbourProc];
  int numVtxPerProcRecv[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProcSend[i] = 0;
    numVtxPerProcRecv[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	/* increment count */
	numVtxPerProcSend[j]++;
	/* end loop */
	break;
      }
    }
    /* did loop end without finding the processor */
    if (j == num_neighbourProc) {
      printf("[%d]i%d WARNING: did not find exportProc[%d] = %d in list of %d processors.\n",rank,iteration,i,exportProc[i],num_neighbourProc);
    }
  }

  /* tell each processor how many vertices they can expect from us */
  /* also ask each processor how many vertices we can expect from them */

  int num_neighbourProcSend = 0;
  int num_neighbourProcRecv = 0;
  MPI_Request requests_send[num_neighbourProc];
  MPI_Request requests_recv[num_neighbourProc];
  for (j = 0; j < num_neighbourProc; j++) {
    /* send the vertex count to this processor */
    if (0) {printf("[%d]i%d Announce send %d vertices from proc %d\n",rank,iteration,numVtxPerProcSend[j],neighbourProc[j]); fflush(NULL);}
    MPI_Isend(&(numVtxPerProcSend[j]),
	      1,
	      MPI_INT,
	      neighbourProc[j],
	      VERTEX_HEADER+iteration,
	      communicator,
	      &(requests_send[j]));
    num_neighbourProcSend++;

    /* recv the vertex count for this processor */
    if (0) {printf("[%d]i%d Listen  recv %d vertices from proc %d\n",rank,iteration,numVtxPerProcRecv[j],neighbourProc[j]); fflush(NULL);}
    MPI_Irecv(&(numVtxPerProcRecv[j]),
             1,
             MPI_INT,
             neighbourProc[j],
             VERTEX_HEADER+iteration,
             communicator,
	      &(requests_recv[j]));
    num_neighbourProcRecv++;
  }

  /* set up memory for the outgoing vertex data blocks */

  VertexPack* vertex_pack_export = (VertexPack*)malloc(sizeof(VertexPack)*(num_exportVtx+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices_export[num_neighbourProc];
  packed_vertices_export[0] = vertex_pack_export;
  for (i = 1; i < num_neighbourProc; i++) {
    packed_vertices_export[i] = packed_vertices_export[i-1] + numVtxPerProcSend[i-1];
  }

  /* place vertex data going to the same processor into consecutive memory space */

  int numVtxPerProcSendPACKED[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProcSendPACKED[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	VertexPack* packing_vertex = packed_vertices_export[j] + numVtxPerProcSendPACKED[j];
	numVtxPerProcSendPACKED[j]++;
	Mesquite::MsqVertex coordinates;
	mesh->vertices_get_coordinates(&(*part_vertices)[exportVtxLIDs[i]],&coordinates,1,msq_err);
	packing_vertex->x = coordinates[0];
	packing_vertex->y = coordinates[1];
	packing_vertex->z = coordinates[2];
	packing_vertex->glob_id = exportVtxGIDs[i];
      }
    }
  }

  /* wait until we have heard from all processors how many vertices we will receive from them */

  MPI_Status status[num_neighbourProc];

  if (num_neighbourProcRecv) {
    MPI_Waitall(num_neighbourProc, requests_recv, status);
  }

  /* how many vertices will we receive */

  int numVtxImport = 0;
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxImport += numVtxPerProcRecv[i];
  }

  /* set up memory for the incoming vertex data blocks */

  VertexPack* vertex_pack_import = (VertexPack*)malloc(sizeof(VertexPack)*(numVtxImport+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices_import[num_neighbourProc];
  packed_vertices_import[0] = vertex_pack_import;
  for (i = 1; i < num_neighbourProc; i++) {
    packed_vertices_import[i] = packed_vertices_import[i-1] + numVtxPerProcRecv[i-1];
  }

  /* post receives for all processors that have something for us */
  /* post sends for all processors that we have something for */

  num_neighbourProcRecv = 0;

  for (i = 0; i < num_neighbourProc; i++) {
    if (numVtxPerProcRecv[i]) {
      if (0) {printf("[%d]i%d Will recv %d vertices from proc %d\n",rank,iteration,numVtxPerProcRecv[i],neighbourProc[i]); fflush(NULL);}
      MPI_Irecv(packed_vertices_import[i],
		4*numVtxPerProcRecv[i],
		MPI_DOUBLE_PRECISION,
		neighbourProc[i],
		VERTEX_BLOCK+iteration,
		communicator,
		&(requests_recv[i]));
      num_neighbourProcRecv++;
    }
    else {
      requests_recv[i] = MPI_REQUEST_NULL;
    }
    if (numVtxPerProcSend[i]) {
      if (0) {printf("[%d]i%d Will send %d vertices to proc %d\n",rank,iteration,numVtxPerProcSend[i],neighbourProc[i]); fflush(NULL);}
      MPI_Isend(packed_vertices_export[i], 
		4*numVtxPerProcSend[i],
		MPI_DOUBLE_PRECISION,
		neighbourProc[i],
		VERTEX_BLOCK+iteration,
		communicator,
		&(requests_send[i]));
    }
    else {
      requests_send[i] = MPI_REQUEST_NULL;
    }
  }

  /* wait for some receive to arrive */ 
  
  int local_id;
  while (num_neighbourProcRecv) {
    MPI_Waitany(num_neighbourProc, requests_recv, &k, &status[0]);
    /* unpack all vertices */
    for (i = 0; i < numVtxPerProcRecv[k]; i++) {
      local_id = vertex_hash_find(vid_hash,(int)(packed_vertices_import[k][i].glob_id), neighbourProc[k]);
      if (local_id) {
	Mesquite::Vector3D coordinates;
	coordinates.set(packed_vertices_import[k][i].x, packed_vertices_import[k][i].y, packed_vertices_import[k][i].z);
	mesh->vertex_set_coordinates((*part_vertices)[local_id],coordinates,msq_err);
	assert(part_smoothed_flag[local_id] == 0);
	part_smoothed_flag[local_id] = 1;
      }
      else {
	printf("[%d]i%d communicate vertex with global_id %d not in mesh\n", rank,iteration,(int)(packed_vertices_import[k][i].glob_id));	  
      }
    }
    num_neighbourProcRecv--;
  }

  /* all receives have completed. it is save to release the memory */
  free(vertex_pack_import);

  /* wait until the sends have completed */

  if (num_neighbourProcSend) {
    MPI_Waitall(num_neighbourProc, requests_send, status);
  }

  /* all sends have completed. it is save to release the memory */
  free(vertex_pack_export);

  return numVtxImport;
}

int ParallelHelperImpl::comm_smoothed_vtx_tnb_no_all()
{
  int i,j,k;
  Mesquite::MsqError msq_err;

  // printf("[%d]i%d truly non blocking avoid reduce all\n", rank, iteration); fflush(NULL);

  /* compute how many vertices we send to each processor */

  int numVtxPerProcSend[num_neighbourProc];
  int numVtxPerProcRecv[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProcSend[i] = 0;
    numVtxPerProcRecv[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	/* increment count */
	numVtxPerProcSend[j]++;
	/* end loop */
	break;
      }
    }
    /* did loop end without finding the processor */
    if (j == num_neighbourProc) {
      printf("[%d]i%d WARNING: did not find exportProc[%d] = %d in list of %d processors.\n",rank,iteration,i,exportProc[i],num_neighbourProc);
    }
  }

  /* tell each processor how many vertices they can expect from us */
  /* also ask each processor how many vertices we can expect from them */

  int num_neighbourProcSend = 0;
  int num_neighbourProcRecv = 0;
  MPI_Request requests_send[num_neighbourProc];
  MPI_Request requests_recv[num_neighbourProc];
  for (j = 0; j < num_neighbourProc; j++) {
    if (neighbourProcSendRemain[j]) {
      /* send the vertex count to this processor */
      if (0) {printf("[%d]i%d Announce send %d vertices to proc %d\n",rank,iteration,numVtxPerProcSend[j],neighbourProc[j]); fflush(NULL);}
      MPI_Isend(&(numVtxPerProcSend[j]),
		1,
		MPI_INT,
		neighbourProc[j],
		VERTEX_HEADER+iteration,
		communicator,
		&(requests_send[j]));
      num_neighbourProcSend++;
    } else {
      requests_send[j] = MPI_REQUEST_NULL;
    }
    if (neighbourProcRecvRemain[j]) {
      /* recv the vertex count for this processor */
      if (0) {printf("[%d]i%d Listen recv xx vertices from proc %d\n",rank,iteration,neighbourProc[j]); fflush(NULL);}
      MPI_Irecv(&(numVtxPerProcRecv[j]),
		1,
		MPI_INT,
		neighbourProc[j],
		VERTEX_HEADER+iteration,
		communicator,
		&(requests_recv[j]));
      num_neighbourProcRecv++;
    } else {
      requests_recv[j] = MPI_REQUEST_NULL;
    }
  }

  /* set up memory for the outgoing vertex data blocks */

  VertexPack* vertex_pack_export = (VertexPack*)malloc(sizeof(VertexPack)*(num_exportVtx+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices_export[num_neighbourProc];
  packed_vertices_export[0] = vertex_pack_export;
  for (i = 1; i < num_neighbourProc; i++) {
    packed_vertices_export[i] = packed_vertices_export[i-1] + numVtxPerProcSend[i-1];
  }

  /* place vertex data going to the same processor into consecutive memory space */

  int numVtxPerProcSendPACKED[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProcSendPACKED[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	VertexPack* packing_vertex = packed_vertices_export[j] + numVtxPerProcSendPACKED[j];
	numVtxPerProcSendPACKED[j]++;
	Mesquite::MsqVertex coordinates;
	mesh->vertices_get_coordinates(&(*part_vertices)[exportVtxLIDs[i]],&coordinates,1,msq_err);
	packing_vertex->x = coordinates[0];
	packing_vertex->y = coordinates[1];
	packing_vertex->z = coordinates[2];
	packing_vertex->glob_id = exportVtxGIDs[i];
	if (0) printf("[%d]i%d vertex %d packed %g %g %g\n", rank,iteration,exportVtxGIDs[i],packing_vertex->x, packing_vertex->y, packing_vertex->z);
      }
    }
  }

  /* wait until we have heard from all processors how many vertices we will receive from them */

  MPI_Status status[num_neighbourProc];

  if (num_neighbourProcRecv) {
    MPI_Waitall(num_neighbourProc, requests_recv, status);
  }

  /* how many vertices will we receive */

  int numVtxImport = 0;
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxImport += numVtxPerProcRecv[i];
    neighbourProcRecvRemain[i] -= numVtxPerProcRecv[i];
    neighbourProcSendRemain[i] -= numVtxPerProcSend[i];
  }

  /* set up memory for the incoming vertex data blocks */

  VertexPack* vertex_pack_import = (VertexPack*)malloc(sizeof(VertexPack)*(numVtxImport+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices_import[num_neighbourProc];
  packed_vertices_import[0] = vertex_pack_import;
  for (i = 1; i < num_neighbourProc; i++) {
    packed_vertices_import[i] = packed_vertices_import[i-1] + numVtxPerProcRecv[i-1];
  }

  /* post receives for all processors that have something for us */
  /* post sends for all processors that we have something for */

  num_neighbourProcRecv = 0;

  for (i = 0; i < num_neighbourProc; i++) {
    if (0) {printf("[%d]i%d Will recv %d vertices from proc %d\n",rank,iteration,numVtxPerProcRecv[i],neighbourProc[i]); fflush(NULL);}
    if (numVtxPerProcRecv[i]) {
      MPI_Irecv(packed_vertices_import[i],
		4*numVtxPerProcRecv[i],
		MPI_DOUBLE_PRECISION,
		neighbourProc[i],
		VERTEX_BLOCK+iteration,
		communicator,
		&(requests_recv[i]));
      num_neighbourProcRecv++;
    }
    else {
      requests_recv[i] = MPI_REQUEST_NULL;
    }
    if (0) {printf("[%d]i%d Will send %d vertices to proc %d\n",rank,iteration,numVtxPerProcSend[i],neighbourProc[i]); fflush(NULL);}
    if (numVtxPerProcSend[i]) {
      MPI_Isend(packed_vertices_export[i], 
		4*numVtxPerProcSend[i],
		MPI_DOUBLE_PRECISION,
		neighbourProc[i],
		VERTEX_BLOCK+iteration,
		communicator,
		&(requests_send[i]));
    }
    else {
      requests_send[i] = MPI_REQUEST_NULL;
    }
  }

  /* wait for some receive to arrive */ 
  
  int local_id;
  while (num_neighbourProcRecv) {
    MPI_Waitany(num_neighbourProc, requests_recv, &k, &status[0]);
    /* unpack all vertices */
    for (i = 0; i < numVtxPerProcRecv[k]; i++) {
      local_id = vertex_hash_find(vid_hash,(int)(packed_vertices_import[k][i].glob_id), neighbourProc[k]);
      if (local_id) {
        Mesquite::Vector3D coordinates;
        coordinates.set(packed_vertices_import[k][i].x, packed_vertices_import[k][i].y, packed_vertices_import[k][i].z);
	if (0) printf("[%d]i%d vertex %d becomes %g %g %g\n", rank,iteration,local_id,packed_vertices_import[k][i].x, packed_vertices_import[k][i].y, packed_vertices_import[k][i].z);
        mesh->vertex_set_coordinates((*part_vertices)[local_id],coordinates,msq_err);
	assert(part_smoothed_flag[local_id] == 0);
        part_smoothed_flag[local_id] = 1;
      }
      else {
	printf("[%d]i%d communicate vertex from %d with global_id %d not in mesh\n", rank,iteration,neighbourProc[k],(int)(packed_vertices_import[k][i].glob_id));	  
      }
    }
    num_neighbourProcRecv--;
  }

  /* all receives have completed. it is save to release the memory */
  free(vertex_pack_import);

  /* wait until the sends have completed */

  if (num_neighbourProcSend) {
    MPI_Waitall(num_neighbourProc, requests_send, status);
  }

  /* all sends have completed. it is save to release the memory */
  free(vertex_pack_export);

  return numVtxImport;
}

int ParallelHelperImpl::comm_smoothed_vtx_nb()
{
  int i,j,k;
  Mesquite::MsqError msq_err;

  // printf("[%d] %d %d non blocking\n",rank, iteration, pass);fflush(NULL);

  /* how many vertices will we receive */

  int numVtxPerProcSend[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProcSend[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	/* increment count */
	numVtxPerProcSend[j]++;
	/* end loop */
	break;
      }
    }
    /* assert loop did not end without finding processor */
    assert(j != num_neighbourProc);
  }

  /* tell each processor how many vertices to expect */

  for (j = 0; j < num_neighbourProc; j++) {
    MPI_Send(&(numVtxPerProcSend[j]),
	     1,
	     MPI_INT,
	     neighbourProc[j],
	     VERTEX_HEADER+iteration,
	     communicator);

    //    printf("[%d]i%d Announcing %d vertices to proc %d\n",rank,iteration,numVtxPerProcSend[j],neighbourProc[j]); fflush(NULL);
  }

  /* place vertex data going to the same processor into consecutive memory space */
  VertexPack* vertex_pack_export = (VertexPack*)malloc(sizeof(VertexPack)*(num_exportVtx+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices_export[num_neighbourProc];
  packed_vertices_export[0] = vertex_pack_export;
  for (i = 1; i < num_neighbourProc; i++) {
    packed_vertices_export[i] = packed_vertices_export[i-1] + numVtxPerProcSend[i-1];
  }

  int* numVtxPerProcSendPACKED = new int[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProcSendPACKED[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	VertexPack* packing_vertex = packed_vertices_export[j] + numVtxPerProcSendPACKED[j];
	numVtxPerProcSendPACKED[j]++;
	Mesquite::MsqVertex coordinates;
	mesh->vertices_get_coordinates(&(*part_vertices)[exportVtxLIDs[i]],&coordinates,1,msq_err);
	packing_vertex->x = coordinates[0];
	packing_vertex->y = coordinates[1];
	packing_vertex->z = coordinates[2];
	packing_vertex->glob_id = exportVtxGIDs[i];
      }
    }
  }
  delete [] numVtxPerProcSendPACKED;

  /* now ask each processor how many vertices to expect */

  int num;
  int proc;
  int numVtxImport = 0;
  int num_neighbourProcRecv = 0;
  int numVtxPerProcRecv[num_neighbourProc];
  MPI_Status status;
  for (j = 0; j < num_neighbourProc; j++) {
    numVtxPerProcRecv[j] = 0;
  }

  for (j = 0; j < num_neighbourProc; j++) {
    /* get the vertex count for some processor */
    MPI_Recv(&num,              /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* of type int */
             MPI_ANY_SOURCE,    /* receive from any sender */
             VERTEX_HEADER+iteration,     /* receive only VERTEX HEADERs from this iteration */
             communicator,    /* default communicator */
             &status);          /* info about the received message */
    proc = status.MPI_SOURCE;
    /* will we import vertices from this processor */

    //    printf("[%d]i%d Heard we will receive %d vertices from proc %d\n",rank,iteration,num,proc); fflush(NULL);

    if (num) {
      /* increase number of processors we will receive from */
      num_neighbourProcRecv++;
      /* add number of vertices we will receive to the import total */
      numVtxImport += num;
    }
    /* find this processor in the list */
    for (i = 0; i < num_neighbourProc; i++) {
      if (neighbourProc[i] == proc) {
	numVtxPerProcRecv[i] = num;
	break;
      }
    }
    /* assert loop did not end without finding processor */
    assert(i != num_neighbourProc);
    if (0) printf("[%d]i%d Will receive %d vertices from proc %d\n",rank,iteration,num,proc); fflush(NULL);
  }

  /* create list of processors we receive from */
  int neighbourProcRecv[num_neighbourProcRecv];
  int numVtxPerProcRecvRecv[num_neighbourProcRecv];
  for (i = 0, k = 0; i < num_neighbourProc; i++) {
    if (0) printf("[%d]i%d Will receive %d vertices from proc %d\n",rank,iteration,numVtxPerProcRecv[i],neighbourProc[i]); fflush(NULL);
    if (numVtxPerProcRecv[i]) {
      neighbourProcRecv[k] = neighbourProc[i];
      numVtxPerProcRecvRecv[k] = numVtxPerProcRecv[i];
      k++;
    }
  }

  /* set up memory for the incoming vertex data blocks */
  VertexPack* vertex_pack_import = (VertexPack*)malloc(sizeof(VertexPack)*(numVtxImport+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices_import[num_neighbourProcRecv];
  packed_vertices_import[0] = vertex_pack_import;
  for (i = 1; i < num_neighbourProcRecv; i++) {
    packed_vertices_import[i] = packed_vertices_import[i-1] + numVtxPerProcRecvRecv[i-1];
  }

  /* receive from all processors that have something for us */
  MPI_Request request[num_neighbourProcRecv];
  for (j = 0; j < num_neighbourProcRecv; j++) {
    MPI_Irecv(packed_vertices_import[j],
	      4*numVtxPerProcRecvRecv[j],
	      MPI_DOUBLE_PRECISION,
	      neighbourProcRecv[j],
	      VERTEX_BLOCK+iteration,
	      communicator,
	      &(request[j]));
    if (0) {printf("[%d]i%d Scheduling receipt of %d vertices to proc %d\n",rank,iteration,numVtxPerProcRecvRecv[j],neighbourProcRecv[j]); fflush(NULL);}
  }

  /* now send the data blocks */

  MPI_Request requests_send[num_neighbourProc];
  for (j = 0; j < num_neighbourProc; j++) {
    if (numVtxPerProcSend[j]) {
      MPI_Isend(packed_vertices_export[j], 
		4*numVtxPerProcSend[j],
		MPI_DOUBLE_PRECISION,
		neighbourProc[j],
		VERTEX_BLOCK+iteration,
		communicator,
		&(requests_send[j]));
      if (0) {printf("[%d]i%d Scheduling send of %d vertices to proc %d\n",rank,iteration,numVtxPerProcSend[j],neighbourProc[j]); fflush(NULL);}
    } else {
      requests_send[j] = MPI_REQUEST_NULL;
    }
  }
  
  /* process messages as they arrive */ 

  int local_id;
  for (j = 0; j < num_neighbourProcRecv; j++) {
    MPI_Waitany(num_neighbourProcRecv, request, &k, &status);

    /* unpack messages */
    proc = status.MPI_SOURCE;
    int count;
    MPI_Get_count(&status, MPI_INT, &count);    
    if (0) printf("[%d]i%d Received %d (%d) vertices from proc %d (%d)\n",rank,iteration,numVtxPerProcRecvRecv[k],count,neighbourProcRecv[k],proc); fflush(NULL);
    for (i = 0; i < numVtxPerProcRecvRecv[k]; i++) {
      local_id = vertex_hash_find(vid_hash,(int)(packed_vertices_import[k][i].glob_id), neighbourProcRecv[k]);
      if (local_id) {
	Mesquite::Vector3D coordinates;
	coordinates.set(packed_vertices_import[k][i].x, packed_vertices_import[k][i].y, packed_vertices_import[k][i].z);
	mesh->vertex_set_coordinates((*part_vertices)[local_id],coordinates,msq_err);
	assert(part_smoothed_flag[local_id] == 0);
	part_smoothed_flag[local_id] = 1;
	if (0) printf("[%d]i%d updating vertex with global_id %d to %g %g %g \n", rank, iteration, (int)(packed_vertices_import[k][i].glob_id), packed_vertices_import[k][i].x, packed_vertices_import[k][i].y, packed_vertices_import[k][i].z);
      }
      else {
	printf("[%d]i%d communicate vertex with global_id %d not in mesh\n", rank,iteration,(int)(packed_vertices_import[k][i].glob_id));	  
      }
    }
  }

  /* all receives have completed. it is save to release the memory */
  free(vertex_pack_import);
  /* wait until the sends have completed */
  MPI_Status stati[num_neighbourProc];
  MPI_Waitall(num_neighbourProc, requests_send, stati);
  /* all sends have completed. it is save to release the memory */
  free(vertex_pack_export);

  return numVtxImport;
}

int ParallelHelperImpl::comm_smoothed_vtx_nb_no_all()
{
  int i,j,k;
  Mesquite::MsqError msq_err;

  // printf("[%d] %d %d non blocking avoid reduce all\n",rank, iteration, pass);fflush(NULL);

  /* how many vertices will we receive */

  int numVtxPerProcSend[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProcSend[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	/* increment count */
	numVtxPerProcSend[j]++;
	/* end loop */
	break;
      }
    }
    /* assert loop did not end without finding processor */
    assert(j != num_neighbourProc);
  }

  /* tell each processor how many vertices to expect */

  for (j = 0; j < num_neighbourProc; j++) {
    if (neighbourProcSendRemain[j]) {
      assert(neighbourProcSendRemain[j] >= numVtxPerProcSend[j]);
      neighbourProcSendRemain[j] -= numVtxPerProcSend[j];
    MPI_Send(&(numVtxPerProcSend[j]),
	     1,
	     MPI_INT,
	     neighbourProc[j],
	     VERTEX_HEADER+iteration,
	     communicator);

    //    printf("[%d]i%d Announcing %d vertices to proc %d\n",rank,iteration,numVtxPerProcSend[j],neighbourProc[j]); fflush(NULL);
    }
    else {
      assert(numVtxPerProcSend[j] == 0);
    }
  }

  /* place vertex data going to the same processor into consecutive memory space */
  VertexPack* vertex_pack_export = (VertexPack*)malloc(sizeof(VertexPack)*(num_exportVtx+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices_export[num_neighbourProc];
  packed_vertices_export[0] = vertex_pack_export;
  for (i = 1; i < num_neighbourProc; i++) {
    packed_vertices_export[i] = packed_vertices_export[i-1] + numVtxPerProcSend[i-1];
  }

  int* numVtxPerProcSendPACKED = new int[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProcSendPACKED[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	VertexPack* packing_vertex = packed_vertices_export[j] + numVtxPerProcSendPACKED[j];
	numVtxPerProcSendPACKED[j]++;
	Mesquite::MsqVertex coordinates;
	mesh->vertices_get_coordinates(&(*part_vertices)[exportVtxLIDs[i]],&coordinates,1,msq_err);
	packing_vertex->x = coordinates[0];
	packing_vertex->y = coordinates[1];
	packing_vertex->z = coordinates[2];
	packing_vertex->glob_id = exportVtxGIDs[i];
      }
    }
  }
  delete [] numVtxPerProcSendPACKED;

  /* now ask each processor how many vertices to expect */

  int num;
  int proc;
  int numVtxImport = 0;
  int num_neighbourProcRecv = 0;
  int numVtxPerProcRecv[num_neighbourProc];
  MPI_Status status;
  for (j = 0; j < num_neighbourProc; j++) {
    numVtxPerProcRecv[j] = 0;
  }

  int num_neighbourProcRecvRemain = 0;
  for (j = 0; j < num_neighbourProc; j++) {
    if (neighbourProcRecvRemain[j]) {
      num_neighbourProcRecvRemain++;
    }
  }
  for (j = 0; j < num_neighbourProcRecvRemain; j++) {
    /* get the vertex count for some processor */
    MPI_Recv(&num,              /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* of type int */
             MPI_ANY_SOURCE,    /* receive from any sender */
             VERTEX_HEADER+iteration,     /* receive only VERTEX HEADERs from this iteration */
             communicator,    /* default communicator */
             &status);          /* info about the received message */
    proc = status.MPI_SOURCE;
    /* will we import vertices from this processor */

    //    printf("[%d]i%d Heard we will receive %d vertices from proc %d\n",rank,iteration,num,proc); fflush(NULL);

    if (num) {
      /* increase number of processors we will receive from */
      num_neighbourProcRecv++;
      /* add number of vertices we will receive to the import total */
      numVtxImport += num;
    }
    /* find this processor in the list */
    for (i = 0; i < num_neighbourProc; i++) {
      if (neighbourProc[i] == proc) {
	numVtxPerProcRecv[i] = num;
	assert(neighbourProcRecvRemain[i] >= num);
	neighbourProcRecvRemain[i] -= num;
	break;
      }
    }
    /* assert loop did not end without finding processor */
    assert(i != num_neighbourProc);
    if (0) printf("[%d]i%d Will receive %d vertices from proc %d\n",rank,iteration,num,proc); fflush(NULL);
  }

  /* create list of processors we receive from */
  int neighbourProcRecv[num_neighbourProcRecv];
  int numVtxPerProcRecvRecv[num_neighbourProcRecv];
  for (i = 0, k = 0; i < num_neighbourProc; i++) {
    if (0) printf("[%d]i%d Will receive %d vertices from proc %d\n",rank,iteration,numVtxPerProcRecv[i],neighbourProc[i]); fflush(NULL);
    if (numVtxPerProcRecv[i]) {
      neighbourProcRecv[k] = neighbourProc[i];
      numVtxPerProcRecvRecv[k] = numVtxPerProcRecv[i];
      k++;
    }
  }

  /* set up memory for the incoming vertex data blocks */
  VertexPack* vertex_pack_import = (VertexPack*)malloc(sizeof(VertexPack)*(numVtxImport+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices_import[num_neighbourProcRecv];
  packed_vertices_import[0] = vertex_pack_import;
  for (i = 1; i < num_neighbourProcRecv; i++) {
    packed_vertices_import[i] = packed_vertices_import[i-1] + numVtxPerProcRecvRecv[i-1];
  }

  /* receive from all processors that have something for us */
  MPI_Request request[num_neighbourProcRecv];
  for (j = 0; j < num_neighbourProcRecv; j++) {
    MPI_Irecv(packed_vertices_import[j],
	      4*numVtxPerProcRecvRecv[j],
	      MPI_DOUBLE_PRECISION,
	      neighbourProcRecv[j],
	      VERTEX_BLOCK+iteration,
	      communicator,
	      &(request[j]));
    if (0) {printf("[%d]i%d Scheduling receipt of %d vertices to proc %d\n",rank,iteration,numVtxPerProcRecvRecv[j],neighbourProcRecv[j]); fflush(NULL);}
  }

  /* now send the data blocks */

  MPI_Request requests_send[num_neighbourProc];
  for (j = 0; j < num_neighbourProc; j++) {
    if (numVtxPerProcSend[j]) {
      MPI_Isend(packed_vertices_export[j], 
		4*numVtxPerProcSend[j],
		MPI_DOUBLE_PRECISION,
		neighbourProc[j],
		VERTEX_BLOCK+iteration,
		communicator,
		&(requests_send[j]));
      if (0) {printf("[%d]i%d Scheduling send of %d vertices to proc %d\n",rank,iteration,numVtxPerProcSend[j],neighbourProc[j]); fflush(NULL);}
    } else {
      requests_send[j] = MPI_REQUEST_NULL;
    }
  }
  
  /* process messages as they arrive */ 

  int local_id;
  for (j = 0; j < num_neighbourProcRecv; j++) {
    MPI_Waitany(num_neighbourProcRecv, request, &k, &status);

    /* unpack messages */
    proc = status.MPI_SOURCE;
    int count;
    MPI_Get_count(&status, MPI_INT, &count);    
    if (0) printf("[%d]i%d Received %d (%d) vertices from proc %d (%d)\n",rank,iteration,numVtxPerProcRecvRecv[k],count,neighbourProcRecv[k],proc); fflush(NULL);
    for (i = 0; i < numVtxPerProcRecvRecv[k]; i++) {
      local_id = vertex_hash_find(vid_hash,(int)(packed_vertices_import[k][i].glob_id), neighbourProcRecv[k]);
      if (local_id) {
	Mesquite::Vector3D coordinates;
	coordinates.set(packed_vertices_import[k][i].x, packed_vertices_import[k][i].y, packed_vertices_import[k][i].z);
	mesh->vertex_set_coordinates((*part_vertices)[local_id],coordinates,msq_err);
	assert(part_smoothed_flag[local_id] == 0);
	part_smoothed_flag[local_id] = 1;
	if (0) printf("[%d]i%d updating vertex with global_id %d to %g %g %g \n", rank, iteration, (int)(packed_vertices_import[k][i].glob_id), packed_vertices_import[k][i].x, packed_vertices_import[k][i].y, packed_vertices_import[k][i].z);
      }
      else {
	printf("[%d]i%d communicate vertex with global_id %d not in mesh\n", rank,iteration,(int)(packed_vertices_import[k][i].glob_id));	  
      }
    }
  }

  /* all receives have completed. it is save to release the memory */
  free(vertex_pack_import);
  /* wait until the sends have completed */
  MPI_Status stati[num_neighbourProc];
  MPI_Waitall(num_neighbourProc, requests_send, stati);
  /* all sends have completed. it is save to release the memory */
  free(vertex_pack_export);

  return numVtxImport;
}

int ParallelHelperImpl::comm_smoothed_vtx_b()
{
  int i,j;
  Mesquite::MsqError msq_err;

  // printf("[%d] %d %d blocking\n",rank, iteration, pass);fflush(NULL);

  /* how many vertices per processor */

  int numVtxPerProc[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProc[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	/* increment count */
	numVtxPerProc[j]++;
	/* end loop */
	break;
      }
      /* assert loop did not end without finding the processor */
      assert (j != num_neighbourProc);
    }
  }

  /* place vertices going to the same processor into consecutive memory space */

  VertexPack* vertex_pack = (VertexPack*)malloc(sizeof(VertexPack)*(num_exportVtx+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices[num_neighbourProc];
  VertexPack* packing_vertex;
  packed_vertices[0] = vertex_pack;
  for (i = 1; i < num_neighbourProc; i++) {
    packed_vertices[i] = packed_vertices[i-1] + numVtxPerProc[i-1];
  }

  int* numVtxPackedPerProc = new int[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPackedPerProc[i] = 0;
  }

  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	packing_vertex = packed_vertices[j] + numVtxPackedPerProc[j];
	numVtxPackedPerProc[j]++;
	Mesquite::MsqVertex coordinates;
	mesh->vertices_get_coordinates(&(*part_vertices)[exportVtxLIDs[i]],&coordinates,1,msq_err);
	packing_vertex->x = coordinates[0];
	packing_vertex->y = coordinates[1];
	packing_vertex->z = coordinates[2];
	packing_vertex->glob_id = exportVtxGIDs[i];
      }
    }
  }

  delete [] numVtxPackedPerProc;

  /* send each block so the corresponding processor preceeded by the number of vertices */
  
  for (j = 0; j < num_neighbourProc; j++) {

    //    printf("[%d]i%dp%d Announcing %d vertices to proc %d\n",rank,iteration,pass,numVtxPerProc[j],neighbourProc[j]); fflush(NULL);

    MPI_Send(&(numVtxPerProc[j]),
	     1,
	     MPI_INT,
	     neighbourProc[j],
	     VERTEX_HEADER+iteration,
	     communicator);

    // printf("[%d]i%dp%d Sending %d vertices to proc %d\n",rank,iteration,pass,numVtxPerProc[j],neighbourProc[j]); fflush(NULL);

    /* is there any vertex data to be sent */

    if (numVtxPerProc[j]) {
      MPI_Send(packed_vertices[j],
	       4*numVtxPerProc[j],
	       MPI_DOUBLE_PRECISION,
	       neighbourProc[j],
	       VERTEX_BLOCK+iteration,
	       communicator);

      // printf("[%d]i%dp%d Sent %d vertices to proc %d\n",rank,iteration,pass,numVtxPerProc[j],neighbourProc[j]); fflush(NULL);
    }
  }

  int num;
  int proc;
  int tag;
  int count;
  int numVtxImport = 0;
  MPI_Status status;
  int local_id;

  //  printf("[%d]i%dp%d Waiting to receive vertices from %d processors ... \n",rank,iteration,pass,num_neighbourProc); fflush(NULL);

  /* receiving blocks from other processors */

  for (j = 0; j < num_neighbourProc; j++) {
    MPI_Recv(&num,              /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* of type int */
             MPI_ANY_SOURCE,    /* receive from any sender */
             VERTEX_HEADER+iteration,     /* receive only VERTEX HEADERs */
             communicator,    /* default communicator */
             &status);          /* info about the received message */
    proc = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_INT, &count);

    //    printf("[%d]i%dp%d Receiving %d vertices from proc %d/%d/%d\n",rank,iteration,pass,num,proc,tag,count); fflush(NULL);

    /* is there any vertex data to be received */

    if (num) {

      numVtxImport += num;

      /* do we have enough space allocated */

      if (num_exportVtx + 10 < num) {
	if (vertex_pack) free(vertex_pack);
	num_exportVtx = num;
	vertex_pack = (VertexPack*)malloc(sizeof(VertexPack)*(num_exportVtx+10));
      }

      MPI_Recv(vertex_pack,          /* message buffer */
	       4*num,                /* num data's item with 4 doubles each */
	       MPI_DOUBLE_PRECISION, /* of type double */
	       proc,                 /* receive from this procesor only */
	       VERTEX_BLOCK+iteration,         /* receive only VERTEX BLOCKs */
	       communicator,       /* default communicator */
	       &status);             /* info about the received message */

      proc = status.MPI_SOURCE;
      tag = status.MPI_TAG;
      MPI_Get_count(&status, MPI_DOUBLE_PRECISION, &count);

      if (count != 4*num) printf("[%d]i%d WARNING: expected %d vertices = %d bytes from proc %d but only got %d bytes\n",rank,iteration,num,num*4,proc,count); fflush(NULL);

      //      printf("[%d]i%d Received %d vertices from proc %d/%d/%d\n",rank,iteration,num,proc,tag,count); fflush(NULL);

      /* update the received vertices in our boundary mesh */
      for (i = 0; i < num; i++) {
	/*	printf("[%d]i%d updating vertex %d with global_id %d\n",rank,iteration,i,(int)(vertex_pack[i].glob_id)); fflush(NULL); */
	local_id = vertex_hash_find(vid_hash,(int)(vertex_pack[i].glob_id), proc);
	if (local_id)
	{
	  Mesquite::Vector3D coordinates;
	  coordinates.set(vertex_pack[i].x, vertex_pack[i].y, vertex_pack[i].z);
	  mesh->vertex_set_coordinates((*part_vertices)[local_id],coordinates,msq_err);
	  assert(part_smoothed_flag[local_id] == 0);
	  part_smoothed_flag[local_id] = 1;
	  if (0) printf("[%d]i%d updating vertex with global_id %d to %g %g %g \n", rank,iteration, (int)(vertex_pack[i].glob_id), vertex_pack[i].x, vertex_pack[i].y, vertex_pack[i].z);
	}
	else {
	  printf("[%d]i%d communicate vertex with global_id %d not in mesh\n", rank,iteration, (int)(vertex_pack[i].glob_id));	  
	}
      }
    }
  }
  if (vertex_pack) free(vertex_pack);

  return numVtxImport;
}

int ParallelHelperImpl::comm_smoothed_vtx_b_no_all()
{
  int i,j;
  Mesquite::MsqError msq_err;

  // printf("[%d] %d %d blocking avoid reduce all\n",rank, iteration, pass);fflush(NULL);

  /* how many vertices per processor */

  int numVtxPerProc[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPerProc[i] = 0;
  }
  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	/* increment count */
	numVtxPerProc[j]++;
	/* end loop */
	break;
      }
      /* assert loop did not end without finding the processor */
      assert (j != num_neighbourProc);
    }
  }

  /* place vertices going to the same processor into consecutive memory space */

  VertexPack* vertex_pack = (VertexPack*)malloc(sizeof(VertexPack)*(num_exportVtx+10)); /* add 10 to have enough memory */
  VertexPack* packed_vertices[num_neighbourProc];
  VertexPack* packing_vertex;
  packed_vertices[0] = vertex_pack;
  for (i = 1; i < num_neighbourProc; i++) {
    packed_vertices[i] = packed_vertices[i-1] + numVtxPerProc[i-1];
  }

  int* numVtxPackedPerProc = new int[num_neighbourProc];
  for (i = 0; i < num_neighbourProc; i++) {
    numVtxPackedPerProc[i] = 0;
  }

  for (i = 0; i < num_exportVtx; i++) {
    for (j = 0; j < num_neighbourProc; j++) {
      if (exportProc[i] == neighbourProc[j]) {
	packing_vertex = packed_vertices[j] + numVtxPackedPerProc[j];
	numVtxPackedPerProc[j]++;
	Mesquite::MsqVertex coordinates;
	mesh->vertices_get_coordinates(&(*part_vertices)[exportVtxLIDs[i]],&coordinates,1,msq_err);
	packing_vertex->x = coordinates[0];
	packing_vertex->y = coordinates[1];
	packing_vertex->z = coordinates[2];
	packing_vertex->glob_id = exportVtxGIDs[i];
      }
    }
  }

  delete [] numVtxPackedPerProc;

  /* send each block so the corresponding processor preceeded by the number of vertices */
  
  for (j = 0; j < num_neighbourProc; j++) {

    if (neighbourProcSendRemain[j])
    {
      assert(neighbourProcSendRemain[j] >= numVtxPerProc[j]);
      neighbourProcSendRemain[j] -= numVtxPerProc[j];
      
      // printf("[%d]i%dp%d Announcing %d vertices to proc %d\n",rank,iteration,pass,numVtxPerProc[j],neighbourProc[j]); fflush(NULL);

      MPI_Send(&(numVtxPerProc[j]),
	       1,
	       MPI_INT,
	       neighbourProc[j],
	       VERTEX_HEADER+iteration,
	       communicator);
      
      // printf("[%d]i%dp%d Sending %d vertices to proc %d\n",rank,iteration,pass,numVtxPerProc[j],neighbourProc[j]); fflush(NULL);

      /* is there any vertex data to be sent */
      
      if (numVtxPerProc[j]) {
	MPI_Send(packed_vertices[j],
		 4*numVtxPerProc[j],
		 MPI_DOUBLE_PRECISION,
		 neighbourProc[j],
		 VERTEX_BLOCK+iteration,
		 communicator);
	
	// printf("[%d]i%dp%d Sent %d vertices to proc %d\n",rank,iteration,pass,numVtxPerProc[j],neighbourProc[j]); fflush(NULL);
      }
    }
    else
    {
      // printf("[%d]i%dp%d no longer sending to %d\n",rank,iteration,neighbourProc[j]); fflush(NULL);
    }
  }

  int num;
  int proc;
  int tag;
  int count;
  int numVtxImport = 0;
  MPI_Status status;
  int local_id;

  //  printf("[%d]i%dp%d Waiting to receive vertices from %d processors ... \n",rank,iteration,pass,num_neighbourProc); fflush(NULL);

  /* receiving blocks in any order from other processors */

  int num_neighbourProcRecvRemain = 0;
  for (j = 0; j < num_neighbourProc; j++) {
    if (neighbourProcRecvRemain[j]) {
      num_neighbourProcRecvRemain++;
    }
  }
  for (j = 0; j < num_neighbourProcRecvRemain; j++) {
    MPI_Recv(&num,              /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* of type int */
             MPI_ANY_SOURCE,    /* receive from any sender */
             VERTEX_HEADER+iteration,     /* receive only VERTEX HEADERs */
             communicator,    /* default communicator */
             &status);          /* info about the received message */
    proc = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_INT, &count);

    // printf("[%d]i%dp%d Receiving %d vertices from proc %d/%d/%d\n",rank,iteration,pass,num,proc,tag,count); fflush(NULL);

    /* is there any vertex data to be received */

    if (num) {

      numVtxImport += num;

      /* do we have enough space allocated */

      if (num_exportVtx + 10 < num) {
	if (vertex_pack) free(vertex_pack);
	num_exportVtx = num;
	vertex_pack = (VertexPack*)malloc(sizeof(VertexPack)*(num_exportVtx+10));
      }

      MPI_Recv(vertex_pack,          /* message buffer */
	       4*num,                /* num data's item with 4 doubles each */
	       MPI_DOUBLE_PRECISION, /* of type double */
	       proc,                 /* receive from this procesor only */
	       VERTEX_BLOCK+iteration,         /* receive only VERTEX BLOCKs */
	       communicator,       /* default communicator */
	       &status);             /* info about the received message */

      proc = status.MPI_SOURCE;
      tag = status.MPI_TAG;
      MPI_Get_count(&status, MPI_DOUBLE_PRECISION, &count);

      if (count != 4*num) printf("[%d]i%d WARNING: expected %d vertices = %d bytes from proc %d but only got %d bytes\n",rank,iteration,num,num*4,proc,count); fflush(NULL);

      // printf("[%d]i%d Received %d vertices from proc %d/%d/%d\n",rank,iteration,num,proc,tag,count); fflush(NULL);

      /* update the received vertices in our boundary mesh */
      for (i = 0; i < num; i++) {
	/*	printf("[%d]i%d updating vertex %d with global_id %d\n",rank,iteration,i,(int)(vertex_pack[i].glob_id)); fflush(NULL); */
	local_id = vertex_hash_find(vid_hash,(int)(vertex_pack[i].glob_id), proc);
	if (local_id)
	{
	  Mesquite::Vector3D coordinates;
	  coordinates.set(vertex_pack[i].x, vertex_pack[i].y, vertex_pack[i].z);
	  mesh->vertex_set_coordinates((*part_vertices)[local_id],coordinates,msq_err);
	  assert(part_smoothed_flag[local_id] == 0);
	  part_smoothed_flag[local_id] = 1;
	  if (0 && rank == 1) printf("[%d]i%d updating vertex with global_id %d to %g %g %g \n", rank,iteration, (int)(vertex_pack[i].glob_id), vertex_pack[i].x, vertex_pack[i].y, vertex_pack[i].z);
	}
	else {
	  printf("[%d]i%d communicate vertex with global_id %d not in mesh\n", rank,iteration, (int)(vertex_pack[i].glob_id));	  
	}
      }
    }
    for (i = 0; i < num_neighbourProc; i++) {
      if (proc == neighbourProc[i]) {
	assert(neighbourProcRecvRemain[i] >= num);
	neighbourProcRecvRemain[i] -= num;
	break;
      }
    }
  }
  if (vertex_pack) free(vertex_pack);

  return numVtxImport;
}

/***********************************************************************
COMPUTING THE INDEPENDENT SET AND EXPORT INFORMATION
***********************************************************************/
/* the only thing that can prevent a vertex from being in the independent
   set is that is has an neighbouring vertex on a different processor that
   has not been smoothed yet and that has a larger random number */

void ParallelHelperImpl::compute_independent_set()
{
  int i,j,k,l;
  int incident_vtx;
  bool done;

  for (i=0;i<num_vtx_partition_boundary_local;i++) in_independent_set[i] = false;

  bool found_more = true;
  int found_iter = 0;

  num_exportVtx = 0;

  while (found_more)
  {  
    found_iter++;
    found_more = false;
    for (i=0;i<num_vtx_partition_boundary_local;i++) {
      /* if this vertex could become part of the independent set */
      if (part_smoothed_flag[i] == 0 && in_independent_set[i] == false) {
	/* assume it's in the independent set */
	done = false;
	/* then loop over the neighbors it has on other processors */
	for (j = 0; !done && j < vtx_off_proc_list_size[i]; j++) {
	  incident_vtx = vtx_off_proc_list[i][j];
	  /* if this neighbour has not yet been smoothed and is not covered and has
	     a higher rand_number than me, then I am not in the independent set */
	  if ( ( part_smoothed_flag[incident_vtx] == 0 ) && ( part_rand_number[i] < part_rand_number[incident_vtx]) ) {
	    done = true;
	  }
	}
	/* if the vertex is in the independent set, add it to the export list */
	if (!done) {
	  found_more = true;
	  //	if (found_iter > 1) printf("[%d]i%d found another one %d in iteration %d\n",rank,iteration,i, found_iter);
	  in_independent_set[i] = true;
	  /* mark vertices with lower random numbers as covered */
	  for (j = 0; j < vtx_off_proc_list_size[i]; j++) {
	    incident_vtx = vtx_off_proc_list[i][j];
	    /* if this neighbour has not yet been smoothed or covered mark it as covered */
	    if ( part_smoothed_flag[incident_vtx] == 0 ) {
	      part_smoothed_flag[incident_vtx] = 2;
	    }
	  }
	  k = num_exportVtx;
	  /* then loop over the neighbors it has on other processors */
	  for (j = 0; j < vtx_off_proc_list_size[i]; j++) {
	    incident_vtx = vtx_off_proc_list[i][j];
	    /* check to see if this processor already on the list */
	    done = false;
	    for (l=k; l < num_exportVtx && !done; l++) {
	      if  (exportProc[l] == part_proc_owner[incident_vtx]) {
		done = true;
	      }
	    }
	    /* if it's not on the list add it */
	    if (!done) {
	      exportVtxLIDs[num_exportVtx] = i;
	      exportVtxGIDs[num_exportVtx] = part_gid[i];
	      exportProc[num_exportVtx] = part_proc_owner[incident_vtx];
	      num_exportVtx++;
	    }
	  }
	}
      }
    }
  }
  
  if (0)
  {
    int in_set  = 0;
    for (i=0;i<num_vtx_partition_boundary_local;i++) if (in_independent_set[i]) in_set++;;
    printf("[%d]i%d independent set has %d of %d vertices sent out %d times\n",rank,iteration, in_set, num_vtx_partition_boundary_local, num_exportVtx);
    fflush(NULL);
  }

  /* unmark the vertices that have been marked as covered */
  for (i=num_vtx_partition_boundary_local; i<num_vtx_partition_boundary; i++) if (part_smoothed_flag[i] == 2) part_smoothed_flag[i] = 0;
}

int ParallelHelperImpl::get_rank() const {
  return rank;
}

bool ParallelHelperImpl::is_our_element(Mesquite::Mesh::ElementHandle element_handle) const {
  int i;
  MsqError err;
  msq_std::vector<Mesh::VertexHandle> vertices;
  msq_std::vector<size_t> junk;
  mesh->elements_get_attached_vertices(&element_handle, 1, vertices, junk, err);
  int num_verts = vertices.size();
  int* proc_ids = new int[num_verts];
  mesh->vertices_get_processor_id(&vertices[0], &proc_ids[0], num_verts, err);
  int max_proc_id = proc_ids[0];
  for (i = 1; i < num_verts; i++)
    if (max_proc_id < proc_ids[i]) max_proc_id = proc_ids[i];
  delete [] proc_ids;
  return (max_proc_id == rank);
}

bool ParallelHelperImpl::is_our_vertex(Mesquite::Mesh::VertexHandle vertex_handle) const {
  int proc_id;
  MsqError err;
  mesh->vertices_get_processor_id(&vertex_handle, &proc_id, 1, err);
  return (proc_id == rank);
}

void ParallelHelperImpl::communicate_min_max_to_all(double* minimum, double* maximum) const {
  double d_min[2];
  double d_min_recv[2];
  d_min[0] = -(*maximum);
  d_min[1] = *minimum;
  MPI_Allreduce(d_min, d_min_recv, 2, MPI_DOUBLE, MPI_MIN, communicator);
  *maximum = -d_min_recv[0];
  *minimum =  d_min_recv[1];
}

void ParallelHelperImpl::communicate_min_max_to_zero(double* minimum, double* maximum) const {
  double d_min[2];
  double d_min_recv[2];
  d_min[0] = -(*maximum);
  d_min[1] = *minimum;
  MPI_Reduce(d_min, d_min_recv, 2, MPI_DOUBLE, MPI_MIN, 0, communicator);
  if (rank == 0) {
    *maximum = -d_min_recv[0];
    *minimum =  d_min_recv[1];
  }
}

void ParallelHelperImpl::communicate_sums_to_zero(int* inverted, int* indeterminate, long unsigned int* count, long unsigned int* invalid, double* sum, double *sqrSum) const {
  double d_sum[6];
  double d_sum_recv[6];

  d_sum[0] = (double)(*inverted);
  d_sum[1] = (double)(*indeterminate);
  d_sum[2] = (double)(*count);
  d_sum[3] = (double)(*invalid);
  d_sum[4] = *sum;
  d_sum[5] = *sqrSum;

  MPI_Reduce(d_sum, d_sum_recv, 6, MPI_DOUBLE, MPI_SUM, 0, communicator);

  if (rank == 0) {
    *inverted = (int)d_sum_recv[0];
    *indeterminate = (int)d_sum_recv[1];
    *count = (long unsigned int)d_sum_recv[2];
    *invalid = (long unsigned int)d_sum_recv[3];
    *sum = d_sum_recv[4];
    *sqrSum = d_sum_recv[5];
  }
}

void ParallelHelperImpl::communicate_power_sum_to_zero(double* pMean) const
{
  double result;
  MPI_Reduce( pMean, &result, 1, MPI_DOUBLE, MPI_SUM, 0, communicator );
  if (rank == 0)
    *pMean = result;
}

void ParallelHelperImpl::communicate_histogram_to_zero(msq_std::vector<int> &histogram) const {
  msq_std::vector<int> histogram_recv(histogram.size());
  MPI_Reduce(&(histogram[0]), &(histogram_recv[0]), histogram.size(), MPI_INT, MPI_SUM, 0, communicator);
  if (rank == 0) {
    histogram.swap( histogram_recv );
  }
}

} // namespace mesquite
