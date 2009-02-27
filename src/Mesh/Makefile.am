SRC_MESH_DIR = src/Mesh
MSQ_INCLUDES += -I$(top_srcdir)/$(SRC_MESH_DIR)
INCLUDED_MAKEFILES += $(SRC_MESH_DIR)/Makefile.am

if WITH_MPI
  MPI_SRCS = $(SRC_MESH_DIR)/ParallelHelper.cpp
  MPI_INCS = $(SRC_MESH_DIR)/ParallelHelper.hpp
else
  MPI_SRCS = 
  MPI_INCS = 
endif

MSQ_SRCS += \
	$(SRC_MESH_DIR)/ArrayMesh.cpp \
	$(SRC_MESH_DIR)/ElementPatches.cpp \
	$(SRC_MESH_DIR)/ExtraData.cpp \
	$(SRC_MESH_DIR)/GlobalPatch.cpp \
	$(SRC_MESH_DIR)/IdealElements.cpp \
	$(SRC_MESH_DIR)/MeshDecorator.cpp \
	$(SRC_MESH_DIR)/MeshImpl.cpp \
	$(SRC_MESH_DIR)/MeshImplData.cpp \
	$(SRC_MESH_DIR)/MeshImplTags.cpp \
	$(SRC_MESH_DIR)/MeshWriter.cpp \
	$(SRC_MESH_DIR)/MsqIBase.cpp \
	$(SRC_MESH_DIR)/MsqIBase.hpp \
	$(SRC_MESH_DIR)/MsqIGeom.cpp \
	$(SRC_MESH_DIR)/MsqIMesh.cpp \
	$(SRC_MESH_DIR)/MsqMeshEntity.cpp \
	$(SRC_MESH_DIR)/MsqVertex.cpp \
	$(SRC_MESH_DIR)/ParallelMeshImpl.cpp \
	$(SRC_MESH_DIR)/PatchData.cpp \
	$(SRC_MESH_DIR)/PatchIterator.cpp \
	$(SRC_MESH_DIR)/TagVertexMesh.cpp \
	$(SRC_MESH_DIR)/TopologyInfo.cpp \
	$(SRC_MESH_DIR)/VertexPatches.cpp \
        $(MPI_SRCS)

MSQ_HDRS += \
	$(SRC_MESH_DIR)/ArrayMesh.hpp \
	$(SRC_MESH_DIR)/ElementPatches.hpp \
	$(SRC_MESH_DIR)/ExtraData.hpp \
	$(SRC_MESH_DIR)/ExtraDataUser.hpp \
	$(SRC_MESH_DIR)/GlobalPatch.hpp \
	$(SRC_MESH_DIR)/IdealElements.hpp \
	$(SRC_MESH_DIR)/MeshDecorator.hpp \
	$(SRC_MESH_DIR)/MeshImpl.hpp \
	$(SRC_MESH_DIR)/MeshImplData.hpp \
	$(SRC_MESH_DIR)/MeshImplTags.hpp \
	$(SRC_MESH_DIR)/MeshWriter.hpp \
	$(SRC_MESH_DIR)/MsqFreeVertexIndexIterator.hpp \
	$(SRC_MESH_DIR)/MsqIGeom.hpp \
	$(SRC_MESH_DIR)/MsqIMesh.hpp \
	$(SRC_MESH_DIR)/MsqMeshEntity.hpp \
	$(SRC_MESH_DIR)/MsqVertex.hpp \
	$(SRC_MESH_DIR)/ParallelMeshImpl.hpp \
	$(SRC_MESH_DIR)/PatchData.hpp \
	$(SRC_MESH_DIR)/PatchIterator.hpp \
	$(SRC_MESH_DIR)/PatchSet.hpp \
	$(SRC_MESH_DIR)/TagVertexMesh.hpp \
	$(SRC_MESH_DIR)/TopologyInfo.hpp \
	$(SRC_MESH_DIR)/VertexPatches.hpp \
	$(MPI_INCS)
