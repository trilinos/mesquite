/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file CachingTargetTest.cpp
 *  \brief Unit tests for CachingTargetCalculator class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "CachingTargetCalculator.hpp"
#include "PatchData.hpp"
#include "UnitUtil.hpp"
#include "PatchDataInstances.hpp"
#include "SamplePoints.hpp"
#include <cppunit/extensions/HelperMacros.h>
#include <algorithm>

using namespace Mesquite;

class CachedTargetCalculator;

class CachingTargetTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(CachingTargetTest);
  CPPUNIT_TEST (test_3d_targets_cached);
  CPPUNIT_TEST (test_2d_targets_cached);
  CPPUNIT_TEST (test_3d_surf_targets_cached);
  CPPUNIT_TEST (test_3d_target_values);
  CPPUNIT_TEST (test_2d_target_values);
  CPPUNIT_TEST (test_3d_surf_target_values);
  CPPUNIT_TEST (test_3d_target_subpatch);
  CPPUNIT_TEST (test_2d_target_subpatch);
  CPPUNIT_TEST (test_3d_surf_target_subpatch);
  CPPUNIT_TEST (test_cache_cleared);
  CPPUNIT_TEST_SUITE_END();
  
  PatchData patch_3d, patch_2d;
  CachedTargetCalculator* cached, *cached_3d_surf;
  CachingTargetCalculator* cacher, *cacher_3d_surf;
  SamplePoints* samples;
  
  unsigned request_all_targets_3d();
  unsigned request_all_targets_2d();
  unsigned request_all_targets_3d_surf();

public:
  
  void setUp();
  void tearDown();

  void test_3d_targets_cached();
  void test_2d_targets_cached();
  void test_3d_surf_targets_cached();
  void test_3d_target_values();
  void test_2d_target_values();
  void test_3d_surf_target_values();
  void test_3d_target_subpatch();
  void test_2d_target_subpatch();
  void test_3d_surf_target_subpatch();
  void test_cache_cleared();
  
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(CachingTargetTest, "CachingTargetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(CachingTargetTest, "Unit");

class CachedTargetCalculator : public TargetCalculator
{
  private:
    unsigned called_3d, called_2d;
    bool surfaces_3d;
  
  public:
    
    static MsqMatrix<3,3> make_3d( size_t elem, Sample sample );
    static MsqMatrix<3,2> make_2d( size_t elem, Sample sample );
    
    CachedTargetCalculator( bool surface_elems_have_3D_targets ) 
      : called_3d(0), called_2d(0), surfaces_3d(surface_elems_have_3D_targets) {}
    
    virtual bool get_3D_target( PatchData&, size_t elem, const SamplePoints*, Sample sample, MsqMatrix<3,3>& result, MsqError& )
      { ++called_3d; result = make_3d( elem, sample); return true; }
    
    virtual bool get_2D_target( PatchData&, size_t elem, const SamplePoints*, Sample sample, MsqMatrix<3,2>& result, MsqError& )
      { ++called_2d; result = make_2d( elem, sample); return true; }
    
    void clear() 
      { called_3d = 0; called_2d = 0; }
      
    unsigned calls_3d() const 
      { return called_3d; }
      
    unsigned calls_2d() const 
      { return called_2d; }
      
    bool surface_targets_are_3D() const
      { return surfaces_3d; }
};

MsqMatrix<3,3> CachedTargetCalculator::make_3d( size_t elem, Sample sample )
{
  double v = 100. * elem + 4 * sample.number + sample.dimension + 1;
  const double values[] = { v, 0, 0,
                            0, v, 0,
                            0, 0, v };
  return MsqMatrix<3,3>(values);
}

MsqMatrix<3,2> CachedTargetCalculator::make_2d( size_t elem, Sample sample )
{
  double v = 100. * elem + 4 * sample.number + sample.dimension + 1;
  const double values[] = { v, 0,
                            0, v,
                            0.5*v, 0.5*v };
  return MsqMatrix<3,2>(values);
}

void CachingTargetTest::setUp()
{
    // make sure these are null so that if we fail within setUp,
    // tearDown doesn't try to delete stale pointers
  samples = 0;
  cached = 0;
  cacher = 0;
  
  MsqError err;
  create_four_quads_patch( patch_2d, err ); CPPUNIT_ASSERT(!err);
  create_qm_two_hex_patch( patch_3d, err ); CPPUNIT_ASSERT(!err);
  
  samples = new SamplePoints( true, false, false, false );
  cached = new CachedTargetCalculator( false );
  cached_3d_surf = new CachedTargetCalculator( true );
  cacher = new CachingTargetCalculator( cached, samples );
  cacher_3d_surf = new CachingTargetCalculator( cached_3d_surf, samples );
}

static void get_samples( EntityTopology type,
                         SamplePoints* pts,
                         msq_std::vector<Sample>& samples )
{
  samples.clear();
  for (unsigned d = 0; d < TopologyInfo::dimension(type); ++d) {
    if (!pts->will_sample_at( type, d ))
      continue;
    unsigned n = TopologyInfo::adjacent( type, d );
    for (unsigned i = 0; i < n; ++i)
      samples.push_back( Sample(d, i) );
  }
  if (pts->will_sample_at( type, TopologyInfo::dimension(type) ))
    samples.push_back( Sample( TopologyInfo::dimension(type), 0 ) );
}

unsigned CachingTargetTest::request_all_targets_3d()
{
  unsigned total = 0;
  MsqMatrix<3,3> W;
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
  for (size_t i = 0; i < patch_3d.num_elements(); ++i)
  {
    MsqMeshEntity& elem = patch_3d.element_by_index(i);
    get_samples( elem.get_element_type(), samples, locations );
    total += locations.size();
    for (unsigned j = 0; j < locations.size(); ++j)
    {
      bool rval = cacher->get_3D_target( patch_3d, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval);
      CPPUNIT_ASSERT(!err);
    }
  }
  return total;
}

unsigned CachingTargetTest::request_all_targets_2d()
{
  unsigned total = 0;
  MsqMatrix<3,2> W;
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
  for (size_t i = 0; i < patch_2d.num_elements(); ++i)
  {
    MsqMeshEntity& elem = patch_2d.element_by_index(i);
    get_samples( elem.get_element_type(), samples, locations );
    total += locations.size();
    for (unsigned j = 0; j < locations.size(); ++j)
    {
      bool rval = cacher->get_2D_target( patch_2d, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval);
      CPPUNIT_ASSERT(!err);
    }
  }
  return total;
}

unsigned CachingTargetTest::request_all_targets_3d_surf()
{
  unsigned total = 0;
  MsqMatrix<3,3> W;
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
  for (size_t i = 0; i < patch_2d.num_elements(); ++i)
  {
    MsqMeshEntity& elem = patch_2d.element_by_index(i);
    get_samples( elem.get_element_type(), samples, locations );
    total += locations.size();
    for (unsigned j = 0; j < locations.size(); ++j)
    {
      bool rval = cacher_3d_surf->get_3D_target( patch_2d, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval);
      CPPUNIT_ASSERT(!err);
    }
  }
  return total;
}

void CachingTargetTest::tearDown()
{
  delete cacher;
  delete cacher_3d_surf;
  delete cached;
  delete cached_3d_surf;
  delete samples;
}

void CachingTargetTest::test_3d_targets_cached()
{
  CPPUNIT_ASSERT_EQUAL( cached->calls_3d(), 0u );
  CPPUNIT_ASSERT_EQUAL( cached->calls_2d(), 0u );
  unsigned count = request_all_targets_3d();
  CPPUNIT_ASSERT_EQUAL( cached->calls_3d(), count );
  CPPUNIT_ASSERT_EQUAL( cached->calls_2d(), 0u );
  cached->clear();
  request_all_targets_3d();
  CPPUNIT_ASSERT_EQUAL( cached->calls_3d(), 0u );
}

void CachingTargetTest::test_2d_targets_cached()
{
  CPPUNIT_ASSERT_EQUAL( cached->calls_2d(), 0u );
  CPPUNIT_ASSERT_EQUAL( cached->calls_3d(), 0u );
  unsigned count = request_all_targets_2d();
  CPPUNIT_ASSERT_EQUAL( cached->calls_2d(), count );
  CPPUNIT_ASSERT_EQUAL( cached->calls_3d(), 0u );
  cached->clear();
  request_all_targets_2d();
  CPPUNIT_ASSERT_EQUAL( cached->calls_2d(), 0u );
}

void CachingTargetTest::test_3d_surf_targets_cached()
{
  CPPUNIT_ASSERT_EQUAL( cached_3d_surf->calls_3d(), 0u );
  CPPUNIT_ASSERT_EQUAL( cached_3d_surf->calls_2d(), 0u );
  unsigned count = request_all_targets_3d_surf();
  CPPUNIT_ASSERT_EQUAL( cached_3d_surf->calls_3d(), count );
  CPPUNIT_ASSERT_EQUAL( cached_3d_surf->calls_2d(), 0u );
  cached_3d_surf->clear();
  request_all_targets_3d_surf();
  CPPUNIT_ASSERT_EQUAL( cached_3d_surf->calls_3d(), 0u );
}

void CachingTargetTest::test_3d_target_values()
{
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
    // evaluate all once to make sure we test the cached values
  request_all_targets_3d();
  
    // test each value
  for (size_t i = 0; i < patch_3d.num_elements(); ++i)
  {
    MsqMeshEntity& elem = patch_3d.element_by_index(i);

    get_samples( elem.get_element_type(), samples, locations );
    for (unsigned j = 0; j < locations.size(); ++j)
    {
      MsqMatrix<3,3> W;
      bool rval = cacher->get_3D_target( patch_3d, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval && !err);
      
      MsqMatrix<3,3> M = CachedTargetCalculator::make_3d( i, locations[j] );
      ASSERT_MATRICES_EQUAL( W, M, DBL_EPSILON );
    }
  }
}

void CachingTargetTest::test_2d_target_values()
{
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
    // evaluate all once to make sure we test the cached values
  request_all_targets_2d();
  
    // test each value
  for (size_t i = 0; i < patch_2d.num_elements(); ++i)
  {
    MsqMeshEntity& elem = patch_2d.element_by_index(i);

    get_samples( elem.get_element_type(), samples, locations );
    for (unsigned j = 0; j < locations.size(); ++j)
    {
      MsqMatrix<3,2> W;
      bool rval = cacher->get_2D_target( patch_2d, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval && !err);
      
      MsqMatrix<3,2> M = CachedTargetCalculator::make_2d( i, locations[j] );
      ASSERT_MATRICES_EQUAL( W, M, DBL_EPSILON );
    }
  }
}

void CachingTargetTest::test_3d_surf_target_values()
{
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
    // evaluate all once to make sure we test the cached values
  request_all_targets_3d_surf();
  
    // test each value
  for (size_t i = 0; i < patch_2d.num_elements(); ++i)
  {
    MsqMeshEntity& elem = patch_2d.element_by_index(i);

    get_samples( elem.get_element_type(), samples, locations );
    for (unsigned j = 0; j < locations.size(); ++j)
    {
      MsqMatrix<3,3> W;
      bool rval = cacher_3d_surf->get_3D_target( patch_2d, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval && !err);
      
      MsqMatrix<3,3> M = CachedTargetCalculator::make_3d( i, locations[j] );
      ASSERT_MATRICES_EQUAL( W, M, DBL_EPSILON );
    }
  }
}

void CachingTargetTest::test_3d_target_subpatch()
{
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
    // cache some values on the main patch
  request_all_targets_3d();
  
    // clear the count so we know if any additional 
    // evalutions of the base target calculator are
    // done during subpatch creation.
  cached->clear();
  
    // create a sub-patch
  CPPUNIT_ASSERT( patch_3d.num_nodes() > 1 );
  PatchData subpatch;
  patch_3d.get_subpatch( 1, 1, subpatch, err );
  CPPUNIT_ASSERT( !err );
  
    // make sure we copied the cached values onto the subpatch
  CPPUNIT_ASSERT_EQUAL( cached->calls_3d(), 0u );
  
    // Test the values for each cached matrix on the subpatch
    // NOTE:  This test takes advantange of the fact that the
    // "handles" in the subpatch are indices into the main patch.
  
    // test each value
  for (size_t i = 0; i < subpatch.num_elements(); ++i)
  {
    MsqMeshEntity& elem = subpatch.element_by_index(i);

    get_samples( elem.get_element_type(), samples, locations );
    for (unsigned j = 0; j < locations.size(); ++j)
   {
      MsqMatrix<3,3> W;
      bool rval = cacher->get_3D_target( subpatch, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval && !err);
      
      Mesh::ElementHandle h = subpatch.get_element_handles_array()[i];
      Mesh::ElementHandle* old_h = patch_3d.get_element_handles_array();
      size_t old_idx = msq_std::find( old_h, old_h + patch_3d.num_elements(), h ) - old_h;
      CPPUNIT_ASSERT(old_idx < patch_3d.num_elements());
      MsqMatrix<3,3> M = CachedTargetCalculator::make_3d( old_idx, locations[j]);
      ASSERT_MATRICES_EQUAL( W, M, DBL_EPSILON );
    }
  }
}

void CachingTargetTest::test_2d_target_subpatch()
{
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
    // cache some values on the main patch
  request_all_targets_2d();
  
    // clear the count so we know if any additional 
    // evalutions of the base target calculator are
    // done during subpatch creation.
  cached->clear();
  
    // create a sub-patch
  CPPUNIT_ASSERT( patch_2d.num_nodes() > 1 );
  PatchData subpatch;
  patch_2d.get_subpatch( 1, 1, subpatch, err );
  CPPUNIT_ASSERT( !err );
  
    // make sure we copied the cached values onto the subpatch
  CPPUNIT_ASSERT_EQUAL( cached->calls_2d(), 0u );
  
    // Test the values for each cached matrix on the subpatch
    // NOTE:  This test takes advantange of the fact that the
    // "handles" in the subpatch are indices into the main patch.
  
    // test each value
  for (size_t i = 0; i < subpatch.num_elements(); ++i)
  {
    MsqMeshEntity& elem = subpatch.element_by_index(i);

    get_samples( elem.get_element_type(), samples, locations );
    for (unsigned j = 0; j < locations.size(); ++j)
    {
      MsqMatrix<3,2> W;
      bool rval = cacher->get_2D_target( subpatch, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval && !err);
      
      Mesh::ElementHandle h = subpatch.get_element_handles_array()[i];
      Mesh::ElementHandle* old_h = patch_2d.get_element_handles_array();
      size_t old_idx = msq_std::find( old_h, old_h + patch_2d.num_elements(), h ) - old_h;
      CPPUNIT_ASSERT(old_idx < patch_2d.num_elements());
      MsqMatrix<3,2> M = CachedTargetCalculator::make_2d( old_idx, locations[j] );
      ASSERT_MATRICES_EQUAL( W, M, DBL_EPSILON );
    }
  }
}

void CachingTargetTest::test_3d_surf_target_subpatch()
{
  MsqPrintError err(msq_stdio::cout);
  msq_std::vector<Sample> locations;
  
    // cache some values on the main patch
  request_all_targets_3d_surf();
  
    // clear the count so we know if any additional 
    // evalutions of the base target calculator are
    // done during subpatch creation.
  cached->clear();
  
    // create a sub-patch
  CPPUNIT_ASSERT( patch_2d.num_nodes() > 1 );
  PatchData subpatch;
  patch_2d.get_subpatch( 1, 1, subpatch, err );
  CPPUNIT_ASSERT( !err );
  
    // make sure we copied the cached values onto the subpatch
  CPPUNIT_ASSERT_EQUAL( cached_3d_surf->calls_2d(), 0u );
  
    // Test the values for each cached matrix on the subpatch
    // NOTE:  This test takes advantange of the fact that the
    // "handles" in the subpatch are indices into the main patch.
  
    // test each value
  for (size_t i = 0; i < subpatch.num_elements(); ++i)
  {
    MsqMeshEntity& elem = subpatch.element_by_index(i);

    get_samples( elem.get_element_type(), samples, locations );
    for (unsigned j = 0; j < locations.size(); ++j)
    {
      MsqMatrix<3,3> W;
      bool rval = cacher_3d_surf->get_3D_target( subpatch, i, samples, locations[j], W, err );
      CPPUNIT_ASSERT(rval && !err);
      
      Mesh::ElementHandle h = subpatch.get_element_handles_array()[i];
      Mesh::ElementHandle* old_h = patch_2d.get_element_handles_array();
      size_t old_idx = msq_std::find( old_h, old_h + patch_2d.num_elements(), h ) - old_h;
      CPPUNIT_ASSERT(old_idx < patch_2d.num_elements());
      MsqMatrix<3,3> M = CachedTargetCalculator::make_3d( old_idx, locations[j] );
      ASSERT_MATRICES_EQUAL( W, M, DBL_EPSILON );
    }
  }
}

void CachingTargetTest::test_cache_cleared()
{
  MsqPrintError err(msq_stdio::cout);
  
    // cache some values on the main patch
  request_all_targets_3d();
  
    // clear the count so we know if any additional 
    // evalutions of the base target calculator are
    // done.
  cached->clear();
  
    // now re-create the patch, which should result in the
    // cached data being notified that the mesh has changed
  create_twelve_hex_patch( patch_3d, err );
  CPPUNIT_ASSERT(!err);
  
    // now get cached values for each element
  unsigned count = request_all_targets_3d();
  
    // and check that they were all re-calculated
  CPPUNIT_ASSERT_EQUAL( cached->calls_3d(), count );
}
  
  
