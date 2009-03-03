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


/** \file ElemSampleQMTest.cpp
 *  \brief Unit tests for ElementSampleQM class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ElemSampleQM.hpp"
#include <cppunit/extensions/HelperMacros.h>

using namespace Mesquite;

class ElemSampleQMTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(ElemSampleQMTest);
  CPPUNIT_TEST (test_handle_from_sample);
  CPPUNIT_TEST (test_handle_from_side);
  CPPUNIT_TEST (test_sample_from_side);
  CPPUNIT_TEST_SUITE_END();

public:
  void test_handle_from_sample();
  void test_handle_from_side();
  void test_sample_from_side();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ElemSampleQMTest, "ElemSampleQMTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ElemSampleQMTest, "Unit");

void test_handle( size_t sample, size_t elem )
{
  size_t handle = ElemSampleQM::handle( sample, elem );
  CPPUNIT_ASSERT_EQUAL( sample, (size_t)ElemSampleQM::sample( handle ) );
  CPPUNIT_ASSERT_EQUAL( elem,   (size_t)ElemSampleQM::elem( handle ) );
}

void ElemSampleQMTest::test_handle_from_sample( )
{
  size_t min_sample_no = 0;
  size_t max_sample_no = ElemSampleQM::MAX_SAMPLES_PER_ELEM - 1;
  size_t min_elem_no = 0;
  size_t max_elem_no = ElemSampleQM::MAX_ELEM_PER_PATCH - 1;
  
  test_handle( min_sample_no, min_elem_no );
  test_handle( min_sample_no, max_elem_no );
  test_handle( max_sample_no, min_elem_no );
  test_handle( max_sample_no, max_elem_no );
}

void test_handle( unsigned side_dim, unsigned side_num, size_t elem )
{
  size_t handle = ElemSampleQM::handle( side_dim, side_num, elem );
  CPPUNIT_ASSERT_EQUAL( side_dim, ElemSampleQM::side_dim_from_handle( handle ) );
  CPPUNIT_ASSERT_EQUAL( side_num, ElemSampleQM::side_num_from_handle( handle ) );
  CPPUNIT_ASSERT_EQUAL( elem, (size_t)ElemSampleQM::elem( handle ) );
}

void ElemSampleQMTest::test_handle_from_side()
{
  size_t min_elem_no = 0;
  size_t max_elem_no = ElemSampleQM::MAX_ELEM_PER_PATCH - 1;
  unsigned min_side_dim = 0;
  unsigned max_side_dim = 3;
  unsigned min_side_no = 0;
  unsigned max_side_no = (1<<ElemSampleQM::SIDE_NUMBER_BITS) - 1;
  
  test_handle( min_side_dim, min_side_no, min_elem_no ); 
  test_handle( min_side_dim, min_side_no, max_elem_no ); 
  test_handle( min_side_dim, max_side_no, min_elem_no ); 
  test_handle( min_side_dim, max_side_no, max_elem_no ); 
  test_handle( max_side_dim, min_side_no, min_elem_no ); 
  test_handle( max_side_dim, min_side_no, max_elem_no ); 
  test_handle( max_side_dim, max_side_no, min_elem_no ); 
  test_handle( max_side_dim, max_side_no, max_elem_no ); 
}

void test_sample( unsigned side_dim, unsigned side_num )
{
  unsigned sample = ElemSampleQM::sample( side_dim, side_num );
  CPPUNIT_ASSERT_EQUAL( side_dim, ElemSampleQM::side_dim_from_sample( sample ) );
  CPPUNIT_ASSERT_EQUAL( side_num, ElemSampleQM::side_num_from_sample( sample ) );
}


void ElemSampleQMTest::test_sample_from_side()
{
  unsigned min_side_dim = 0;
  unsigned max_side_dim = 3;
  unsigned min_side_no = 0;
  unsigned max_side_no = (1<<ElemSampleQM::SIDE_NUMBER_BITS) - 1;
  
  test_sample( min_side_dim, min_side_no ); 
  test_sample( min_side_dim, max_side_no ); 
  test_sample( max_side_dim, min_side_no ); 
  test_sample( max_side_dim, max_side_no ); 
}
