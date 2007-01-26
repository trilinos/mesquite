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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file I_DFT_Test.cpp
 *  \brief unit tests for _IDFT quality metrics
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "I_DFT.hpp"
#include "I_DFT_WeakBarrier.hpp"
#include "I_DFT_NoBarrier.hpp"
#include "I_DFT_StrongBarrier.hpp"
#include "I_DFT_InverseMeanRatio.hpp"
#include "I_DFT_Generalized.hpp"
#include "IdealCornerTarget.hpp"
#include "UnitWeight.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "QualityMetricTester.hpp"

using namespace Mesquite;

class I_DFT_Test : public CppUnit::TestFixture
{
private:  

  CPPUNIT_TEST_SUITE(I_DFT_Test);

  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_indices);
  CPPUNIT_TEST (test_get_fixed_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_one_vertex_gradient);
  CPPUNIT_TEST (test_one_vertex_Hessian);
  
  CPPUNIT_TEST (test_ideal_element_eval_I_DFT);
  CPPUNIT_TEST (test_ideal_element_grad_I_DFT);
  CPPUNIT_TEST (test_ideal_element_hess_I_DFT);
  CPPUNIT_TEST (test_valid_hessian_I_DFT);
  CPPUNIT_TEST (test_measures_quality_I_DFT);
  CPPUNIT_TEST (test_eval_with_gradient_I_DFT);
  CPPUNIT_TEST (test_eval_with_hessian_I_DFT);
  CPPUNIT_TEST (test_gradient_reflects_quality_I_DFT);
  CPPUNIT_TEST (test_domain_deviation_I_DFT);
  //CPPUNIT_TEST (test_inverted_elements_I_DFT);
  //CPPUNIT_TEST (test_degenerate_elements_I_DFT);
  CPPUNIT_TEST (test_location_invariant_I_DFT);
  //CPPUNIT_TEST (test_scale_invariant_I_DFT);
  //CPPUNIT_TEST (test_orient_invariant_I_DFT);
  
  CPPUNIT_TEST (test_ideal_element_eval_no_barrier);
  CPPUNIT_TEST (test_ideal_element_grad_no_barrier);
  CPPUNIT_TEST (test_ideal_element_hess_no_barrier);
  CPPUNIT_TEST (test_valid_hessian_no_barrier);
  CPPUNIT_TEST (test_measures_quality_no_barrier);
  CPPUNIT_TEST (test_eval_with_gradient_no_barrier);
  CPPUNIT_TEST (test_eval_with_hessian_no_barrier);
  CPPUNIT_TEST (test_gradient_reflects_quality_no_barrier);
  //CPPUNIT_TEST (test_domain_deviation_no_barrier);
  CPPUNIT_TEST (test_inverted_elements_no_barrier);
  CPPUNIT_TEST (test_degenerate_elements_no_barrier);
  CPPUNIT_TEST (test_location_invariant_no_barrier);
  //CPPUNIT_TEST (test_scale_invariant_no_barrier);
  //CPPUNIT_TEST (test_orient_invariant_no_barrier);
  
  CPPUNIT_TEST (test_ideal_element_eval_weak_barrier);
  CPPUNIT_TEST (test_ideal_element_grad_weak_barrier);
  CPPUNIT_TEST (test_ideal_element_hess_weak_barrier);
  CPPUNIT_TEST (test_valid_hessian_weak_barrier);
  CPPUNIT_TEST (test_measures_quality_weak_barrier);
  CPPUNIT_TEST (test_eval_with_gradient_weak_barrier);
  CPPUNIT_TEST (test_eval_with_hessian_weak_barrier);
  CPPUNIT_TEST (test_gradient_reflects_quality_weak_barrier);
  CPPUNIT_TEST (test_domain_deviation_weak_barrier);
  //CPPUNIT_TEST (test_inverted_elements_weak_barrier);
  //CPPUNIT_TEST (test_degenerate_elements_weak_barrier);
  CPPUNIT_TEST (test_location_invariant_weak_barrier);
  //CPPUNIT_TEST (test_scale_invariant_weak_barrier);
  //CPPUNIT_TEST (test_orient_invariant_weak_barrier);
  
  CPPUNIT_TEST (test_ideal_element_eval_strong_barrier);
  CPPUNIT_TEST (test_ideal_element_grad_strong_barrier);
  CPPUNIT_TEST (test_ideal_element_hess_strong_barrier);
  CPPUNIT_TEST (test_valid_hessian_strong_barrier);
  CPPUNIT_TEST (test_measures_quality_strong_barrier);
  CPPUNIT_TEST (test_eval_with_gradient_strong_barrier);
  CPPUNIT_TEST (test_eval_with_hessian_strong_barrier);
  CPPUNIT_TEST (test_gradient_reflects_quality_strong_barrier);
  CPPUNIT_TEST (test_domain_deviation_strong_barrier);
  //CPPUNIT_TEST (test_inverted_elements_strong_barrier);
  //CPPUNIT_TEST (test_degenerate_elements_strong_barrier);
  CPPUNIT_TEST (test_location_invariant_strong_barrier);
  //CPPUNIT_TEST (test_scale_invariant_strong_barrier);
  //CPPUNIT_TEST (test_orient_invariant_strong_barrier);
  
  CPPUNIT_TEST (test_ideal_element_eval_imr);
  CPPUNIT_TEST (test_ideal_element_grad_imr);
  CPPUNIT_TEST (test_ideal_element_hess_imr);
  CPPUNIT_TEST (test_valid_hessian_imr);
  CPPUNIT_TEST (test_measures_quality_imr);
  CPPUNIT_TEST (test_eval_with_gradient_imr);
  CPPUNIT_TEST (test_eval_with_hessian_imr);
  CPPUNIT_TEST (test_gradient_reflects_quality_imr);
  CPPUNIT_TEST (test_domain_deviation_imr);
  CPPUNIT_TEST (test_inverted_elements_imr);
  CPPUNIT_TEST (test_degenerate_elements_imr);
  CPPUNIT_TEST (test_location_invariant_imr);
  //CPPUNIT_TEST (test_scale_invariant_imr);
  CPPUNIT_TEST (test_orient_invariant_imr);
  
  CPPUNIT_TEST_SUITE_END();
  
  IdealCornerTarget tc;
  UnitWeight wc;
  I_DFT i_dft_metric;
  I_DFT_WeakBarrier weak_barrier_metric;
  I_DFT_StrongBarrier strong_barrier_metric;
  I_DFT_NoBarrier no_barrier_metric;
  I_DFT_InverseMeanRatio imr_metric;
  QualityMetricTester tester;
  
public:

  I_DFT_Test() 
    : i_dft_metric( &tc, &wc ),
      weak_barrier_metric( &tc, &wc ),
      strong_barrier_metric( &tc, &wc ),
      no_barrier_metric( &tc, &wc ),
      imr_metric( &tc, &wc ),
      tester(QualityMetricTester::ALL_FE_EXCEPT_SEPTAHEDRON)
    { tester.ideal_pyramid_base_equals_height( true ); }

  void test_supported_types()
    { tester.test_supported_element_types( &i_dft_metric ); }

  void test_get_evaluations()
    { tester.test_get_element_evaluations( &i_dft_metric ); }
    
  void test_get_element_indices()
    { tester.test_get_element_indices( &i_dft_metric ); }
  
  void test_get_fixed_indices()
    { tester.test_get_indices_fixed( &i_dft_metric ); }
  
  void test_eval_with_indices()
    { tester.compare_eval_and_eval_with_indices( &i_dft_metric ); }

  void test_one_vertex_gradient()
    { tester.test_gradient_with_fixed_vertex( &i_dft_metric ); }

  void test_one_vertex_Hessian()
    { tester.test_hessian_with_fixed_vertex( &i_dft_metric ); }

    
  void test_ideal_element_eval_I_DFT()
  {
    tester.test_evaluate_unit_edge_element( &i_dft_metric, TRIANGLE, 0.0 );
    tester.test_evaluate_unit_edge_element( &i_dft_metric, QUADRILATERAL, 0.0 );
    tester.test_evaluate_unit_edge_element( &i_dft_metric, TETRAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &i_dft_metric, HEXAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &i_dft_metric, PRISM, 0.0 );
    tester.test_evaluate_unit_edge_element( &i_dft_metric, PYRAMID, 0.0 );
  }
    
  void test_ideal_element_eval_no_barrier()
  {
    tester.test_evaluate_unit_edge_element( &no_barrier_metric, TRIANGLE, 0.0 );
    tester.test_evaluate_unit_edge_element( &no_barrier_metric, QUADRILATERAL, 0.0 );
    tester.test_evaluate_unit_edge_element( &no_barrier_metric, TETRAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &no_barrier_metric, HEXAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &no_barrier_metric, PRISM, 0.0 );
    tester.test_evaluate_unit_edge_element( &no_barrier_metric, PYRAMID, 0.0 );
  }
    
  void test_ideal_element_eval_weak_barrier()
  {
    tester.test_evaluate_unit_edge_element( &weak_barrier_metric, TRIANGLE, 0.0 );
    tester.test_evaluate_unit_edge_element( &weak_barrier_metric, QUADRILATERAL, 0.0 );
    tester.test_evaluate_unit_edge_element( &weak_barrier_metric, TETRAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &weak_barrier_metric, HEXAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &weak_barrier_metric, PRISM, 0.0 );
    tester.test_evaluate_unit_edge_element( &weak_barrier_metric, PYRAMID, 0.0 );
  }
    
  void test_ideal_element_eval_strong_barrier()
  {
    tester.test_evaluate_unit_edge_element( &strong_barrier_metric, TRIANGLE, 0.0 );
    tester.test_evaluate_unit_edge_element( &strong_barrier_metric, QUADRILATERAL, 0.0 );
    tester.test_evaluate_unit_edge_element( &strong_barrier_metric, TETRAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &strong_barrier_metric, HEXAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &strong_barrier_metric, PRISM, 0.0 );
    tester.test_evaluate_unit_edge_element( &strong_barrier_metric, PYRAMID, 0.0 );
  }
    
  void test_ideal_element_eval_imr()
  {
    tester.test_evaluate_unit_edge_element( &imr_metric, TRIANGLE, 1.0 );
    tester.test_evaluate_unit_edge_element( &imr_metric, QUADRILATERAL, 1.0 );
    tester.test_evaluate_unit_edge_element( &imr_metric, TETRAHEDRON, 1.0 );
    tester.test_evaluate_unit_edge_element( &imr_metric, HEXAHEDRON, 1.0 );
    tester.test_evaluate_unit_edge_element( &imr_metric, PRISM, 1.0 );
    tester.test_evaluate_unit_edge_element( &imr_metric, PYRAMID, 1.0 );
  }
  
  void test_ideal_element_grad_I_DFT()
    { tester.test_ideal_element_zero_gradient( &i_dft_metric, false ); }
  void test_ideal_element_grad_no_barrier()
    { tester.test_ideal_element_zero_gradient( &no_barrier_metric, false ); }
  void test_ideal_element_grad_weak_barrier()
    { tester.test_ideal_element_zero_gradient( &weak_barrier_metric, false ); }
  void test_ideal_element_grad_strong_barrier()
    { tester.test_ideal_element_zero_gradient( &strong_barrier_metric, false ); }
  void test_ideal_element_grad_imr()
    { tester.test_ideal_element_zero_gradient( &imr_metric, false ); }
  
  void test_ideal_element_hess_I_DFT()
    { tester.test_ideal_element_positive_definite_Hessian( &i_dft_metric, false ); }
  void test_ideal_element_hess_no_barrier()
    { tester.test_ideal_element_positive_definite_Hessian( &no_barrier_metric, false ); }
  void test_ideal_element_hess_weak_barrier()
    { tester.test_ideal_element_positive_definite_Hessian( &weak_barrier_metric, false ); }
  void test_ideal_element_hess_strong_barrier()
    { tester.test_ideal_element_positive_definite_Hessian( &strong_barrier_metric, false ); }
  void test_ideal_element_hess_imr()
    { tester.test_ideal_element_positive_definite_Hessian( &imr_metric, false ); }
    
  void test_valid_hessian_I_DFT()
    { tester.test_symmetric_Hessian_diagonal_blocks( &i_dft_metric ); }
  void test_valid_hessian_no_barrier()
    { tester.test_symmetric_Hessian_diagonal_blocks( &no_barrier_metric ); }
  void test_valid_hessian_weak_barrier()
    { tester.test_symmetric_Hessian_diagonal_blocks( &weak_barrier_metric ); }
  void test_valid_hessian_strong_barrier()
    { tester.test_symmetric_Hessian_diagonal_blocks( &strong_barrier_metric ); }
  void test_valid_hessian_imr()
    { tester.test_symmetric_Hessian_diagonal_blocks( &imr_metric ); }
  
  void test_measures_quality_I_DFT()
    { tester.test_measures_quality( &i_dft_metric ); }
  void test_measures_quality_no_barrier()
    { tester.test_measures_quality( &no_barrier_metric ); }
  void test_measures_quality_weak_barrier()
    { tester.test_measures_quality( &weak_barrier_metric ); }
  void test_measures_quality_strong_barrier()
    { tester.test_measures_quality( &strong_barrier_metric ); }
  void test_measures_quality_imr()
    { tester.test_measures_quality( &imr_metric ); }
  
  void test_gradient_reflects_quality_I_DFT()
    { tester.test_gradient_reflects_quality( &i_dft_metric ); }
  void test_gradient_reflects_quality_no_barrier()
    { tester.test_gradient_reflects_quality( &no_barrier_metric ); }
  void test_gradient_reflects_quality_weak_barrier()
    { tester.test_gradient_reflects_quality( &weak_barrier_metric ); }
  void test_gradient_reflects_quality_strong_barrier()
    { tester.test_gradient_reflects_quality( &strong_barrier_metric ); }
  void test_gradient_reflects_quality_imr()
    { tester.test_gradient_reflects_quality( &imr_metric ); }
  
  void test_domain_deviation_I_DFT()
  {
    tester.test_domain_deviation_quality( &i_dft_metric );
    tester.test_domain_deviation_gradient( &i_dft_metric );
  }
  void test_domain_deviation_no_barrier()
  {
    tester.test_domain_deviation_quality( &no_barrier_metric );
    tester.test_domain_deviation_gradient( &no_barrier_metric );
  }
  void test_domain_deviation_weak_barrier()
  {
    tester.test_domain_deviation_quality( &weak_barrier_metric );
    tester.test_domain_deviation_gradient( &weak_barrier_metric );
  }
  void test_domain_deviation_strong_barrier()
  {
    tester.test_domain_deviation_quality( &strong_barrier_metric );
    tester.test_domain_deviation_gradient( &strong_barrier_metric );
  }
  void test_domain_deviation_imr()
  {
    tester.test_domain_deviation_quality( &imr_metric );
    tester.test_domain_deviation_gradient( &imr_metric );
  }
  
  void test_inverted_elements_I_DFT()
    { tester.test_evaluate_inverted_element( &i_dft_metric, false ); }
  void test_inverted_elements_no_barrier()
    { tester.test_evaluate_inverted_element( &no_barrier_metric, true ); }
  void test_inverted_elements_weak_barrier()
    { tester.test_evaluate_inverted_element( &weak_barrier_metric, false ); }
  void test_inverted_elements_strong_barrier()
    { tester.test_evaluate_inverted_element( &strong_barrier_metric, false ); }
  void test_inverted_elements_imr()
    { tester.test_evaluate_inverted_element( &imr_metric, false ); }
    
  void test_degenerate_elements_I_DFT()
    { tester.test_evaluate_degenerate_element( &i_dft_metric, false ); }
  void test_degenerate_elements_no_barrier()
    { tester.test_evaluate_degenerate_element( &no_barrier_metric, true ); }
  void test_degenerate_elements_weak_barrier()
    { tester.test_evaluate_degenerate_element( &weak_barrier_metric, false ); }
  void test_degenerate_elements_strong_barrier()
    { tester.test_evaluate_degenerate_element( &strong_barrier_metric, false ); }
  void test_degenerate_elements_imr()
    { tester.test_evaluate_degenerate_element( &imr_metric, false ); }
    
  void test_eval_with_gradient_I_DFT()
  {
    tester.compare_eval_with_indices_and_eval_with_gradient( &i_dft_metric );
    tester.compare_analytical_and_numerical_gradients( &i_dft_metric );
  }
  void test_eval_with_gradient_no_barrier()
  {
    tester.compare_eval_with_indices_and_eval_with_gradient( &no_barrier_metric );
    tester.compare_analytical_and_numerical_gradients( &no_barrier_metric );
  }
  void test_eval_with_gradient_weak_barrier()
  {
    tester.compare_eval_with_indices_and_eval_with_gradient( &weak_barrier_metric );
    tester.compare_analytical_and_numerical_gradients( &weak_barrier_metric );
  }
  void test_eval_with_gradient_strong_barrier()
  {
    tester.compare_eval_with_indices_and_eval_with_gradient( &strong_barrier_metric );
    tester.compare_analytical_and_numerical_gradients( &strong_barrier_metric );
  }
  void test_eval_with_gradient_imr()
  {
    tester.compare_eval_with_indices_and_eval_with_gradient( &imr_metric );
    tester.compare_analytical_and_numerical_gradients( &imr_metric );
  }
  
  void test_eval_with_hessian_I_DFT()
  {
    tester.compare_eval_with_indices_and_eval_with_hessian( &i_dft_metric );
    tester.compare_eval_with_grad_and_eval_with_hessian( &i_dft_metric );
    tester.compare_analytical_and_numerical_hessians( &i_dft_metric );
  }
  void test_eval_with_hessian_no_barrier()
  {
    tester.compare_eval_with_indices_and_eval_with_hessian( &no_barrier_metric );
    tester.compare_eval_with_grad_and_eval_with_hessian( &no_barrier_metric );
    tester.compare_analytical_and_numerical_hessians( &no_barrier_metric );
  }
  void test_eval_with_hessian_weak_barrier()
  {
    tester.compare_eval_with_indices_and_eval_with_hessian( &weak_barrier_metric );
    tester.compare_eval_with_grad_and_eval_with_hessian( &weak_barrier_metric );
    tester.compare_analytical_and_numerical_hessians( &weak_barrier_metric );
  }
  void test_eval_with_hessian_strong_barrier()
  {
    tester.compare_eval_with_indices_and_eval_with_hessian( &strong_barrier_metric );
    tester.compare_eval_with_grad_and_eval_with_hessian( &strong_barrier_metric );
    tester.compare_analytical_and_numerical_hessians( &strong_barrier_metric );
  }
  void test_eval_with_hessian_imr()
  {
    tester.compare_eval_with_indices_and_eval_with_hessian( &imr_metric );
    tester.compare_eval_with_grad_and_eval_with_hessian( &imr_metric );
    tester.compare_analytical_and_numerical_hessians( &imr_metric );
  }
  
  void test_location_invariant_I_DFT()
  {
    tester.test_location_invariant( &i_dft_metric );
    tester.test_grad_location_invariant( &i_dft_metric );
    tester.test_hessian_location_invariant( &i_dft_metric );
  }
  void test_location_invariant_no_barrier()
  {
    tester.test_location_invariant( &no_barrier_metric );
    tester.test_grad_location_invariant( &no_barrier_metric );
    tester.test_hessian_location_invariant( &no_barrier_metric );
  }
  void test_location_invariant_weak_barrier()
  {
    tester.test_location_invariant( &weak_barrier_metric );
    tester.test_grad_location_invariant( &weak_barrier_metric );
    tester.test_hessian_location_invariant( &weak_barrier_metric );
  }
  void test_location_invariant_strong_barrier()
  {
    tester.test_location_invariant( &strong_barrier_metric );
    tester.test_grad_location_invariant( &strong_barrier_metric );
    tester.test_hessian_location_invariant( &strong_barrier_metric );
  }
  void test_location_invariant_imr()
  {
    tester.test_location_invariant( &imr_metric );
    tester.test_grad_location_invariant( &imr_metric );
    tester.test_hessian_location_invariant( &imr_metric );
  }
  
  void test_scale_invariant_I_DFT()
  {
    tester.test_scale_invariant( &i_dft_metric );
  }
  void test_scale_invariant_no_barrier()
  {
    tester.test_scale_invariant( &no_barrier_metric );
  }
  void test_scale_invariant_weak_barrier()
  {
    tester.test_scale_invariant( &weak_barrier_metric );
  }
  void test_scale_invariant_strong_barrier()
  {
    tester.test_scale_invariant( &strong_barrier_metric );
  }
  void test_scale_invariant_imr()
  {
    tester.test_scale_invariant( &imr_metric );
  }
  
  void test_orient_invariant_I_DFT()
  {
    tester.test_orient_invariant( &i_dft_metric );
    tester.test_grad_orient_invariant( &i_dft_metric );
  }
  void test_orient_invariant_no_barrier()
  {
    tester.test_orient_invariant( &no_barrier_metric );
    tester.test_grad_orient_invariant( &no_barrier_metric );
  }
  void test_orient_invariant_weak_barrier()
  {
    tester.test_orient_invariant( &weak_barrier_metric );
    tester.test_grad_orient_invariant( &weak_barrier_metric );
  }
  void test_orient_invariant_strong_barrier()
  {
    tester.test_orient_invariant( &strong_barrier_metric );
    tester.test_grad_orient_invariant( &strong_barrier_metric );
  }
  void test_orient_invariant_imr()
  {
    tester.test_orient_invariant( &imr_metric );
    tester.test_grad_orient_invariant( &imr_metric );
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(I_DFT_Test, "I_DFT_Test");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(I_DFT_Test, "Unit");
