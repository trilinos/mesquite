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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Nov-02 at 18:05:56
//  LAST-MOD: 23-Jul-03 at 17:40:28 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file ObjectiveFunctionTest.cpp

Unit testing of various functions in the ObjectiveFunction class. 
*/
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "ObjectiveFunction.hpp"
#include "LPtoPTemplate.hpp"
#include "MaxTemplate.hpp"
#include "LInfTemplate.hpp"
#include "CompositeOFAdd.hpp"
#include "CompositeOFMultiply.hpp"
#include "CompositeOFScalarMultiply.hpp"
#include "CompositeOFScalarAdd.hpp"
#include "GeneralizedConditionNumberQualityMetric.hpp"
#include "MeanRatioQualityMetric.hpp"
#include "InverseMeanRatioQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "MsqHessian.hpp"

#include "PatchDataInstances.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "MsqFreeVertexIndexIterator.hpp"
#include <list>
#include <iterator>

using namespace Mesquite;
using std::cout;
using std::endl;
using std::cerr;

class ObjectiveFunctionTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(ObjectiveFunctionTest);
  CPPUNIT_TEST (test_get_quality_metric_list);
  CPPUNIT_TEST (test_max_templates);
  CPPUNIT_TEST (test_compute_gradient_3D_LPtoPTemplate_L1_hex);
  CPPUNIT_TEST (test_compute_gradient_3D_LPtoPTemplate_L1_hex_inverse);
  CPPUNIT_TEST (test_compute_gradient_3D_LPtoPTemplate_L2_hex);
  CPPUNIT_TEST (test_compute_gradient_3D_LPtoPTemplate_L2_hex_scaled);
  CPPUNIT_TEST (test_compute_ana_hessian_tet);
  CPPUNIT_TEST (test_compute_ana_hessian_tet_scaled);
  CPPUNIT_TEST (test_compute_gradient3D_composite);
  CPPUNIT_TEST (test_OFval_from_evaluate_and_gradient_LPtoP);
  CPPUNIT_TEST (test_grad_from_gradient_and_hessian_LPtoP);
  CPPUNIT_TEST (test_grad_from_gradient_and_hessian_LPtoP_negate);
  CPPUNIT_TEST (test_LPtoP_negate_flag);

  CPPUNIT_TEST_SUITE_END();
   
private:
   
  PatchData m4Quads;
  PatchData m6Quads;
  PatchData m12Hex;
  PatchData triPatch;
  PatchData tetPatch;

public:
  void setUp()
  {
    MsqError err;
     
    /* our 2D set up: 4 quads, center vertex outcentered by (0,-0.5)
       7____6____5
       |    |    |
       | 2  |  3 |
       8-_  |  _-4       vertex 1 is at (0,0)
       |  -_0_-  |       vertex 5 is at (2,2)
       | 0  |  1 |
       1----2----3
    */
    create_four_quads_patch(m4Quads, err); MSQ_CHKERR(err);

    /*! \fn create_six_quads_patch(PatchData &four_quads, MsqError &err)
      our 2D set up: 6 quads, 1 center vertex outcentered by (0,-0.5), the other centered
      7____6____5___11
      |    |    |    |
      | 2  |  3 | 5  |
      8-_  |  _-4---10       vertex 1 is at (0,0)
      |  -_0_-  |    |       vertex 11 is at (3,2)
      | 0  |  1 | 4  |
      1----2----3----9
    */
    create_six_quads_patch_with_domain(m6Quads, err); MSQ_CHKERR(err);

    /*! \fn create_twelve_hex_patch(PatchData &pd, MsqError &err)
      3D set up: 12 quads, one center vertex outcentered by (0,-0.5),
      the other centered. Vertex 1 is at (0,0,-1). Vertex 35 is at (3,2,1).
     
      7____6____5___11     19___18____17__23     31___30___29___35
      |    |    |    |      |    |    |    |      |    |    |    |
      | 2  |  3 | 5  |      |    |    |    |      | 8  |  9 | 11 |
      8----0----4---10     20-_  |  _16---22     32---24---28---34       
      |    |    |    |      |  -12_-  |    |      |    |    |    |       
      | 0  |  1 | 4  |      |    |    |    |      | 6  |  7 | 10 |
      1----2----3----9     13---14---15---21     25---26---27---33
    */
    create_twelve_hex_patch(m12Hex, err); MSQ_CHKERR(err);

   /*! \fn create_two_tri_patch(PatchData &one_tri_patch, MsqError &err)
            2
           / \      creates a Patch containing two ideal triangles
          / 0 \
         0-----1
          \ 1 /
           \ /
            3
   */
    create_qm_two_tri_patch_with_domain(triPatch,err);MSQ_CHKERR(err);
    
    create_qm_two_tet_patch(tetPatch,err);MSQ_CHKERR(err);
    
  }

  void tearDown()
  {
    destroy_patch_with_domain(triPatch);
    destroy_patch_with_domain(m6Quads);
  }
  
public:
  ObjectiveFunctionTest()
  {}
  
  void test_get_quality_metric_list()
  {
    MsqError err;
      
    // instantiates a couple of QualityMetrics
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    ShapeQualityMetric* condition_nb = new GeneralizedConditionNumberQualityMetric;

    // and creates a composite objective function.
    LPtoPTemplate* LP2_mean_ratio = new LPtoPTemplate(mean_ratio, 2, err); MSQ_CHKERR(err);
    LPtoPTemplate* LP2_condition_nb = new LPtoPTemplate(condition_nb, 2, err); MSQ_CHKERR(err);
    CompositeOFScalarMultiply* LP2_condition_nb_x3 = new CompositeOFScalarMultiply(3,LP2_condition_nb);   
    CompositeOFAdd comp_OF(LP2_mean_ratio, LP2_condition_nb_x3);

    // test the (simple) get_quality_metric function, on OF with 1 QM. 
    QualityMetric* QM_;
    QM_ = LP2_mean_ratio->get_quality_metric();
    CPPUNIT_ASSERT(QM_==mean_ratio);
    // test the (simple) get_quality_metric function, on OF with 2 QMs. 
    QM_ = comp_OF.get_quality_metric();
    CPPUNIT_ASSERT(QM_==NULL);

    // test the get_quality_metric_list function, on OF with 1 QM. 
    std::list<QualityMetric*> QM__;
    QM__ = LP2_mean_ratio->get_quality_metric_list();
    CPPUNIT_ASSERT( QM__.size()==1 );
    CPPUNIT_ASSERT( *(QM__.begin())==mean_ratio );
    // test the get_quality_metric_list function, on OF with 2 QMs. 
    QM__.clear();
    QM__ = comp_OF.get_quality_metric_list();
    CPPUNIT_ASSERT( QM__.size()==2 );
    std::list<QualityMetric*>::const_iterator QM2 = QM__.begin();
    std::list<QualityMetric*>::const_iterator QM1 = QM2++;
    CPPUNIT_ASSERT( *QM1==mean_ratio   || *QM2==mean_ratio &&
                    *QM1==condition_nb || *QM2==condition_nb  );

    delete mean_ratio;
    delete condition_nb;
    delete LP2_mean_ratio;
    delete LP2_condition_nb;
    delete LP2_condition_nb_x3;
  }
    //This function takes a patch data object and compares
    //the max template value and the l_inf value on that patch
    //for objective functions using a mean ratio metric and
    //an edge length metric.  These templates return the same
    //value assuming the metrics are non-negative.
  void compare_l_inf_and_max_templates(PatchData &pd)
     {
       MsqError err;
         // creates a mean ratio quality metric ...
       bool return_bool;
       double max_val=0.0;
       double l_inf_val = 1.0;
       ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
         //creates an edge length metric
       SmoothnessQualityMetric* smooth = new EdgeLengthQualityMetric;
       
         // ... and builds an objective function with it
       LInfTemplate l_inf_mean(mean_ratio);
       MaxTemplate max_mean(mean_ratio);
       LInfTemplate l_inf_smooth(smooth);
       MaxTemplate max_smooth(smooth);
         //check an element based metric
       return_bool=l_inf_mean.evaluate(pd, l_inf_val, err);
       CPPUNIT_ASSERT(return_bool);
       return_bool=max_mean.evaluate(pd, max_val, err);
       CPPUNIT_ASSERT(return_bool);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(l_inf_val, max_val, 1e-12);
         //check on a vertex based metric
       return_bool=l_inf_smooth.evaluate(pd, l_inf_val, err);
       CPPUNIT_ASSERT(return_bool);
       return_bool=max_smooth.evaluate(pd, max_val, err);
       CPPUNIT_ASSERT(return_bool);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(l_inf_val, max_val, 1e-12);
         //std::cout<<"\nL_INF_VAL = "<<l_inf_val<<" MAX_VAL = "<<max_val;
       
       delete mean_ratio;
       delete smooth;
     }
  void test_max_templates()
     {
       /*
           PatchData m6Quads;
           PatchData m12Hex;
           PatchData triPatch;
           PatchData tetPatch;
       */
         //test on tris
       compare_l_inf_and_max_templates(triPatch);
         //test on quads
       compare_l_inf_and_max_templates(m6Quads);
         //test on tets
       compare_l_inf_and_max_templates(tetPatch);
         //test on hexes
       compare_l_inf_and_max_templates(m12Hex);
     } 
  
  void compare_numerical_analytical_gradient(ObjectiveFunction *obj,
                                             PatchData &pd)
     {
       MsqError err;
       bool return_bool;
       double OF_val1, OF_val2;
       MsqFreeVertexIndexIterator free_ind(&pd, err);
       Vector3D* grad_num = new Vector3D[pd.num_vertices()];
       Vector3D* grad_ana = new Vector3D[pd.num_vertices()];     
       
       obj->set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
       return_bool=obj->compute_gradient(pd, grad_num, OF_val1, err);
//       std::cout<<"\nNumerical's value = "<<OF_val1;
       CPPUNIT_ASSERT(return_bool==true);
       
       int grad_pos=0;

//         free_ind.reset();
//         std::cout << "NUMERICAL GRADIENT\n";
//         for (int i=0; i<2; ++i){
//          free_ind.next();
//          grad_pos=free_ind.value();
//          for (int j=0; j<3; ++j){
//            std::cout << grad_num[grad_pos][j] << std::endl;
//          }
//        }
       
       obj->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
       return_bool=obj->compute_gradient(pd, grad_ana, OF_val2, err);
//       std::cout<<"\nAnalytical's value = "<<OF_val2;
       CPPUNIT_ASSERT(return_bool==true);
       
//        std::cout << "ANALYTICAL GRADIENT\n";
//        free_ind.reset();
//        for (int i=0; i<2; ++i){
//          free_ind.next();
//          grad_pos=free_ind.value();
//          for (int j=0; j<3; ++j){
//            std::cout << grad_ana[grad_pos][j] << std::endl;
//          }
//        }

       CPPUNIT_ASSERT_DOUBLES_EQUAL(OF_val1, OF_val2, 1e-12);
       int j;
       free_ind.reset();
       while(free_ind.next()){
         grad_pos=free_ind.value();
         for (j=0; j<3; ++j){
           CPPUNIT_ASSERT_DOUBLES_EQUAL(grad_num[grad_pos][j],
                                        grad_ana[grad_pos][j], 0.01);
         }
       } 
       delete[] grad_num;
       delete[] grad_ana;
     } 

  void test_compute_gradient_3D_LPtoPTemplate_L1_hex()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP1(mean_ratio, 1, err);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
      //    mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    compare_numerical_analytical_gradient(&LP1, m12Hex);
    delete mean_ratio;
  }
  void test_compute_gradient_3D_LPtoPTemplate_L1_hex_inverse()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* i_mean_ratio = new InverseMeanRatioQualityMetric;
    i_mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP1(i_mean_ratio, 1, err);
    i_mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
      //    mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    compare_numerical_analytical_gradient(&LP1, m12Hex);
    delete i_mean_ratio;
  }
  void test_compute_gradient_3D_LPtoPTemplate_L2_hex()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP2(mean_ratio, 2, err);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
      //    mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    compare_numerical_analytical_gradient(&LP2, m12Hex);
    delete mean_ratio;
  }
  
  void test_compute_gradient_3D_LPtoPTemplate_L2_hex_scaled()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP2(mean_ratio, 2, err);
    LP2.set_dividing_by_n(true);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
      //    mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    compare_numerical_analytical_gradient(&LP2, m12Hex);
    delete mean_ratio;
  }
  void test_compute_gradient3D_composite()
     {
       MsqError err;
         // creates a mean ratio quality metric ...
       ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
       mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
       
         // ... and builds an objective function with it
       LPtoPTemplate LP3(mean_ratio, 3, err);
       LPtoPTemplate LP2(mean_ratio,2,err);
         //build four composite objective functions
       CompositeOFScalarAdd csa_of(2,&LP2);
       CompositeOFScalarMultiply csm_of(20,&LP2);
       CompositeOFAdd ca_of(&LP3,&LP2);
       CompositeOFMultiply cm_of(&LP3,&LP2);
       
         //test scalar add
       compare_numerical_analytical_gradient(&csa_of, m12Hex);
         //test scalar multiply
       compare_numerical_analytical_gradient(&csm_of, m12Hex);
         //test add
       compare_numerical_analytical_gradient(&ca_of, m12Hex);
         //test multiply
       compare_numerical_analytical_gradient(&cm_of, m12Hex);

       delete mean_ratio;
     }
  

  void test_compute_ana_hessian_tet()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::SUM, err);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP2(mean_ratio, 2, err);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    mean_ratio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
    
    MsqHessian H;
    Vector3D* g = new Vector3D[tetPatch.num_vertices()];
    double dummy;
    H.initialize(tetPatch, err); MSQ_CHKERR(err);
    LP2.compute_hessian(tetPatch, H, g, dummy, err); MSQ_CHKERR(err);

    Matrix3D mat00(" 2.44444  0.2566   0.181444 "
		   " 0.2566   2.14815  0.104757 "
		   " 0.181444 0.104757 2.07407 ");


    Matrix3D mat13(" 5.47514 3.16659    9.83479 "
		   " -1.11704 -5.29718 -3.67406 "
		   " 10.3635 -13.5358  -15.5638 ");

    Matrix3D* mat; 

    mat = H.get_block(0,0);
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
	CPPUNIT_ASSERT_DOUBLES_EQUAL((*mat)[i][j], mat00[i][j], 1e-4);
    
    mat = H.get_block(1,3);
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
	CPPUNIT_ASSERT_DOUBLES_EQUAL((*mat)[i][j], mat13[i][j], 1e-4);
    
//    cout << H <<endl;

    delete mean_ratio;
    delete [] g;
  }
  
  void test_compute_ana_hessian_tet_scaled()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::SUM, err);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP2(mean_ratio, 2, err);
    LP2.set_dividing_by_n(true);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    mean_ratio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
    
    MsqHessian H;
    Vector3D* g = new Vector3D[tetPatch.num_vertices()];
    double dummy;
    H.initialize(tetPatch, err); MSQ_CHKERR(err);
    LP2.compute_hessian(tetPatch, H, g, dummy, err); MSQ_CHKERR(err);

    Matrix3D mat00(" 2.44444  0.2566   0.181444 "
		   " 0.2566   2.14815  0.104757 "
		   " 0.181444 0.104757 2.07407 ");

    mat00*=.5;
    
    Matrix3D mat13(" 5.47514 3.16659    9.83479 "
		   " -1.11704 -5.29718 -3.67406 "
		   " 10.3635 -13.5358  -15.5638 ");

    mat13*=.5;
    
    Matrix3D* mat; 

    mat = H.get_block(0,0);
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL((*mat)[i][j], mat00[i][j], 1e-4);
    
    mat = H.get_block(1,3);
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL((*mat)[i][j], mat13[i][j], 1e-4);
    
//    cout << H <<endl;
    delete [] g;
    delete mean_ratio;
  }
  //! Tests that the Objective function value returned from evaluate() is the
  //! same as the one returned from the gradient.
  void test_OFval_from_evaluate_and_gradient(ObjectiveFunction* OF,
                                           PatchData &pd)
  {
    MsqError err;
    bool OF_bool;
    Vector3D* grad = new Vector3D[pd.num_vertices()];

    double OF_val1;
    OF_bool = OF->evaluate(pd, OF_val1, err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(OF_bool);

    double OF_val2;
    OF_bool = OF->compute_gradient(pd, grad, OF_val2, err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(OF_bool);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(OF_val1, OF_val2, 1e-12);

    delete[] grad;
  }

  //! Calls test_OFval_from_evaluate_and_gradient() for \f$ \ell_4^4 \f$
  void test_OFval_from_evaluate_and_gradient_LPtoP()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP4(mean_ratio, 4, err);
    test_OFval_from_evaluate_and_gradient(&LP4, m12Hex);

    // ... and builds an objective function with it
    LPtoPTemplate LP1(mean_ratio, 1, err);
    test_OFval_from_evaluate_and_gradient(&LP1, m12Hex);

    SmoothnessQualityMetric* edge = new EdgeLengthQualityMetric;
    LPtoPTemplate LP5(edge, 1, err);

    test_OFval_from_evaluate_and_gradient(&LP5, m12Hex);

      //now test with scaling by N
    LP4.set_dividing_by_n(true);
    LP5.set_dividing_by_n(true);
    test_OFval_from_evaluate_and_gradient(&LP4, m12Hex);
    test_OFval_from_evaluate_and_gradient(&LP5, m12Hex);
    
    delete edge;
    delete mean_ratio;
    
  }
  
  //! Tests that the gradient value returned from compute_gradient() is the
  //! same as the one returned from the compute_hessian.
  void test_grad_from_gradient_and_hessian(ObjectiveFunction* OF,
                                           PatchData &pd)
  {
    MsqError err;
    bool OF_bool;
    double OF_val1;
    double OF_val2;
    Vector3D* grad1 = new Vector3D[pd.num_vertices()];
    Vector3D* grad2 = new Vector3D[pd.num_vertices()];
    MsqHessian hessian;
    hessian.initialize(pd, err); MSQ_CHKERR(err);

    OF_bool = OF->compute_gradient(pd, grad1, OF_val1, err);
    MSQ_CHKERR(err);
    CPPUNIT_ASSERT(OF_bool);

//     cout << "compute_gradient(pd, grad1 ...) " << endl;
//     for (int i=0; i<pd.num_vertices(); ++i)
//       cout << grad1[i];
      
    OF_bool = OF->compute_hessian(pd, hessian, grad2, OF_val2, err);
    MSQ_CHKERR(err);
    CPPUNIT_ASSERT(OF_bool);

//     cout << "\ncompute_hessian(pd, grad2 ...) " << endl;
//     for (int i=0; i<pd.num_vertices(); ++i)
//       cout << grad2[i];
    
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(OF_val1, OF_val2, 1e-8);
    for (size_t i=0; i<pd.num_vertices(); ++i){
      for (int j=0; j<3; ++j){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grad1[i][j], grad2[i][j], 5e-5);
      }
    }
    delete[] grad1;
    delete[] grad2;
  }

  //! Calls test_grad_from_gradient_and_hessian() for \f$ \ell_4^4 \f$
  void test_grad_from_gradient_and_hessian_LPtoP()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP4(mean_ratio, 4, err);
    test_grad_from_gradient_and_hessian(&LP4, m12Hex);

    // ... and builds an objective function with it
    LPtoPTemplate LP1(mean_ratio, 1, err);
    test_grad_from_gradient_and_hessian(&LP1, m12Hex);

      //test scaled versions
    LP4.set_dividing_by_n(true);
    LP1.set_dividing_by_n(true);
    test_grad_from_gradient_and_hessian(&LP4, m12Hex);
    test_grad_from_gradient_and_hessian(&LP1, m12Hex);
    delete mean_ratio;
  }
    //! Calls test_grad_from_gradient_and_hessian() for \f$ \ell_4^4 \f$
  void test_grad_from_gradient_and_hessian_LPtoP_negate()
  {
    MsqError err;
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* i_mean_ratio = new InverseMeanRatioQualityMetric;
    i_mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
    i_mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP4(i_mean_ratio, 4, err);
    test_grad_from_gradient_and_hessian(&LP4, m12Hex);

    // ... and builds an objective function with it
    LPtoPTemplate LP1(i_mean_ratio, 1, err);
    test_grad_from_gradient_and_hessian(&LP1, m12Hex);

      //test scaled versions
    LP4.set_dividing_by_n(true);
    LP1.set_dividing_by_n(true);
      //test_grad_from_gradient_and_hessian(&LP4, m12Hex);
      //test_grad_from_gradient_and_hessian(&LP1, m12Hex);
    delete i_mean_ratio;
  }
// ----------------------------------------------------------- 
// numerical objective function hessian does not work for now. 
// It will only be used for test purposes anyway. 
// ----------------------------------------------------------- 
  
//   void test_compute_hessian(PatchData &pd)
//   {
//     MsqError err;

//     MsqHessian OF_hessian_num;
//     MsqHessian OF_hessian_ana;
//     OF_hessian_num.initialize(pd, err); MSQ_CHKERR(err);
//     OF_hessian_ana.initialize(pd, err); MSQ_CHKERR(err);
    
//     // creates a mean ratio quality metric ...
//     ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
// //    mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
//     mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
//     mean_ratio->set_averaging_method(QualityMetric::SUM, err); MSQ_CHKERR(err);
// //    mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
//       //mean_ratio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);

//     // Creates an L1 objective function.
//     LPtoPTemplate L_1(mean_ratio, 1, err); MSQ_CHKERR(err);

//     // Compute numerical hessian.
//     L_1.set_gradient_type(ObjectiveFunction::NUMERICAL_GRADIENT);
//       //L_1.set_hessian_type(ObjectiveFunction::NUMERICAL_HESSIAN);
//     L_1.compute_hessian(pd, OF_hessian_num, err); MSQ_CHKERR(err);

//     cout << "Numerical OF Hessian:\n";
//     cout << OF_hessian_num << "\n\n\n";

//     // Compute analytical hessian
//     L_1.set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
//       //L_1.set_hessian_type(ObjectiveFunction::ANALYTICAL_HESSIAN);
//     L_1.compute_hessian(pd, OF_hessian_ana, err); MSQ_CHKERR(err);

//     cout << "Analytical OF Hessian:\n";
//     cout << OF_hessian_ana << endl;

//     // test
//     Matrix3D* block_num;
//     Matrix3D* block_ana;
//     CPPUNIT_ASSERT(OF_hessian_num.size() == OF_hessian_ana.size());
//     for (size_t m=0; m<OF_hessian_ana.size(); ++m) {
//       for (size_t n=m; n<OF_hessian_ana.size(); ++n) {
//         block_num = OF_hessian_num.get_block(m,n);
//         block_ana = OF_hessian_ana.get_block(m,n);
//         for (int i=0; i<3; ++i)
//           for (int j=0; j<3; ++j)
//             CPPUNIT_ASSERT_DOUBLES_EQUAL( (*block_num)[i][j], (*block_ana)[i][j], 0.001);
//       }
//     }

//     delete mean_ratio;
//   }

  
//   void test_compute_hessian_tri_patch() {
//     test_compute_hessian(triPatch);
//   }

//   void test_compute_hessian_tet_patch() {
//     test_compute_hessian(tetPatch);
//   }

    //This function tests to make sure that LPtoP handles the negate
    // flag correctly for the analytical Hessian, gradient, and evaluate.
    // It does this by creating two objective functions that are
    // identical except that the metric in one has negate flag of -1
    // and the metric in the other has a negate flag of 1.  Thus,
    // the Hessians, gradients, and function values should be the
    // the same except for a negative sign.
  void test_LPtoP_negate_flag()
  {
    size_t i;
    int j;
    bool valid;
    MsqError err;
    MsqHessian Hpos;
    Vector3D* gpos = new Vector3D[tetPatch.num_vertices()];
    double fpos;
    MsqHessian Hneg;
    Vector3D* gneg = new Vector3D[tetPatch.num_vertices()];
    double fneg;
    Hpos.initialize(tetPatch, err); MSQ_CHKERR(err);
    Hneg.initialize(tetPatch, err); MSQ_CHKERR(err);
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::LINEAR, err);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    ShapeQualityMetric* mean_ratio_neg = new MeanRatioQualityMetric;
    mean_ratio_neg->set_averaging_method(QualityMetric::LINEAR, err);
    mean_ratio_neg->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    mean_ratio_neg->set_negate_flag(-1);
    // ... and builds an objective function with it
      //FIRST TEST L4 **********************************************
    LPtoPTemplate LP4(mean_ratio, 4, err);
    LPtoPTemplate LP4_neg(mean_ratio_neg, 4, err);
      //test evaluate
    valid=LP4.evaluate(tetPatch, fpos, err);MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);
    valid=LP4_neg.evaluate(tetPatch, fneg, err);MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);
      //std::cout<<"\nFrom eval Orig fpos = "<<fpos<<" Mod fneg = "<<fneg;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(fpos,-fneg,1.e-12);

    valid=LP4.compute_gradient(tetPatch, gpos, fpos, err);
    MSQ_CHKERR(err); CPPUNIT_ASSERT(valid);
    valid=LP4_neg.compute_gradient(tetPatch, gneg, fneg, err);
    MSQ_CHKERR(err);CPPUNIT_ASSERT(valid);
      //std::cout<<"\nOrig fpos = "<<fpos<<" Mod fneg = "<<fneg;
      //test function value
    CPPUNIT_ASSERT_DOUBLES_EQUAL(fpos,-fneg,1.e-12);
      //test gradient;
    for(i=0;i<tetPatch.num_vertices(); ++i){
      for(j=0;j<3;++j){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gpos[i][j], -gneg[i][j], 1.e-12);
      }
    }
    
    valid=LP4.compute_hessian(tetPatch, Hpos, gpos, fpos, err);
    MSQ_CHKERR(err); CPPUNIT_ASSERT(valid);
    valid=LP4_neg.compute_hessian(tetPatch, Hneg, gneg, fneg, err);
    MSQ_CHKERR(err);CPPUNIT_ASSERT(valid);
      //std::cout<<"\nOrig fpos = "<<fpos<<" Mod fneg = "<<fneg;
      //test function value
    CPPUNIT_ASSERT_DOUBLES_EQUAL(fpos,-fneg,1.e-12);
      //test gradient;
    for(i=0;i<tetPatch.num_vertices(); ++i){
      for(j=0;j<3;++j){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gpos[i][j], -gneg[i][j], 1.e-12);
      }
    }
    
      // test hessian
    Matrix3D* block_pos;
    Matrix3D* block_neg;
    CPPUNIT_ASSERT(Hpos.size() == Hneg.size());
      //std::cout<<"\nHessian size "<< Hneg.size();
      //std::cout<<"Hpos"<<Hpos;
      //std::cout<<"Hneg"<<Hneg;
    for (size_t m=0; m<Hpos.size(); ++m) {
      for (size_t n=m; n<3; ++n) {
         if(!(m==0 && n==4)){ 
            block_pos = Hpos.get_block(m,n);
            block_neg = Hneg.get_block(m,n);
            for (i=0; i<3; ++i)
               for (j=0; j<3; ++j)
                  CPPUNIT_ASSERT_DOUBLES_EQUAL( (*block_pos)[i][j],
                                                -((*block_neg)[i][j]),
                                                1.e-12);
         }
      }
    }
      // THEN TRY L1 **************************************************
        // ... and builds an objective function with it
    LPtoPTemplate LP1(mean_ratio, 1, err);
    LPtoPTemplate LP1_neg(mean_ratio_neg, 1, err);
      //test evaluate
    valid=LP1.evaluate(tetPatch, fpos, err);MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);
    valid=LP1_neg.evaluate(tetPatch, fneg, err);MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);
      //std::cout<<"\nFrom eval Orig fpos = "<<fpos<<" Mod fneg = "<<fneg;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(fpos,-fneg,1.e-12);

    valid=LP1.compute_gradient(tetPatch, gpos, fpos, err);
    MSQ_CHKERR(err); CPPUNIT_ASSERT(valid);
    valid=LP1_neg.compute_gradient(tetPatch, gneg, fneg, err);
    MSQ_CHKERR(err);CPPUNIT_ASSERT(valid);
      //std::cout<<"\nOrig fpos = "<<fpos<<" Mod fneg = "<<fneg;
      //test function value
    CPPUNIT_ASSERT_DOUBLES_EQUAL(fpos,-fneg,1.e-12);
      //test gradient;
    for(i=0;i<tetPatch.num_vertices(); ++i){
      for(j=0;j<3;++j){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gpos[i][j], -gneg[i][j], 1.e-12);
      }
    }
    
    valid=LP1.compute_hessian(tetPatch, Hpos, gpos, fpos, err);
    MSQ_CHKERR(err); CPPUNIT_ASSERT(valid);
    valid=LP1_neg.compute_hessian(tetPatch, Hneg, gneg, fneg, err);
    MSQ_CHKERR(err);CPPUNIT_ASSERT(valid);
      //std::cout<<"\nOrig fpos = "<<fpos<<" Mod fneg = "<<fneg;
      //test function value
    CPPUNIT_ASSERT_DOUBLES_EQUAL(fpos,-fneg,1.e-12);
      //test gradient;
    for(i=0;i<tetPatch.num_vertices(); ++i){
      for(j=0;j<3;++j){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gpos[i][j], -gneg[i][j], 1.e-12);
      }
    }
    
      // test hessian
    CPPUNIT_ASSERT(Hpos.size() == Hneg.size());
      //std::cout<<"\nHessian size "<< Hneg.size();
      //std::cout<<"Hpos"<<Hpos;
      //std::cout<<"Hneg"<<Hneg;
    for (size_t m=0; m<Hpos.size(); ++m) {
      for (size_t n=m; n<3; ++n) {
         if(!(m==0 && n==4)){ 
            block_pos = Hpos.get_block(m,n);
            block_neg = Hneg.get_block(m,n);
            for (i=0; i<3; ++i)
               for (j=0; j<3; ++j)
                  CPPUNIT_ASSERT_DOUBLES_EQUAL( (*block_pos)[i][j], -((*block_neg)[i][j]), 1.e-12);
         }
      }
    }
    
    delete[] gpos;
    delete[] gneg;
    delete mean_ratio;
    delete mean_ratio_neg;
    
  }
  

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ObjectiveFunctionTest, "ObjectiveFunctionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ObjectiveFunctionTest, "Unit");
