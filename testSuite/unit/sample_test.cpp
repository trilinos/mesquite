#ifndef __VECTOR3DTEST_H__
#define __VECTOR3DTEST_H__

#include "Vector3D.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

class Vector3DTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Vector3DTest);
  CPPUNIT_TEST (test_default_constructor);
  CPPUNIT_TEST (test_double_constructor);
  CPPUNIT_TEST (test_copy_constructor);
  CPPUNIT_TEST_EXCEPTION (throw_exception, CppUnit::SignalException);
  CPPUNIT_TEST (throw_exception);
  CPPUNIT_TEST_EXCEPTION (throw_exception, CppUnit::SignalException);
  CPPUNIT_TEST (throw_exception);
  CPPUNIT_TEST (test_default_constructor);
  CPPUNIT_TEST (test_double_constructor);
  CPPUNIT_TEST (test_copy_constructor);
  CPPUNIT_TEST_SUITE_END();

public:
  Vector3DTest()
    {}
  
  void test_default_constructor()
    {
      Mesquite::Vector3D v;
      CPPUNIT_ASSERT(v.x() == 0);
      CPPUNIT_ASSERT(v.y() == 0);
      CPPUNIT_ASSERT(v.z() == 0);
    }
  void test_double_constructor()
    {
      Mesquite::Vector3D v(3, 2, 1);
      CPPUNIT_ASSERT(v.x() == 3);
      CPPUNIT_ASSERT(v.y() == 2);
      CPPUNIT_ASSERT(v.z() == 1);
    }
  void test_copy_constructor()
    {
      Mesquite::Vector3D v2(3, 2, 1);
      Mesquite::Vector3D v(v2);
      
      CPPUNIT_ASSERT(v.x() == 3);
      CPPUNIT_ASSERT(v.y() == 2);
      CPPUNIT_ASSERT(v.z() == 1);
    }
  void test_equals()
    {
      Mesquite::Vector3D v1(1,4,5.2);
      Mesquite::Vector3D v2(1,4,5.2);
      CPPUNIT_ASSERT(v1 == v2);
    }
  void throw_exception()
    {
      Mesquite::Vector3D* v = NULL;
      double d = v->x();
      std::cout << d << std::endl;
    }
};

class Vector3DTest2 : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Vector3DTest);
  CPPUNIT_TEST (test_equals);
  CPPUNIT_TEST_SUITE_END();

public:
  Vector3DTest2()
    {}
  
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(Vector3DTest, "Misc");

#endif
