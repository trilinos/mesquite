// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-

#include "Mesquite.hpp"

// When changing version number/type, change these #defines, as well
// as what is returned from Mesquite::release_type().
#define MSQ_MAJOR_VERSION 1
#define MSQ_MINOR_VERSION 0
#define MSQ_BUILD_NUMBER 0
#define MSQ_VERSION_STRING "1.0 Alpha"
#define MSQ_BUILD_STRING "Build Number 0"

const char* Mesquite::version_string(bool include_build_number)
{
  if (include_build_number)
    return MSQ_VERSION_STRING MSQ_BUILD_STRING;
  return MSQ_VERSION_STRING;
}

unsigned int Mesquite::major_version_number()
{
  return MSQ_MAJOR_VERSION;
}

unsigned int Mesquite::minor_version_number()
{
  return MSQ_MINOR_VERSION;
}

unsigned int Mesquite::build_number()
{
  return MSQ_BUILD_NUMBER;
}

Mesquite::ReleaseType Mesquite::release_type()
{  
  return Mesquite::ALPHA;
}

