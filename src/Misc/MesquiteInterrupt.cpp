// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-

#include "MesquiteInterrupt.hpp"

void Mesquite::signal_interrupt()
{ Mesquite::interruptFlag = true; }
