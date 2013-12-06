// ---------------------------------------------------------------------
// $Id: paralution_vector_01.cc 31349 2013-10-20 19:07:06Z maier $
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test copy constructor and access elements of ParalutionWrappers::Vector

#include "../tests.h"
#include <deal.II/lac/paralution_vector.h>
#include "paralution.hpp"

void check()
{
  ParalutionWrappers::Vector<double> vector(10);

  for (unsigned int i=0; i<10; ++i)
    vector[i] = i;

  ParalutionWrappers::Vector<double> copy_vector(vector);

  for (unsigned int i=0; i<10; ++i)
    deallog << vector[i] << std::endl;

  for (unsigned int i=0; i<10; ++i)
    deallog << vector(i)-copy_vector(i) <<std::endl;
}

int main()
{
  paralution::init_paralution();

  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check();

  paralution::stop_paralution();
}
