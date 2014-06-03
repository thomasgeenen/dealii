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


// Test ParalutionWrappers::Vector constructor

#include "../tests.h"
#include <deal.II/lac/paralution_vector.h>
#include "paralution.hpp"

template <typename Number>
void check()
{
  ParalutionWrappers::Vector<Number> vector_1;
  ParalutionWrappers::Vector<Number> vector_2(100);

  deallog << vector_1.size() << std::endl;
  deallog << vector_2.size() << std::endl;
  vector_2.clear();
  deallog << vector_2.size() << std::endl;
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

  check<int>();
  check<float>();
  check<double>();

  paralution::stop_paralution();
}
