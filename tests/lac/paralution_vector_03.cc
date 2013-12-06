// ---------------------------------------------------------------------
// $Id: paralution_vector_03.cc 31349 2013-10-20 19:07:06Z maier $
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


// Test ParalutionWrappers::Vector iterator

#include "../tests.h"
#include <deal.II/lac/paralution_vector.h>
#include "paralution.hpp"

void check()
{
  ParalutionWrappers::Vector<double> vector(10);

  ParalutionWrappers::Vector<double>::iterator it(vector.begin());
  ParalutionWrappers::Vector<double>::iterator end_it(vector.end());

  for (double i=0.; it!=end_it; ++it, i+=1.)
    *it = i;

  ParalutionWrappers::Vector<double>::iterator cit(vector.begin());
  ParalutionWrappers::Vector<double>::iterator cend_it(vector.end());

  for (; cit!=cend_it; ++cit)
    deallog << *cit << std::endl;
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
