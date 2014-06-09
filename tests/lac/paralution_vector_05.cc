// ---------------------------------------------------------------------
// $Id: paralution_vector_04.cc 31349 2013-10-20 19:07:06Z maier $
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


// Test for ParalutionWrappers::Vector operator*= operator/= operator+=
// operator-=

#include "../tests.h"
#include <deal.II/base/paralution.h>
#include <deal.II/lac/paralution_vector.h>

void check_1()
{
  ParalutionWrappers::Vector<double> vector(10);

  for (unsigned int i=0; i<10; ++i)
    vector[i] = i;
  
  vector *= 10.;
  vector /= 5.;

  for (unsigned int i=0; i<10; ++i)
    deallog << vector[i] << std::endl;
}

void check_2()
{
  ParalutionWrappers::Vector<double> vector_1(10);
  ParalutionWrappers::Vector<double> vector_2(10);
  ParalutionWrappers::Vector<double> vector_3(10);

  for (unsigned int i=0; i<10; ++i)
  {
    vector_1[i] = i;
    vector_2[i] = 10*i;
    vector_3[i] = i;
  }

  vector_1 += vector_2;
  vector_1 -= vector_3;

  for (unsigned int i=0; i<10; ++i)
    deallog << vector_1[i] << std::endl;
}

int main()
{
  Utilities::Paralution::Paralution_InitFinalize paralution(1);

  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check_1();
  check_2();
}
