// ---------------------------------------------------------------------
// $Id: paralution_sparse_matrix_01.cc 31349 2013-10-20 19:07:06Z maier $
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


// Test ParalutionWrappers::SparseMatrix

#include "../tests.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/paralution_sparse_matrix.h>
#include "paralution.hpp"

template <typename Number>
void check()
{
  ParalutionWrappers::SparseMatrix<Number> matrix;
  SparsityPattern pattern(4,5,2);
  pattern.add(0,2);
  pattern.add(0,0);
  pattern.add(1,0);
  pattern.add(1,3);
  pattern.add(2,4);
  pattern.add(2,2);
  pattern.add(3,0);
  pattern.add(3,4);
  pattern.compress();

  matrix.reinit(sparsity_pattern);
  matrix.add(0,2,3.5);
  matrix.add(0,0,1.);
  matrix.add(1,0,-2.);
  matrix.add(1,3,1.5);
  matrix.add(2,4,-2.25);
  matrix.add(2,2,-0.5);
  matrix.add(3,0,2.);
  matrix.add(3,4,0.);

  matrix.convert_to_paralution_csr();

  deallog << matrix.m() <<std::endl;
  deallog << matrix.n() <<std::endl;
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

  check<float>();
  check<double>();

  paralution::stop_paralution();
}

