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
#include <deal.II/base/paralution.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/paralution_sparse_matrix.h>

template <typename Number>
void check()
{
  SparsityPattern sparsity_pattern(4,5,2);
  sparsity_pattern.add(0,2);
  sparsity_pattern.add(0,0);
  sparsity_pattern.add(1,0);
  sparsity_pattern.add(1,3);
  sparsity_pattern.add(2,4);
  sparsity_pattern.add(2,2);
  sparsity_pattern.add(3,0);
  sparsity_pattern.add(3,4);
  sparsity_pattern.compress();

  ParalutionWrappers::SparseMatrix<Number> matrix;
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
  Utilities::Paralution::Paralution_InitFinalize paralution(1);

  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<float>();
  check<double>();
}

