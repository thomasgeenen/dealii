// ---------------------------------------------------------------------
// $Id: paralution_solver_01.cc 31349 2013-10-20 19:07:06Z maier $
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


// Create a test for ParalutionWrappers::SolverCG

#include "../tests.h"
#include "testmatrix.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/paralution_sparse_matrix.h>
#include <deal.II/lac/paralution_vector.h>
#include <deal.II/lac/paralution_solver.h>

int main()
{
  paralution::init_paralution();

  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  SolverControl control(100,1.e-3);
  ParalutionWrappers::SolverCG solver(control);

  for (unsigned int size=4; size<=30; size*=3)
  {
    unsigned int dim = (size-1)*(size-1);

    deallog << "Size "<< size <<" Unknowns "<< dim << std::endl;

    // Make matrix 
    FDMatrix testproblem(size, size);
    SparsityPattern structure(dim, dim, 5);
    testproblem.five_point_structure(structure);
    structure.compress();
    ParalutionWrappers::SparseMatrix<double> A(structure);
    testproblem.five_point(A);
    A.convert_to_paralution_csr();

    ParalutionWrappers::Vector<double> u(dim);
    ParalutionWrappers::Vector<double> f(dim);
    for (unsigned int i=0; i<dim; ++i)
    {
      u[i] = 0.;
      f[i] = 1.;
    }
    solver.solve(A,u,f);
  }

  paralution::stop_paralution();
}
