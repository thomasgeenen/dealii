// ---------------------------------------------------------------------
// $Id: paralution_solver.cc 31349 2013-10-20 19:07:06Z maier $
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

#include <deal.II/lac/paralution_solver.h>

#ifdef DEAL_II_WITH_PARALUTION

#include "paralution.hpp"

DEAL_II_NAMESPACE_OPEN

namespace ParalutionWrappers
{
  SolverBase::SolverBase (SolverControl &cn)
    :
    solver_control (cn)
  {}



  SolverControl &SolverBase::control() const
  {
    return solver_control;
  }



  /* ---------------------- SolverCG ------------------------ */

  SolverCG::SolverCG (SolverControl &cn)
    :
    SolverBase (cn)
  {}



  template <typename Number>
  void SolverCG::solve(const SparseMatrix<Number> &A,
                       Vector<Number>             &x,
                       const Vector<Number>       &b)
  {
    paralution::CG<paralution::LocalMatrix<Number>,
               paralution::LocalVector<Number>,Number> solver;
    solver.SetOperator(A.paralution_matrix());
    // Set absolute tolerance, relative tolerance, divergence tolerance,
    // maximum number of iterations.
    solver.Init(solver_control.tolerance(),1e10,1e10,
                solver_control.max_steps());
    solver.Build();
    solver.Solve(b.paralution_vector(),&(x.paralution_vector()));
  }



  /* -------------------- SolverBicgstab --------------------- */

  SolverBicgstab::SolverBicgstab (SolverControl &cn)
    :
    SolverBase (cn)
  {}



  template <typename Number>
  void SolverBicgstab::solve(const SparseMatrix<Number> &A,
                             Vector<Number>             &x,
                             const Vector<Number>       &b)
  {
    paralution::BiCGStab<paralution::LocalMatrix<Number>,
               paralution::LocalVector<Number>,Number> solver;
    // Set absolute tolerance, relative tolerance, divergence tolerance,
    // maximum number of iterations.
    solver.Init(solver_control.tolerance(),1e10,1e10,
                solver_control.max_steps());
    solver.Build();
    solver.Solve(b.paralution_vector(),&(x.paralution_vector()));
  }



  /* --------------------- SolverGMRES ----------------------- */

  SolverGMRES::SolverGMRES (SolverControl        &cn,
                            const AdditionalData &data)
    :
    SolverBase (cn),
    additional_data (data)
  {}



  template <typename Number>
  void SolverGMRES::solve(const SparseMatrix<Number> &A,
                          Vector<Number>             &x,
                          const Vector<Number>       &b)
  {
    paralution::GMRES<paralution::LocalMatrix<Number>,
               paralution::LocalVector<Number>,Number> solver;
    // Set absolute tolerance, relative tolerance, divergence tolerance,
    // maximum number of iterations.
    solver.Init(solver_control.tolerance(),1e10,1e10,
                solver_control.max_steps());
    solver.SetBasisSize(additional_data.restart_parameter);
    solver.Build();
    solver.Solve(b.paralution_vector(),&(x.paralution_vector()));
  }


}



// Explicit instantiations
namespace ParalutionWrappers
{
  template void SolverCG::solve<float>(const SparseMatrix<float> &A,
                                       Vector<float>             &x,
                                       const Vector<float>       &b);

  template void SolverCG::solve<double>(const SparseMatrix<double> &A,
                                        Vector<double>             &x,
                                        const Vector<double>       &b);

  template void SolverBicgstab::solve<float>(const SparseMatrix<float> &A,
                                             Vector<float>             &x,
                                             const Vector<float>       &b);

  template void SolverBicgstab::solve<double>(const SparseMatrix<double> &A,
                                              Vector<double>             &x,
                                              const Vector<double>       &b);

  template void SolverGMRES::solve<float>(const SparseMatrix<float> &A,
                                          Vector<float>             &x,
                                          const Vector<float>       &b);

  template void SolverGMRES::solve<double>(const SparseMatrix<double> &A,
                                           Vector<double>             &x,
                                           const Vector<double>       &b);
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION
