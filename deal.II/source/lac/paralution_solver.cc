// ---------------------------------------------------------------------
// $Id: paralution_solver.cc 31349 2013-10-20 19:07:06Z maier $
//
// Copyright (C) 2013, 2014 by the deal.II authors
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

  template <typename Number>
  void SolverBase::execute_solve(std_cxx1x::shared_ptr<paralution::IterativeLinearSolver<paralution::
                                 LocalMatrix<Number>,paralution::LocalVector<Number>,Number> > solver,
                                 const SparseMatrix<Number>     &A,
                                 Vector<Number>                 &x,
                                 const Vector<Number>           &b,
                                 const PreconditionBase<Number> &preconditioner,
                                 bool                            move_to_accelerator)
  {
    // Set the preconditioner.
    solver->SetPreconditioner(*(preconditioner.preconditioner));

    // Set the system to solve
    solver->SetOperator(A.paralution_matrix());

    // Set absolute tolerance, relative tolerance, divergence tolerance,
    // maximum number of iterations.
    solver->Init(solver_control.tolerance(),0.,1.e100,
                 solver_control.max_steps());

    // Move the solver to the accelerator if necessary.
    if (move_to_accelerator==true)
      solver->MoveToAccelerator();

    solver->Build();
    solver->Solve(b.paralution_vector(),&(x.paralution_vector()));
  }



  /* ---------------------- SolverCG ------------------------ */

  SolverCG::SolverCG (SolverControl &cn)
    :
    SolverBase (cn)
  {}



  template <typename Number>
  void SolverCG::solve(const SparseMatrix<Number>     &A,
                       Vector<Number>                 &x,
                       const Vector<Number>           &b,
                       const PreconditionBase<Number> &preconditioner,
                       bool                            move_to_accelerator)
  {
    std_cxx1x::shared_ptr<paralution::CG<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> > solver(new  paralution::
                                                               CG<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,Number>);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* -------------------- SolverBicgstab --------------------- */

  SolverBicgstab::SolverBicgstab (SolverControl &cn)
    :
    SolverBase (cn)
  {}



  template <typename Number>
  void SolverBicgstab::solve(const SparseMatrix<Number>     &A,
                             Vector<Number>                 &x,
                             const Vector<Number>           &b,
                             const PreconditionBase<Number> &preconditioner,
                             bool                            move_to_accelerator)
  {
    std_cxx1x::shared_ptr<paralution::BiCGStab<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> > solver(new  paralution::
                                                               BiCGStab<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,Number>);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* --------------------- SolverGMRES ----------------------- */

  SolverGMRES::AdditionalData::AdditionalData(const unsigned int restart_parameter)
    :
    restart_parameter(restart_parameter)
  {}




  SolverGMRES::SolverGMRES (SolverControl        &cn,
                            const AdditionalData &data)
    :
    SolverBase (cn),
    additional_data (data)
  {}



  template <typename Number>
  void SolverGMRES::solve(const SparseMatrix<Number>     &A,
                          Vector<Number>                 &x,
                          const Vector<Number>           &b,
                          const PreconditionBase<Number> &preconditioner,
                          bool                            move_to_accelerator)
  {
    std_cxx1x::shared_ptr<paralution::GMRES<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> > solver(new  paralution::
                                                               GMRES<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,Number>);

    // Set the restart parameter.
    solver->SetBasisSize(additional_data.restart_parameter);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }
}



// Explicit instantiations
namespace ParalutionWrappers
{
  template void SolverCG::solve<float>(const SparseMatrix<float>     &A,
                                       Vector<float>                 &x,
                                       const Vector<float>           &b,
                                       const PreconditionBase<float> &preconditioner,
                                       bool                           move_to_accelerator);

  template void SolverCG::solve<double>(const SparseMatrix<double>     &A,
                                        Vector<double>                 &x,
                                        const Vector<double>           &b,
                                        const PreconditionBase<double> &preconditioner,
                                        bool                            move_to_accelerator);

  template void SolverBicgstab::solve<float>(const SparseMatrix<float>     &A,
                                             Vector<float>                 &x,
                                             const Vector<float>           &b,
                                             const PreconditionBase<float> &preconditioner,
                                             bool                           move_to_accelerator);

  template void SolverBicgstab::solve<double>(const SparseMatrix<double>     &A,
                                              Vector<double>                 &x,
                                              const Vector<double>           &b,
                                              const PreconditionBase<double> &preconditioner,
                                              bool                            move_to_accelerator);

  template void SolverGMRES::solve<float>(const SparseMatrix<float>     &A,
                                          Vector<float>                 &x,
                                          const Vector<float>           &b,
                                          const PreconditionBase<float> &preconditioner,
                                          bool                           move_to_accelerator);

  template void SolverGMRES::solve<double>(const SparseMatrix<double>     &A,
                                           Vector<double>                 &x,
                                           const Vector<double>           &b,
                                           const PreconditionBase<double> &preconditioner,
                                           bool                            move_to_accelerator);
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION
