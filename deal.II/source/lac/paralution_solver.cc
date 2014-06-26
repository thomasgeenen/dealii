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



  SolverControl& SolverBase::control() const
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



  /* --------------------- SolverRichardson----------------- */

  SolverRichardson::AdditionalData::AdditionalData(const unsigned int verbose,
                                                   const double relaxation_parameter)
    :
    verbose(verbose),
    relaxation_parameter(relaxation_parameter)
  {}



  SolverRichardson::SolverRichardson (SolverControl        &cn,
                                      const AdditionalData &data)
    :
    SolverBase (cn),
    additional_data (data)
  {}



  template <typename Number>
  void SolverRichardson::solve(const SparseMatrix<Number>     &A,
                               Vector<Number>                 &x,
                               const Vector<Number>           &b,
                               const PreconditionBase<Number> &preconditioner,
                               bool                            move_to_accelerator)
  {
    std_cxx1x::shared_ptr<paralution::FixedPoint<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> > solver(new  paralution::
                                                               FixedPoint<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,Number>);

    // Set the verbosity of the solver.
    solver->Verbose(additional_data.verbose);

    // Set the relaxation parameter.
    solver->SetRelaxation(additional_data.relaxation_parameter);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* ---------------------- SolverCG ------------------------- */

  SolverCG::AdditionalData::AdditionalData(const unsigned int verbose)
    :
    verbose(verbose)
  {}



  SolverCG::SolverCG (SolverControl &cn, const AdditionalData &data)
    :
    SolverBase (cn),
    additional_data (data)
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

    // Set the verbosity of the solver.
    solver->Verbose(additional_data.verbose);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* ---------------------- SolverCR ------------------------- */

  SolverCR::AdditionalData::AdditionalData(const unsigned int verbose)
    :
    verbose(verbose)
  {}



  SolverCR::SolverCR (SolverControl &cn, const AdditionalData &data)
    :
    SolverBase (cn),
    additional_data (data)
  {}



  template <typename Number>
  void SolverCR::solve(const SparseMatrix<Number>     &A,
                       Vector<Number>                 &x,
                       const Vector<Number>           &b,
                       const PreconditionBase<Number> &preconditioner,
                       bool                            move_to_accelerator)
  {
    std_cxx1x::shared_ptr<paralution::CR<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> > solver(new  paralution::
                                                               CR<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,Number>);

    // Set the verbosity of the solver.
    solver->Verbose(additional_data.verbose);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* ---------------------- SolverDPCG ----------------------- */

  SolverDPCG::AdditionalData::AdditionalData(const unsigned int verbose,
                                             const const unsigned int n_deflated_vectors)
    :
    verbose (verbose),
    n_deflated_vectors(n_deflated_vectors)
  {}



  SolverDPCG::SolverDPCG (SolverControl &cn, const AdditionalData &data)
    :
    SolverBase (cn),
    additional_data (data)
  {}



  template <typename Number>
  void SolverDPCG::solve(const SparseMatrix<Number>     &A,
                         Vector<Number>                 &x,
                         const Vector<Number>           &b,
                         const PreconditionBase<Number> &preconditioner,
                         bool                            move_to_accelerator)
  {
    std_cxx1x::shared_ptr<paralution::DPCG<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> > solver(new  paralution::
                                                               DPCG<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,Number>);

    // Set the verbosity of the solver.
    solver->Verbose(additional_data.verbose);

    // Set the number of deflated vectors.
    solver->SetNVectors(additional_data.n_deflated_vectors);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* -------------------- SolverBicgstab --------------------- */

  SolverBicgstab::AdditionalData::AdditionalData(const unsigned int verbose)
    :
    verbose(verbose)
  {}



  SolverBicgstab::SolverBicgstab (SolverControl &cn, const AdditionalData &data)
    :
    SolverBase (cn),
    additional_data (data)
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

    // Set the verbosity of the solver.
    solver->Verbose(additional_data.verbose);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* --------------------- SolverGMRES ----------------------- */

  SolverGMRES::AdditionalData::AdditionalData(const unsigned int verbose,
                                              const unsigned int restart_parameter)
    :
    verbose(verbose),
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

    // Set the verbosity of the solver.
    solver->Verbose(additional_data.verbose);

    // Set the restart parameter.
    solver->SetBasisSize(additional_data.restart_parameter);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* --------------------- SolverFGMRES ---------------------- */

  SolverFGMRES::AdditionalData::AdditionalData(const unsigned int verbose,
                                               const unsigned int restart_parameter)
    :
    verbose(verbose),
    restart_parameter(restart_parameter)
  {}



  SolverFGMRES::SolverFGMRES (SolverControl        &cn,
                              const AdditionalData &data)
    :
    SolverBase (cn),
    additional_data (data)
  {}



  template <typename Number>
  void SolverFGMRES::solve(const SparseMatrix<Number>     &A,
                           Vector<Number>                 &x,
                           const Vector<Number>           &b,
                           const PreconditionBase<Number> &preconditioner,
                           bool                            move_to_accelerator)
  {
    std_cxx1x::shared_ptr<paralution::FGMRES<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> > solver(new  paralution::
                                                               FGMRES<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,Number>);

    // Set the verbosity of the solver.
    solver->Verbose(additional_data.verbose);

    // Set the restart parameter.
    solver->SetBasisSize(additional_data.restart_parameter);

    this->execute_solve<Number>(solver,A,x,b,preconditioner,move_to_accelerator);
  }



  /* ---------------------- SolverAMG------------------------- */

  SolverAMG::AdditionalData::AdditionalData(const unsigned int verbose,
                                            const mg_solver coarse_solver,
                                            const unsigned int n_unknowns_coarse_level,
                                            const mg_solver smoother,
                                            const mg_preconditioner preconditioner,
                                            const mg_cycle cycle,
                                            const mg_interpolation interpolation,
                                            const double coupling_strength,
                                            const double relaxation_parameter,
                                            const double over_interpolation)
    :
    verbose(verbose),
    coarse_solver(coarse_solver),
    n_unknowns_coarse_level(n_unknowns_coarse_level),
    smoother(smoother),
    preconditioner(preconditioner),
    cycle(cycle),
    interpolation(interpolation),
    coupling_strength(coupling_strength),
    relaxation_parameter(relaxation_parameter),
    over_interpolation(over_interpolation)
  {}



  SolverAMG::SolverAMG (SolverControl        &cn,
                        const AdditionalData &data)
    :
    solver_control (cn),
    additional_data (data)
  {}



  template <typename Number>
  void SolverAMG::solve(SparseMatrix<Number>           &A,
                        Vector<Number>                 &x,
                        Vector<Number>                 &b,
                        bool                            move_to_accelerator)
  {
    paralution::AMG<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,Number> solver;

    // Set the system to solve.
    solver.SetOperator(A.paralution_matrix());
    solver.SetOperatorFormat(CSR);
    
    // Set the verbosity of the multigrid.
    solver.Verbose(additional_data.verbose);
    
    // Set number of unknowns on coarsest level.
    solver.SetCoarsestLevel(additional_data.n_unknowns_coarse_level);

    // Set the interpolation type.
    switch (additional_data.interpolation)
    {
      case aggregation :
        {
          solver.SetInterpolation(paralution::Aggregation);
          solver.SetOverInterp(additional_data.over_interpolation);
          break;
        }
      case smoothed_aggregation :
        {
          solver.SetInterpolation(paralution::SmoothedAggregation);
          solver.SetInterpRelax(additional_data.relaxation_parameter);
          break;
        }
      default :
        {
          AssertThrow(false,ExcMessage("Unknown interpolation type for PreconditionAMG."));
        }
    }

    // Set the coupling strength.
    solver.SetCouplingStrength(additional_data.coupling_strength);

    // Set the type of cycle
    switch (additional_data.cycle)
    {
      case V_cycle :
        {
          solver.SetCycle(paralution::Vcycle);
          break;
        }
      case  W_cycle :
        {
          solver.SetCycle(paralution::Wcycle);
          break;
        }
      case K_cycle :
        {
          solver.SetCycle(paralution::Kcycle);
          break;
        }
      case F_cycle :
        {
          solver.SetCycle(paralution::Fcycle);
        }
      default :
        {
          AssertThrow(false,ExcMessage("Unknown cycle type for PreconditionAMG."));
        }
    }
    
    // Manual smoothers
    solver.SetManualSmoothers(true);
    
    // Manual course grid solver
    solver.SetManualSolver(true);
    
    // Set grid transfer scaling
    solver.SetScaling(true);

    // Build the grid hierarchy.
    solver.BuildHierarchy();

    unsigned int n_levels = solver.GetNumLevels();

    // Smoother for each level
    paralution::IterativeLinearSolver<paralution::LocalMatrix<Number>,
      paralution::LocalVector<Number>,Number> **smoothers = NULL;
    smoothers = new paralution::IterativeLinearSolver<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number>* [n_levels-1];

    // Preconditioner for each smoother
    paralution::Preconditioner<paralution::LocalMatrix<Number>,
      paralution::LocalVector<Number>,Number> **preconds = NULL;
    preconds = new paralution::Preconditioner<paralution::LocalMatrix<Number>,
             paralution::LocalVector<Number>,Number>* [n_levels-1];

    // Coarse Grid Solver
    paralution::IterativeLinearSolver<paralution::LocalMatrix<Number>,
      paralution::LocalVector<Number>,Number> *coarse_solver;

    // Smoother
    switch (additional_data.smoother)
    {
      case richardson :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::FixedPoint<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> *richardson = NULL;
            richardson = new paralution::FixedPoint<paralution::LocalMatrix<Number>,
                       paralution::LocalVector<Number>,Number>;
            smoothers[i] = richardson;
            smoothers[i]->Verbose(additional_data.verbose);
          }

          break;
        }
      case cg :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::CG<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *cg = NULL;
            cg = new paralution::CG<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
               Number>;
            smoothers[i] = cg;
            smoothers[i]->Verbose(additional_data.verbose);
          }

          break;
        }
      case cr :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::CR<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *cr = NULL;
            cr = new paralution::CR<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
               Number>;
            smoothers[i] = cr;
            smoothers[i]->Verbose(additional_data.verbose);
          }

          break;
        }
      case bicgstab :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::BiCGStab<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *bicgstab = NULL;
            bicgstab = new paralution::BiCGStab<paralution::LocalMatrix<Number>,
                     paralution::LocalVector<Number>,Number>;
            smoothers[i] = bicgstab;
            smoothers[i]->Verbose(additional_data.verbose);
          }

          break;
        }
      case gmres :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::GMRES<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *gmres = NULL;
            gmres = new paralution::GMRES<paralution::LocalMatrix<Number>,
                  paralution::LocalVector<Number>,Number>;
            smoothers[i] = gmres;
            smoothers[i]->Verbose(additional_data.verbose);
          }

          break;
        }
      case fgmres :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::FGMRES<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *fgmres = NULL;
            fgmres = new paralution::FGMRES<paralution::LocalMatrix<Number>,
                  paralution::LocalVector<Number>,Number>;
            smoothers[i] = fgmres;
            smoothers[i]->Verbose(additional_data.verbose);
          }

          break;
        }
      default :
        {
          AssertThrow(false,ExcMessage("Unknown smoother for AMG."));
        }
    }

    // Preconditioner
    switch (additional_data.preconditioner)
    {
      case jacobi :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::Jacobi<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *jacobi = NULL;
            jacobi = new paralution::Jacobi<paralution::LocalMatrix<Number>,
                   paralution::LocalVector<Number>,Number>;
            preconds[i] = jacobi;
            smoothers[i]->SetPreconditioner(*preconds[i]);
          }

          break;
        }
      case sgs :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::SGS<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *sgs = NULL;
            sgs = new paralution::SGS<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
                Number>;
            preconds[i] = sgs;
            smoothers[i]->SetPreconditioner(*preconds[i]);
          }

          break;
        }
      case multicolored_sgs :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::MultiColoredSGS<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *multicolored_sgs = NULL;
            multicolored_sgs = new paralution::MultiColoredSGS<paralution::LocalMatrix<Number>,
                             paralution::LocalVector<Number>,Number>;
            multicolored_sgs->SetPrecondMatrixFormat(ELL);
            preconds[i] = multicolored_sgs;
            smoothers[i]->SetPreconditioner(*preconds[i]);
          }

          break;
        }
      case multicolored_sor :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::MultiColoredGS<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *multicolored_sor = NULL;
            multicolored_sor = new paralution::MultiColoredGS<paralution::LocalMatrix<Number>,
                     paralution::LocalVector<Number>,Number>;
            multicolored_sor->SetPrecondMatrixFormat(ELL);
            preconds[i] = multicolored_sor;
            smoothers[i]->SetPreconditioner(*preconds[i]);
          }

          break;
        }
      case ilu :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::ILU<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *ilu = NULL;
            ilu = new paralution::ILU<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
                Number>;
            preconds[i] = ilu;
            smoothers[i]->SetPreconditioner(*preconds[i]);
          }

          break;
        }
      case ilut :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::ILUT<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *ilut = NULL;
            ilut = new paralution::ILUT<paralution::LocalMatrix<Number>,
                  paralution::LocalVector<Number>,Number>;
            preconds[i] = ilut;
            smoothers[i]->SetPreconditioner(*preconds[i]);
          }

          break;
        }
      case multicolored_ilu :
        {
          for (unsigned int i=0; i<n_levels-1; ++i)
          {
            paralution::MultiColoredILU<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
              Number> *multicolored_ilu = NULL;
            multicolored_ilu = new paralution::MultiColoredILU<paralution::LocalMatrix<Number>,
                     paralution::LocalVector<Number>,Number>;
            multicolored_ilu->SetPrecondMatrixFormat(ELL);
            preconds[i] = multicolored_ilu;
            smoothers[i]->SetPreconditioner(*preconds[i]);
          }

          break;
        }
      default :
        {
          AssertThrow(false,ExcMessage("Unknown preconditioner for the smoother for AMG."));
        }
    }

    // Set the smoothers.
    solver.SetSmoother(smoothers);
    solver.SetSmootherPreIter(1);
    solver.SetSmootherPostIter(2);

    // Coarse solver
    switch (additional_data.coarse_solver)
    {
      case richardson :
        {
          paralution::FixedPoint<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,
            Number> *cs = new paralution::FixedPoint<paralution::LocalMatrix<Number>,
            paralution::LocalVector<Number>, Number>;
          coarse_solver = cs;
          break;
        }
      case cg :
        {
          paralution::CG<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,
            Number> *cs = new paralution::CG<paralution::LocalMatrix<Number>,
            paralution::LocalVector<Number>, Number>;
          coarse_solver = cs;
          break;
        }
      case cr :
        {
          paralution::CR<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,
            Number> *cs = new paralution::CR<paralution::LocalMatrix<Number>,
            paralution::LocalVector<Number>, Number>;
          coarse_solver = cs;
          break;
        }
      case bicgstab :
        {
          paralution::BiCGStab<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,
            Number> *cs =  new paralution::BiCGStab<paralution::LocalMatrix<Number>,
            paralution::LocalVector<Number>, Number>;
          coarse_solver = cs;
          break;
        }
      case gmres :
        {
          paralution::GMRES<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,
            Number> *cs = new paralution::GMRES<paralution::LocalMatrix<Number>,
            paralution::LocalVector<Number>, Number>;
          coarse_solver = cs;
          break;
        }
      case fgmres :
        {
          paralution::FGMRES<paralution::LocalMatrix<Number>, paralution::LocalVector<Number>,
            Number> *cs = new paralution::FGMRES<paralution::LocalMatrix<Number>,
            paralution::LocalVector<Number>, Number>;
          coarse_solver = cs;
          break;
        }
      default :
        {
          AssertThrow(false,ExcMessage("Unknown coarse solver for PreconditionAMG."));
        }
    }
    solver.SetSolver(*coarse_solver);

    // Set absolute tolerance, relative tolerance, divergence tolerance,
    // maximum number of iterations.
    solver.Init(solver_control.tolerance(),0.,1.e100,solver_control.max_steps());

    // Build the solver.
    solver.Build();

    // Move the solver to the accelerator if necessary.
    if (move_to_accelerator==true)
      {
        A.move_to_accelerator();
        x.move_to_accelerator();
        b.move_to_accelerator();
        solver.MoveToAccelerator();
      }
    solver.Solve(b.paralution_vector(),&(x.paralution_vector()));

    // Let the deal.II SolverControl object know what has happened. If the solve
    // succeeded, the status of the solver control will turn into
    // SolverControl::success. solver_control only exists on the host so the solver
    // needs to be on the most. 
    solver.MoveToHost();

    solver_control.check (solver.GetIterationCount(), solver.GetCurrentResidual());

    if (solver_control.last_check() != SolverControl::success)
      AssertThrow(false, SolverControl::NoConvergence (solver_control.last_step(),
                                                       solver_control.last_value()));

    // Free the coarse solver, the preconditioners and the smoothers.
    if (coarse_solver!=NULL)
    {
      delete coarse_solver;
      coarse_solver = NULL;
    }
    if (preconds!=NULL)
    {
      for (unsigned int i=0; i<n_levels-1; ++i)
        delete [] preconds[i];
      delete [] preconds;
      preconds = NULL;
    }
    if (smoothers!=NULL)
    {
      for (unsigned int i=0; i<n_levels-1; ++i)
        delete [] smoothers[i];
      delete [] smoothers;
      smoothers = NULL;
    }
  }



  SolverControl& SolverAMG::control() const
  {
    return solver_control;
  }
}



// Explicit instantiations
namespace ParalutionWrappers
{
  template void SolverCG::solve<float>(SparseMatrix<float>     &A,
                                       Vector<float>           &x,
                                       Vector<float>           &b,
                                       PreconditionBase<float> &preconditioner,
                                       bool                     move_to_accelerator);

  template void SolverCG::solve<double>(SparseMatrix<double>     &A,
                                        Vector<double>           &x,
                                        Vector<double>           &b,
                                        PreconditionBase<double> &preconditioner,
                                        bool                      move_to_accelerator);

  template void SolverCR::solve<float>(SparseMatrix<float>     &A,
                                       Vector<float>           &x,
                                       Vector<float>           &b,
                                       PreconditionBase<float> &preconditioner,
                                       bool                     move_to_accelerator);

  template void SolverCR::solve<double>(SparseMatrix<double>     &A,
                                        Vector<double>           &x,
                                        Vector<double>           &b,
                                        PreconditionBase<double> &preconditioner,
                                        bool                      move_to_accelerator);

  template void SolverDPCG::solve<float>(SparseMatrix<float>     &A,
                                         Vector<float>           &x,
                                         Vector<float>           &b,
                                         PreconditionBase<float> &preconditioner,
                                         bool                     move_to_accelerator);

  template void SolverDPCG::solve<double>(SparseMatrix<double>     &A,
                                          Vector<double>           &x,
                                          Vector<double>           &b,
                                          PreconditionBase<double> &preconditioner,
                                          bool                      move_to_accelerator);

  template void SolverBicgstab::solve<float>(SparseMatrix<float>     &A,
                                             Vector<float>           &x,
                                             Vector<float>           &b,
                                             PreconditionBase<float> &preconditioner,
                                             bool                     move_to_accelerator);

  template void SolverBicgstab::solve<double>(SparseMatrix<double>     &A,
                                              Vector<double>           &x,
                                              Vector<double>           &b,
                                              PreconditionBase<double> &preconditioner,
                                              bool                      move_to_accelerator);

  template void SolverGMRES::solve<float>(SparseMatrix<float>     &A,
                                          Vector<float>           &x,
                                          Vector<float>           &b,
                                          PreconditionBase<float> &preconditioner,
                                          bool                     move_to_accelerator);

  template void SolverGMRES::solve<double>(SparseMatrix<double>     &A,
                                           Vector<double>           &x,
                                           Vector<double>           &b,
                                           PreconditionBase<double> &preconditioner,
                                           bool                      move_to_accelerator);

  template void SolverFGMRES::solve<float>(SparseMatrix<float>     &A,
                                           Vector<float>           &x,
                                           Vector<float>           &b,
                                           PreconditionBase<float> &preconditioner,
                                           bool                     move_to_accelerator);

  template void SolverFGMRES::solve<double>(SparseMatrix<double>     &A,
                                            Vector<double>           &x,
                                            Vector<double>           &b,
                                            PreconditionBase<double> &preconditioner,
                                            bool                      move_to_accelerator);

  template void SolverAMG::solve<float>(SparseMatrix<float> &A,
                                        Vector<float>       &x,
                                        Vector<float>       &b,
                                        bool                 move_to_accelerator);

  template void SolverAMG::solve<double>(SparseMatrix<double> &A,
                                         Vector<double>       &x,
                                         Vector<double>       &b,
                                         bool                  move_to_accelerator);
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION
