// ---------------------------------------------------------------------
// $Id: paralution_solver.h 30040 2013-07-18 17:06:48Z maier $
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

#ifndef __deal2__paralution_solver_h
#define __deal2__paralution_solver_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PARALUTION

#include <deal.II/base/std_cxx1x/shared_ptr.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/paralution_precondition.h>
#include <deal.II/lac/paralution_vector.h>
#include <deal.II/lac/paralution_sparse_matrix.h>

#include <paralution.hpp>

DEAL_II_NAMESPACE_OPEN

namespace ParalutionWrappers
{
  /**
   * Base class for solver classes using the Paralution solvers except AMG which
   * is not derived from this class.
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin 2013
   */
  class SolverBase
  {
  public:
    /**
     * Constructor. Takes the solver control object and creates the solver.
     */
    SolverBase (SolverControl &cn);

    /**
     * Access to object that controls convergence.
     */
    SolverControl &control() const;

  protected:
    /**
     * Initialize the solver and solve the system of equations.
     */
    template <typename Number>
    void execute_solve(std_cxx1x::shared_ptr<paralution::IterativeLinearSolver<paralution::
                       LocalMatrix<Number>,paralution::LocalVector<Number>,Number> > solver,
                       const SparseMatrix<Number>     &A,
                       Vector<Number>                 &x,
                       const Vector<Number>           &b,
                       const PreconditionBase<Number> &preconditioner,
                       bool                            move_to_accelerator);

    /**
     * Reference to the object that controls convergence of the iterative
     * solver. In fact, for these Paralution wrappers, Paralution does so
     * itself, but we copy the data from this object before starting solution
     * process, and copy the data back into it afterwards.
     */
    SolverControl &solver_control;
  };



  /**
   * An implementation of the solver interface using the Paralution Richardson solver
   * (paralution::fixed_point_iteration).
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin, 2014
   */
  class SolverRichardson : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the verbosity of the solver to zero and
       * the relaxation parameter to one.
       */
      AdditionalData (const unsigned int verbose = 0,
                      const double relaxation_parameter = 1.);

      /**
       * Verbosity of the solver: 0 =  no output, 1 = print information at the
       * start and at the end, 3 = print residual at each iteration.
       */
      unsigned int verbose;

      /**
       * Relaxation parameter.
       */
      double relaxation_parameter;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particular solver.
     */
    SolverRichardson (SolverControl        &cn,
                      const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the CG solver of
     * Paralution. If the flag @p move_to_accelerator is set to true, the
     * solver is built on the accelerator.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number>     &A,
                Vector<Number>                 &x,
                const Vector<Number>           &b,
                const PreconditionBase<Number> &preconditioner,
                bool                            move_to_accelerator=false);

  private:
    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Paralution CG solver.
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin, 2013
   */
  class SolverCG : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the verbosity of the solver to zero.
       */
      AdditionalData (const unsigned int verbose = 0);

      /**
       * Verbosity of the solver: 0 =  no output, 1 = print information at the
       * start and at the end, 3 = print residual at each iteration.
       */
      unsigned int verbose;

      /**
       * Relaxation parameter.
       */
      double relaxation_parameter;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particular solver.
     */
    SolverCG (SolverControl        &cn,
              const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the CG solver of
     * Paralution. If the flag @p move_to_accelerator is set to true, the
     * solver is built on the accelerator.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number>     &A,
                Vector<Number>                 &x,
                const Vector<Number>           &b,
                const PreconditionBase<Number> &preconditioner,
                bool                            move_to_accelerator=false);

  private:
    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Paralution CR solver.
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin, 2014
   */
  class SolverCR : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the verbosity of the solver to zero.
       */
      AdditionalData (const unsigned int verbose = 0);

      /**
       * Verbosity of the solver: 0 =  no output, 1 = print information at the
       * start and at the end, 3 = print residual at each iteration.
       */
      unsigned int verbose;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particular solver.
     */
    SolverCR (SolverControl        &cn,
              const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the CR solver of
     * Paralution. If the flag @p move_to_accelerator is set to true, the
     * solver is built on the accelerator.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number>     &A,
                Vector<Number>                 &x,
                const Vector<Number>           &b,
                const PreconditionBase<Number> &preconditioner,
                bool                            move_to_accelerator=false);

  private:
    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Paralution DPCG solver.
   * DPCG stands for Deflated Preconditioned Conjugate Gradient. It is a
   * two-level preconditioned CG algorithm to iteratively solve and
   * ill-conditioned linear system. Deflation tries to remove the bad
   * eigenvalues from the spectrum of the preconditioned system.
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin, 2014
   */
  class SolverDPCG : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the verbosity of the solver to zero and
       * set the number of deflated vectors to two.
       */
      AdditionalData (const unsigned int verbose = 0,
                      const unsigned int n_deflated_vectors = 2);

      /**
       * Verbosity of the solver: 0 =  no output, 1 = print information at the
       * start and at the end, 3 = print residual at each iteration.
       */
      unsigned int verbose;

      /**
       * Number of deflated vectors.
       */
      unsigned int n_deflated_vectors;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particular solver.
     */
    SolverDPCG (SolverControl        &cn,
                const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the DPCG solver of
     * Paralution. If the flag @p move_to_accelerator is set to true, the
     * solver is built on the accelerator.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number>     &A,
                Vector<Number>                 &x,
                const Vector<Number>           &b,
                const PreconditionBase<Number> &preconditioner,
                bool                            move_to_accelerator=false);

  private:
    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Paralution BiCGStab
   * solver.
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin, 2013
   */
  class SolverBicgstab : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the verbosity of the solver to zero.
       */
      AdditionalData (const unsigned int verbose = 0);

      /**
       * Verbosity of the solver: 0 =  no output, 1 = print information at the
       * start and at the end, 3 = print residual at each iteration.
       */
      unsigned int verbose;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particular solver.
     */
    SolverBicgstab (SolverControl &cn,
                    const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the BiCGStab solver of
     * Paralution. If the flag @p move_to_accelerator is set to true, the
     * solver is built on the accelerator.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number>     &A,
                Vector<Number>                 &x,
                const Vector<Number>           &b,
                const PreconditionBase<Number> &preconditioner,
                bool                            move_to_accelerator=false);

  private:
    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Paralution GMRES
   * solver.
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin, 2013
   */
  class SolverGMRES : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the verbosity of the solver to zero and
       * set the size of the Krylov space to 30, i.e. do a restart every 30
       * iterations.
       */
      AdditionalData (const unsigned int verbose = 0,
                      const unsigned int restart_parameter = 30);

      /**
       * Verbosity of the solver: 0 =  no output, 1 = print information at the
       * start and at the end, 3 = print residual at each iteration.
       */
      unsigned int verbose;

      /**
       * Size of the Krylov space.
       */
      unsigned int restart_parameter;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particular solver.
     */
    SolverGMRES (SolverControl        &cn,
                 const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the GMRES solver of
     * Paralution. If the flag @p move_to_accelerator is set to true, the
     * solver is built on the accelerator.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number>     &A,
                Vector<Number>                 &x,
                const Vector<Number>           &b,
                const PreconditionBase<Number> &preconditioner,
                bool                            move_to_accelerator=false);

  private:
    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Paralution GMRES
   * solver.
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin, 2014
   */
  class SolverFGMRES : public SolverBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the verbosity of the solver to zero and
       * set the size of the Krylov space to 30, i.e. do a restart every 30
       * iterations.
       */
      AdditionalData (const unsigned int verbose = 0,
                      const unsigned int restart_parameter = 30);

      /**
       * Verbosity of the solver: 0 =  no output, 1 = print information at the
       * start and at the end, 3 = print residual at each iteration.
       */
      unsigned int verbose;

      /**
       * Size of the Krylov space.
       */
      unsigned int restart_parameter;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particular solver.
     */
    SolverFGMRES (SolverControl        &cn,
                  const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the GMRES solver of
     * Paralution. If the flag @p move_to_accelerator is set to true, the
     * solver is built on the accelerator.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number>     &A,
                Vector<Number>                 &x,
                const Vector<Number>           &b,
                const PreconditionBase<Number> &preconditioner,
                bool                            move_to_accelerator=false);

  private:
    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };



  /**
   * An implementation of the solver interface using the Paralution GMRES
   * solver.
   *
   * @ingroup ParalutionWrappers
   * @author Bruno Turcksin, 2014
   */
  class SolverAMG
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. @p relaxation_parameter is used only with smoothed
       * aggregation and is disregarded otherwise. @p over_interpolation is used
       * only with aggregation and is disregared otherwise.
       */
      AdditionalData (const unsigned int verbose = 0,
                      const mg_solver coarse_solver = gmres,
                      const unsigned int n_unknowns_coarse_level = 300,
                      const mg_solver smoother = richardson,
                      const mg_preconditioner = multicolored_sor,
                      const mg_cycle cycle = V_cycle,
                      const mg_interpolation interpolation = smoothed_aggregation,
                      const double coupling_strength = 0.01,
                      const double relaxation_parameter = 2./3.,
                      const double over_interpolation = 1.5);

      /**
       * Verbosity: 0 =  no output, 1 = print information at the
       * start and at the end, 3 = print residual at each iteration.
       */
      unsigned int verbose;

      /**
       * Solver of the coarse system.
       */
      mg_solver coarse_solver;

      /**
       * Number of unknowns on the coarsest level.
       */
      unsigned int n_unknowns_coarse_level;

      /**
       * Smoother (pre- and post-smoothers).
       */
      mg_solver smoother;

      /**
       * Preconditioner of the smoother.
       */
      mg_preconditioner preconditioner;

      /**
       * Cycle used by the multigrid preconditioner: V cycle, W cycle, K
       * cycle, or F cycle.
       */
      mg_cycle cycle;

      /**
       * Interpolation used by the multigrid preconditioner: smoothed
       * aggregation or aggregation.
       */
      mg_interpolation interpolation;

      /**
       * Coupling strength.
       */
      double coupling_strength;

      /**
       * Relaxation parameter for smoothed interpolation aggregation.
       */
      double relaxation_parameter;

      /**
       * Over-interpolation parameter for aggregation.
       */
      double over_interpolation;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particular solver.
     */
    SolverAMG (SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the AMG solver of
     * Paralution. If the flag @p move_to_accelerator is set to true, the
     * solver, the right-hand side, the matrix, and the solution vector are
     * moved on the accelerator once the solver is built. The multigrid must be
     * build on the CPU.
     */
    template <typename Number>
    void solve (SparseMatrix<Number>     &A,
                Vector<Number>           &x,
                Vector<Number>           &b,
                bool                      move_to_accelerator=false);

    /**
     * Access to object that controls convergence.
     */
    SolverControl &control() const;

  private:
    /**
     * Reference to the object that controls convergence of the iterative
     * solver. In fact, for these Paralution wrappers, Paralution does so
     * itself, but we copy the data from this object before starting solution
     * process, and copy the data back into it afterwards.
     */
    SolverControl &solver_control;

    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION

/*----------------------------   paralution_solver.h     ---------------------------*/

#endif
/*----------------------------   paralution_solver.h     ---------------------------*/
