// ---------------------------------------------------------------------
// $Id: paralution_solver.h 30040 2013-07-18 17:06:48Z maier $
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

#ifndef __deal2__paralution_solver_h
#define __deal2__paralution_solver_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PARALUTION

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/paralution_vector.h>
#include <deal.II/lac/paralution_sparse_matrix.h>


DEAL_II_NAMESPACE_OPEN

namespace ParalutionWrappers
{
  /**
   * Base class for solver classes using the Paralution solvers.
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
     * Reference to the object that controls convergence of the iteratove
     * solver. In fact, for these Paralution wrappers, Paralution does so
     * itself, but we copy the data from this object solition process, and
     * copy the data back into it afterwards.
     */
    SolverControl &solver_control;
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
     * Constructor.
     */
    SolverCG (SolverControl &cn);

    /**
     * Solve the linear system <tt>Ax=b</tt> using the CG solver of Paralution.
     */
    //TODO: add a flag to move to the accelerator
    template <typename Number>
    void solve (const SparseMatrix<Number> &A,
                Vector<Number>             &x,
                const Vector<Number>       &b);
  };



  /**
   * An implementation of the solver interface using the Paralution BiCGStab
   * solver.
   */
  class SolverBicgstab : public SolverBase
  {
  public:
    /**
     * Constructor.
     */
    SolverBicgstab (SolverControl &cn);

    /**
     * Solve the linear system <tt>Ax=b</tt> using the BiCGStab solver of
     * Paralution.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number> &A,
                Vector<Number>             &x,
                const Vector<Number>       &b);
  };



  /**
   * An implementation of the solver interface using the Paralution GMRES
   * solver.
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
       * Constructor. By default, set the size of the Krylov space to 30,
       * i.e. do a restart every 30 iterations.
       */
      AdditionalData (const unsigned int restart_parameter = 30);

      /**
       * Size of the Krylov space.
       */
      unsigned int restart_parameter;
    };

    /**
     * Constructor. AdditionalData is a structure that contains additional
     * flags for tuning this particulat solver.
     */
    SolverGMRES (SolverControl &cn,
                 const AdditionalData &data = AdditionalData());

    /**
     * Solve the linear system <tt>Ax=b</tt> using the GMRES solver of
     * Paralution.
     */
    template <typename Number>
    void solve (const SparseMatrix<Number> &A,
                Vector<Number>             &x,
                const Vector<Number>       &b);

  private:
    /**
     * Store a copy of the flags for this solver.
     */
    const AdditionalData additional_data;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION

/*----------------------------   trilinos_solver.h     ---------------------------*/

#endif
/*----------------------------   trilinos_solver.h     ---------------------------*/
