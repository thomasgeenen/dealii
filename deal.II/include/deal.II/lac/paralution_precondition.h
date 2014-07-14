// ---------------------------------------------------------------------
// $Id: paralution_precondition.h 30040 2013-07-18 17:06:48Z maier $
//
// Copyright (C) 2014 by the deal.II authors
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

#ifndef __deal2__paralution_precondition_h
#define __deal2__paralution_precondition_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PARALUTION

#include <deal.II/base/std_cxx1x/shared_ptr.h>

#include <paralution.hpp>

DEAL_II_NAMESPACE_OPEN


namespace ParalutionWrappers
{
  /**
   * Cycle for multigrid.
   */
  enum mg_cycle {V_cycle, W_cycle, K_cycle, F_cycle};

  /**
   * Interpolation for grid transfer operator.
   */
  enum mg_interpolation {aggregation, smoothed_aggregation};

  /**
   * Preconditioner for the smoother of the multigrid.
   */
  enum mg_preconditioner {jacobi, sgs, multicolored_sgs, multicolored_sor, ilu,
                          ilut, multicolored_ilu
                         };

  /**
   * Solver used as smoother and to solve the coarse system of the multigrid.
   */
  enum mg_solver {richardson, cg, cr, bicgstab, gmres, fgmres};

  /**
   * The base class for the preconditioner classes using the Paralution
   * functionality. In Paralution, the preconditioners are not built at the same
   * time as the solver. Therefore, only the declaration is done in the
   * preconditiner classes.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionBase
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. Does not do anything.
     */
    PreconditionBase ();

    friend class SolverBase;

  protected :
    /**
     * This is a pointer to the preconditioner object that is used when
     * applying the preconditioner.
     */
    std_cxx1x::shared_ptr<paralution::Preconditioner<paralution::LocalMatrix<Number>,
              paralution::LocalVector<Number>,Number> > preconditioner;
  };



  /**
   * A wrapper for Jacobi preconditioner for Paralution matrices.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionJacobi : public PreconditionBase<Number>
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     */
    PreconditionJacobi ();

  private :
    /**
     * Store a copy of the flags for this
     * particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A wrapper for Symmetric Gauss-Seidel for Paralution matrices.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionSGS : public PreconditionBase<Number>
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     */
    PreconditionSGS ();

  private :
    /**
     * Store a copy of the flags for this
     * particular preconditioner.
     */
    AdditionalData additional_data;

  };



  /**
   * A wrapper for multi-colored Symmetric Gauss-Seidel for Paralution matrices.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionMultiColoredSGS : public PreconditionBase<Number>
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {};

    /**
     * Constructor.
     */
    PreconditionMultiColoredSGS ();

  private :
    /**
     * Store a copy of the flags for this
     * particular preconditioner.
     */
    AdditionalData additional_data;

  };



  /**
   * A wrapper for multi-colored SOR for Paralution matrices.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionMultiColoredSOR : public PreconditionBase<Number>
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the damping parameter to one.
       */
      AdditionalData (const double omega = 1);

      /**
       * Relaxation parameter.
       */
      double omega;
    };

    /**
     * Constructor. Take additional flags if there are any.
     */
    PreconditionMultiColoredSOR (const AdditionalData &additional_data = AdditionalData());

    /**
     * This function changes the value of the additional flags.
     */
    void initialize (const AdditionalData &additional_data = AdditionalData());

  private :
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A wrapper for ILU(p) for Paralution matrices.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionILU : public PreconditionBase<Number>
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the fill-in parameter to zero.
       */
      AdditionalData (const unsigned int levels = 0);

      /**
       * Fill-in parameter.
       */
      unsigned int levels;
    };

    /**
     * Constructor. Take additional flags if there are any.
     */
    PreconditionILU (const AdditionalData &additional_data = AdditionalData());

    /**
     * This function changes the value of the additional flags.
     */
    void initialize (const AdditionalData &additional_data = AdditionalData());

  private :
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A wrapper for ILUT(t,m) factorization based on threshold (t) and maximum
   * number of elements per row (m) for Paralution matrices. The preconditioner
   * can be initialized using only the threshild value or both the threshold
   * value and the maximum number of elements per row.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionILUT : public PreconditionBase<Number>
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the threshold to zero.
       */
      AdditionalData (const Number threshold = 0.);

      /**
       * Constructor.
       */
      AdditionalData (const Number threshold, const unsigned int max_row);

      /**
       * Threshold.
       */
      Number threshold;

      /**
       * Maximum number of elements per row.
       */
      unsigned int max_row;
    };

    /**
     * Constructor. Take additional flags if there are any.
     */
    PreconditionILUT (const AdditionalData &additional_data = AdditionalData());

    /**
     * This function changes the value of the additional flags.
     */
    void initialize (const AdditionalData &additional_data = AdditionalData());

  private :
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A wrapper for ILU(p,q) for Paralution matrices. ILU(p,q) is based on ILU(p)
   * with a power(q)-pattern method. Compared to the standard ILU(p), ILU(p,q)
   * provides a higher degree of parallelism of forward and backward
   * substitution.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionMultiColoredILU : public PreconditionBase<Number>
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the fill-in parameter to zero and the
       * power-pattern to one.
       */
      AdditionalData (const unsigned int levels = 0, const unsigned int power = 1);

      /**
       * Fill-in parameter.
       */
      unsigned int levels;

      /**
       * Power parameter.
       */
      unsigned int power;
    };

    /**
     * Constructor. Take additional flags if there are any.
     */
    PreconditionMultiColoredILU (const AdditionalData &additional_data = AdditionalData());

    /**
     * This function changes the value of additional flags.
     */
    void initialize (const AdditionalData &additional_data = AdditionalData());
  };

  /**
   * A simplified wrapper for AMG for Paralution matrices. To use other parameter
   * combinations use the Paralution directly. AMG can only be built on the host.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2014
   */
  template <typename Number>
  class PreconditionAMG : public PreconditionBase<Number>
  {
  public :
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
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
                      const unsigned int max_iter = 2,
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
       * Maximum number of iterations.
       */
      unsigned int max_iter;

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
     * Constructor. Take additional flags if there are any.
     */
    PreconditionAMG (const AdditionalData &additional_data = AdditionalData());

    /**
     * Destructor.
     */
    ~PreconditionAMG();

    /**
     * This function changes the value of the additional flags.
     */
    void initialize (const AdditionalData &additional_data = AdditionalData());

  private :
    /**
     * Free the coarse level solver.
     */
    void clear();

    /**
     * Coarse grid solver.
     */
    AdditionalData additional_data;
  };
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION

/*----------------------------   paralution_precondition.h     ---------------------------*/

#endif
/*----------------------------   paralution_precondition.h     ---------------------------*/
