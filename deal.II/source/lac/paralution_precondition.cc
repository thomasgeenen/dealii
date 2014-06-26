// ---------------------------------------------------------------------
// $Id: paralution_precondition.cc 31349 2013-10-20 19:07:06Z maier $
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

#include <deal.II/lac/paralution_precondition.h>

#ifdef DEAL_II_WITH_PARALUTION

DEAL_II_NAMESPACE_OPEN

namespace ParalutionWrappers
{
  template <typename Number>
  PreconditionBase<Number>::PreconditionBase()
  {}



  /* -------------------------- PreconditionJacobi -------------------------- */

  template <typename Number>
  PreconditionJacobi<Number>::PreconditionJacobi()
  {
    this->preconditioner.reset(new paralution::Jacobi<paralution::LocalMatrix<Number>,
                               paralution::LocalVector<Number>,Number>);
  }



  /* -------------------------- PreconditionSGS ----------------------------- */

  template <typename Number>
  PreconditionSGS<Number>::PreconditionSGS()
  {
    this->preconditioner.reset(new paralution::SGS<paralution::LocalMatrix<Number>,
                               paralution::LocalVector<Number>,Number>);
  }



  /* -------------------------- PreconditionMultiColoredSGS ----------------- */

  template <typename Number>
  PreconditionMultiColoredSGS<Number>::PreconditionMultiColoredSGS()
  {
    this->preconditioner.reset(new paralution::MultiColoredSGS<paralution::LocalMatrix<Number>,
                               paralution::LocalVector<Number>,Number>);
  }



  /* -------------------------- PreconditionMulticoloredSOR ----------------- */

  template <typename Number>
  PreconditionMultiColoredSOR<Number>::AdditionalData::
  AdditionalData(const double omega)
    :
    omega(omega)
  {}



  template <typename Number>
  PreconditionMultiColoredSOR<Number>::PreconditionMultiColoredSOR(
    const AdditionalData &additional_data)
    :
    additional_data(additional_data)
  {
    this->preconditioner.reset(new paralution::MultiColoredGS<paralution::LocalMatrix<Number>,
                               paralution::LocalVector<Number>,Number>);

    initialize(additional_data);
  }



  template <typename Number>
  void PreconditionMultiColoredSOR<Number>::initialize(const AdditionalData &additional_data)
  {
    // Downcast the preconditioner pointer to use SetRelaxation
    paralution::MultiColoredGS<paralution::LocalMatrix<Number>,paralution::
    LocalVector<Number>,Number>* downcasted_ptr = static_cast<paralution::
                                                  MultiColoredGS<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
                                                  Number>* >(this->preconditioner.get());

    downcasted_ptr->SetRelaxation(additional_data.omega);
  }



  /* -------------------------- PreconditionILU ----------------------------- */

  template <typename Number>
  PreconditionILU<Number>::AdditionalData::AdditionalData(const unsigned int levels)
    :
    levels(levels)
  {}



  template <typename Number>
  PreconditionILU<Number>::PreconditionILU(const AdditionalData &additional_data)
    :
    additional_data(additional_data)
  {
    this->preconditioner.reset(new paralution::ILU<paralution::LocalMatrix<Number>,
                               paralution::LocalVector<Number>,Number>);

    initialize(additional_data);
  }



  template <typename Number>
  void PreconditionILU<Number>::initialize(const AdditionalData &additional_data)
  {
    // Downcast the preconditioner pointer to use Set
    paralution::ILU<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
               Number>* downcasted_ptr = static_cast<paralution::ILU<paralution::
                                         LocalMatrix<Number>,paralution::LocalVector<Number>,Number>* >(this->preconditioner.get());

    downcasted_ptr->Set(additional_data.levels);
  }



  /* -------------------------- PreconditionILUT ---------------------------- */

  template <typename Number>
  PreconditionILUT<Number>::AdditionalData::AdditionalData(const Number threshold)
    :
    threshold(threshold),
    max_row(0)
  {}



  template <typename Number>
  PreconditionILUT<Number>::AdditionalData::AdditionalData(const Number threshold,
                                                           const unsigned int max_row)
    :
    threshold(threshold),
    max_row(max_row)
  {}



  template <typename Number>
  PreconditionILUT<Number>::PreconditionILUT(const AdditionalData &additional_data)
    :
    additional_data(additional_data)
  {
    this->preconditioner.reset(new paralution::ILUT<paralution::LocalMatrix<Number>,
                               paralution::LocalVector<Number>,Number>);

    initialize(additional_data);
  }



  template <typename Number>
  void PreconditionILUT<Number>::initialize(const AdditionalData &additional_data)
  {
    // Downcast the preconditioner pointer to use Set
    paralution::ILUT<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
               Number>* downcasted_ptr = static_cast<paralution::ILUT<paralution::
                                         LocalMatrix<Number>,paralution::LocalVector<Number>,Number>* >(this->preconditioner.get());

    if (additional_data.max_row==0)
      downcasted_ptr->Set(additional_data.threshold);
    else
      downcasted_ptr->Set(additional_data.threshold,additional_data.max_row);
  }



  /* -------------------------- PreconditionMulticoloredILU --------------- */

  template <typename Number>
  PreconditionMultiColoredILU<Number>::AdditionalData::
  AdditionalData(const unsigned int levels, const unsigned int power)
    :
    levels(levels),
    power(power)
  {}



  template <typename Number>
  PreconditionMultiColoredILU<Number>::
  PreconditionMultiColoredILU(const AdditionalData &additional_data)
    :
    additional_data(additional_data)
  {
    this->preconditioner.reset(new paralution::MultiColoredILU<paralution::LocalMatrix<Number>,
                               paralution::LocalVector<Number>,Number>);

    initialize(additional_data);
  }



  template <typename Number>
  void PreconditionMultiColoredILU<Number>::initialize(const AdditionalData &additional_data)
  {
    // Downcast the preconditioner pointer to use Set
    paralution::MultiColoredILU<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
               Number>* downcasted_ptr = static_cast<paralution::MultiColoredILU<paralution::
                                         LocalMatrix<Number>,paralution::LocalVector<Number>,Number>* >(this->preconditioner.get());

    downcasted_ptr->Set(additional_data.levels,additional_data.power);
  }



  /* -------------------------- PreconditionAMG --------------------------- */

  template <typename Number>
  PreconditionAMG<Number>::AdditionalData::AdditionalData(const unsigned int verbose,
                                                          const mg_solver coarse_solver,
                                                          const unsigned int n_unknowns_coarse_level,
                                                          const unsigned int max_iter,
                                                          const mg_cycle cycle,
                                                          const mg_interpolation interpolation,
                                                          const double coupling_strength,
                                                          const double relaxation_parameter,
                                                          const double over_interpolation)
    :
    verbose(verbose),
    coarse_solver(coarse_solver),
    n_unknowns_coarse_level(n_unknowns_coarse_level),
    max_iter(max_iter),
    cycle(cycle),
    interpolation(interpolation),
    coupling_strength(coupling_strength),
    relaxation_parameter(relaxation_parameter),
    over_interpolation(over_interpolation)
  {}



  template <typename Number>
  PreconditionAMG<Number>::PreconditionAMG(const AdditionalData &additional_data)
    :
    coarse_solver(NULL)
  {
    this->preconditioner.reset(new paralution::AMG<paralution::LocalMatrix<Number>,
                               paralution::LocalVector<Number>,Number>);

    initialize(additional_data);
  }



  template <typename Number>
  PreconditionAMG<Number>::~PreconditionAMG()
  {
    clear();
  }



  template <typename Number>
  void PreconditionAMG<Number>::initialize(const AdditionalData &additional_data)
  {
    // Downcast the preconditioner pointer
    paralution::AMG<paralution::LocalMatrix<Number>,paralution::LocalVector<Number>,
      Number>* downcasted_ptr = static_cast<paralution::AMG<paralution::
        LocalMatrix<Number>,paralution::LocalVector<Number>,Number>* >(this->preconditioner.get());

    // Set the maximum number of iterations.
    downcasted_ptr->InitMaxIter(additional_data.max_iter);

    // Set the verbosity of the multigrid.
    downcasted_ptr->Verbose(additional_data.verbose);

    // Free the coarse level solver if necessary.
    clear();

    // Set the number of unknowns on the coarsest level.
    downcasted_ptr->SetCoarsestLevel(additional_data.n_unknowns_coarse_level);

    // Set the interpolation type.
    switch (additional_data.interpolation)
    {
      case aggregation :
        {
          downcasted_ptr->SetInterpolation(paralution::Aggregation);
          downcasted_ptr->SetOverInterp(additional_data.over_interpolation);
          break;
        }
      case smoothed_aggregation :
        {
          downcasted_ptr->SetInterpolation(paralution::SmoothedAggregation);
          downcasted_ptr->SetInterpRelax(additional_data.relaxation_parameter);
          break;
        }
      default :
        {
          AssertThrow(false,ExcMessage("Unknown interpolation type for PreconditionAMG."));
        }
    }

    // Set the coupling strength
    downcasted_ptr->SetCouplingStrength(additional_data.coupling_strength);

    // Set the type of cycle
    switch (additional_data.cycle)
    {
      case V_cycle :
        {
          downcasted_ptr->SetCycle(paralution::Vcycle);
          break;
        }
      case  W_cycle :
        {
          downcasted_ptr->SetCycle(paralution::Wcycle);
          break;
        }
      case K_cycle :
        {
          downcasted_ptr->SetCycle(paralution::Kcycle);
          break;
        }
      case F_cycle :
        {
          downcasted_ptr->SetCycle(paralution::Fcycle);
        }
      default :
        {
          AssertThrow(false,ExcMessage("Unknown cycle type for PreconditionAMG."));
        }
    }

    // Use manual coarse grid solver.
    downcasted_ptr->SetManualSolver(true);


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

    // Set the coarse solver and the verbosity of the coarse solver.
    coarse_solver->Verbose(additional_data.verbose);
    downcasted_ptr->SetSolver(*coarse_solver);
  }



  template <typename Number>
  void PreconditionAMG<Number>::clear()
  {
    // Free the coarse solver.
    if (coarse_solver!=NULL)
      {
        delete coarse_solver;
        coarse_solver = NULL;
      }
  }

}

// Explicit instantiations
namespace ParalutionWrappers
{
  template class PreconditionBase<float>;
  template class PreconditionBase<double>;
  template class PreconditionJacobi<float>;
  template class PreconditionJacobi<double>;
  template class PreconditionSGS<float>;
  template class PreconditionSGS<double>;
  template class PreconditionMultiColoredSGS<float>;
  template class PreconditionMultiColoredSGS<double>;
  template class PreconditionMultiColoredSOR<float>;
  template class PreconditionMultiColoredSOR<double>;
  template class PreconditionILU<float>;
  template class PreconditionILU<double>;
  template class PreconditionILUT<float>;
  template class PreconditionILUT<double>;
  template class PreconditionMultiColoredILU<float>;
  template class PreconditionMultiColoredILU<double>;
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION
