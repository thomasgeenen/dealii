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
