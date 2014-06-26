// ---------------------------------------------------------------------
// $Id: paralution_vector.cc 31349 2013-10-20 19:07:06Z maier $
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

#include <deal.II/lac/paralution_vector.h>

#ifdef DEAL_II_WITH_PARALUTION

DEAL_II_NAMESPACE_OPEN

namespace ParalutionWrappers
{

  template <typename Number>
  Vector<Number>::Vector(const Vector <Number> &v, bool copy_backend)
  {
    if (copy_backend)
      local_vector.CloneFrom(v.paralution_vector());
    else
      local_vector.CopyFrom(v.paralution_vector());
  }


  template <typename Number>
  Number Vector<Number>::mean_value() const
  {
    Number mean(0.);
    unsigned int i_max(size());
    for (unsigned int i=0; i<i_max; ++i)
      mean += local_vector[i];
    mean /= i_max;

    return mean;
  }



  template <typename Number>
  void Vector<Number>::extract_subvector_to (const std::vector<size_type> &indices,
                                             std::vector<Number>          &values) const
  {
    for (size_type i=0; i<indices.size(); ++i)
      values[i] = operator()(indices[i]);
  }



  template <typename Number>
  template <typename Number2>
  void Vector<Number>::add (const std::vector<size_type> &indices,
                            const std::vector<Number2>   &values)
  {
    Assert(indices.size()==values.size(),
           ExcDimensionMismatch(indices.size(),values.size()));
    add(indices.size(),&indices[0],&values[0]);
  }



  template <typename Number>
  template <typename Number2>
  void Vector<Number>::add (const std::vector<size_type>    &indices,
                            const ::dealii::Vector<Number2> &values)
  {
    Assert(indices.size()==values.size(),
           ExcDimensionMismatch(indices.size(),values.size()));
    add(indices.size(),&indices[0],values.begin());
  }



  template <typename Number>
  template <typename Number2>
  void Vector<Number>::add (const size_type n_indices,const size_type *indices,
                            const Number2 *values)
  {
    for (size_type i=0; i<n_indices; ++i)
      {
        Assert(indices[i]<size(),ExcIndexRange(indices[i],0,size()));
        Assert(numbers::is_finite(values[i]),ExcNumberNotFinite());

        local_vector[indices[i]] += values[i];
      }
  }



  template <typename Number>
  void Vector<Number>::add(const Number s)
  {
    Assert(numbers::is_finite(s),ExcNumberNotFinite());

    size_type size = local_vector.get_size();
    for (size_type i=0; i<size; ++i)
      local_vector[i] += s;
  }



  template <typename Number>
  void Vector<Number>::add(const Vector<Number> &V)
  {
    Assert(size()==V.size(),ExcDimensionMismatch(size(),V.size()));
    local_vector.AddScale(V.paralution_vector(),1.);
  }



  template <typename Number>
  void Vector<Number>::add(const Number a,const Vector<Number> &V)
  {
    Assert(numbers::is_finite(a),ExcNumberNotFinite());
    Assert(size()==V.size(),ExcDimensionMismatch(size(),V.size()));

    local_vector.ScaleAddScale(1.,V.paralution_vector(),a);
  }



  template <typename Number>
  void Vector<Number>::add(const Number a,const Vector<Number> &V,
                           const Number b,const Vector<Number> &W)
  {
    Assert(numbers::is_finite(a),ExcNumberNotFinite());
    Assert(numbers::is_finite(b),ExcNumberNotFinite());
    Assert(size()==V.size(),ExcDimensionMismatch(size(),V.size()));
    Assert(size()==W.size(),ExcDimensionMismatch(size(),W.size()));

    local_vector.ScaleAdd2(1.,V.paralution_vector(),a,W.paralution_vector(),b);
  }



  template <typename Number>
  void Vector<Number>::sadd(const Number s, const Vector<Number> &V)
  {
    Assert(numbers::is_finite(s),ExcNumberNotFinite());
    Assert(size()==V.size(),ExcDimensionMismatch(size(),V.size()));

    local_vector.ScaleAdd(s,V.paralution_vector());
  }



  template <typename Number>
  void Vector<Number>::sadd(const Number s, const Number a, const Vector<Number> &V)
  {
    Assert(numbers::is_finite(s),ExcNumberNotFinite());
    Assert(numbers::is_finite(a),ExcNumberNotFinite());
    Assert(size()==V.size(),ExcDimensionMismatch(size(),V.size()));

    local_vector.ScaleAddScale(s,V.paralution_vector(),a);
  }



  template <typename Number>
  void Vector<Number>::sadd(const Number s, const Number a, const Vector<Number> &V,
                            const Number b, const Vector<Number> &W)
  {
    Assert(numbers::is_finite(s),ExcNumberNotFinite());
    Assert(numbers::is_finite(a),ExcNumberNotFinite());
    Assert(numbers::is_finite(b),ExcNumberNotFinite());
    Assert(size()==V.size(),ExcDimensionMismatch(size(),V.size()));
    Assert(size()==W.size(),ExcDimensionMismatch(size(),W.size()));

    local_vector.ScaleAdd2(s,V.paralution_vector(),a,W.paralution_vector(),b);
  }



  template <typename Number>
  void Vector<Number>::sadd(const Number s, const Number a, const Vector<Number> &V,
                            const Number b, const Vector<Number> &W, const Number c,
                            const Vector<Number> &X)
  {
    Assert(numbers::is_finite(s),ExcNumberNotFinite());
    Assert(numbers::is_finite(a),ExcNumberNotFinite());
    Assert(numbers::is_finite(b),ExcNumberNotFinite());
    Assert(numbers::is_finite(c),ExcNumberNotFinite());
    Assert(size()==V.size(),ExcDimensionMismatch(size(),V.size()));
    Assert(size()==W.size(),ExcDimensionMismatch(size(),W.size()));
    Assert(size()==X.size(),ExcDimensionMismatch(size(),X.size()));

    local_vector.ScaleAdd2(s,V.paralution_vector(),a,W.paralution_vector(),b);
    local_vector.AddScale(X.paralution_vector(),c);
  }



  template <typename Number>
  void Vector<Number>::equ (const Number a, const Vector<Number> &u)
  {
    Assert(numbers::is_finite(a), ExcNumberNotFinite());
    Assert(size()==u.size(),ExcDimensionMismatch(size(),u.size()));

    local_vector.ScaleAddScale(0.,u.paralution_vector(),a);
  }



  template <typename Number>
  void Vector<Number>::equ (const Number a, const Vector<Number> &u,
                            const Number b, const Vector<Number> &v)
  {
    Assert(numbers::is_finite(a), ExcNumberNotFinite());
    Assert(numbers::is_finite(b), ExcNumberNotFinite());
    Assert(size()==u.size(),ExcDimensionMismatch(size(),u.size()));
    Assert(size()==v.size(),ExcDimensionMismatch(size(),v.size()));

    local_vector.ScaleAdd2(0.,u.paralution_vector(),a,v.paralution_vector(),b);
  }



  template <typename Number>
  void Vector<Number>::equ (const Number a, const Vector<Number> &u,
                            const Number b, const Vector<Number> &v,
                            const Number c, const Vector<Number> &w)
  {
    Assert(numbers::is_finite(a), ExcNumberNotFinite());
    Assert(numbers::is_finite(b), ExcNumberNotFinite());
    Assert(numbers::is_finite(c), ExcNumberNotFinite());
    Assert(size()==u.size(),ExcDimensionMismatch(size(),u.size()));
    Assert(size()==v.size(),ExcDimensionMismatch(size(),v.size()));
    Assert(size()==w.size(),ExcDimensionMismatch(size(),w.size()));

    local_vector.ScaleAdd2(0.,u.paralution_vector(),a,v.paralution_vector(),b);
    local_vector.AddScale(w.paralution_vector(),c);
  }



  template <typename Number>
  void Vector<Number>::print (std::ostream      &out,
                              const unsigned int precision,
                              const bool         scientific,
                              const bool         across) const
  {
    AssertThrow (out, ExcIO());

    // get a representation of the vector and loop over all the elements 
    out.precision (precision);
    if (scientific)
      out.setf (std::ios::scientific, std::ios::floatfield);
    else
      out.setf (std::ios::fixed, std::ios::floatfield);

    if (across)
      for (size_type i=0; i<size(); ++i)
        out << static_cast<double>(local_vector[i]) << ' ';
    else
      for (size_type i=0; i<size(); ++i)
        out << static_cast<double>(local_vector[i]) << std::endl;
    out << std::endl;

    // restore the representation
    // of the vector
    AssertThrow (out, ExcIO());
  }
}

// Explicit instantiations
namespace ParalutionWrappers
{
  template class Vector<float>;
  template class Vector<double>;
  template void Vector<float>::add<float> (const std::vector<size_type> &indices,
                                           const std::vector<float>     &values);
  template void Vector<float>::add<float> (const std::vector<size_type>  &indices,
                                           const ::dealii::Vector<float> &values);
  template void Vector<float>::add<float> (const size_type  n_elements,
                                           const size_type *indices,
                                           const float     *values);
  template void Vector<float>::add<double> (const std::vector<size_type> &indices,
                                            const std::vector<double>    &values);
  template void Vector<float>::add<double> (const std::vector<size_type>   &indices,
                                            const ::dealii::Vector<double> &values);
  template void Vector<float>::add<double> (const size_type  n_elements,
                                            const size_type *indices,
                                            const double    *values);
  template void Vector<double>::add<float> (const std::vector<size_type> &indices,
                                            const std::vector<float>     &values);
  template void Vector<double>::add<float> (const std::vector<size_type>  &indices,
                                            const ::dealii::Vector<float> &values);
  template void Vector<double>::add<float> (const size_type  n_elements,
                                            const size_type *indices,
                                            const float     *values);
  template void Vector<double>::add<double> (const std::vector<size_type> &indices,
                                             const std::vector<double>    &values);
  template void Vector<double>::add<double> (const std::vector<size_type>  &indices,
                                             const ::dealii::Vector<double> &values);
  template void Vector<double>::add<double> (const size_type  n_elements,
                                             const size_type *indices,
                                             const double    *values);
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION
