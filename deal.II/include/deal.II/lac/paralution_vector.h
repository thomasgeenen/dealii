// ---------------------------------------------------------------------
// $Id: paralution_vector.h 30040 2013-07-18 17:06:48Z maier $
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

#ifndef __deal2__paralution_vector_h
#define __deal2__paralution_vector_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PARALUTION

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include "paralution.hpp"

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup ParalutionWrappers
 *@{
 */
namespace ParalutionWrappers
{
  /**
   * This class is a wrapper around the Paralution localized vector. This kind
   * of vector is designed for use in either serial implementation or as a
   * localized copy on each processor. The implementation of this class is
   * based on the Paralution vector class LocalVector. The only requirement
   * for this class to work is that Paralution is installed with the same
   * compiler as is used for compilation of deal.II.
   *
   * The interface of this class is modeled after the existing Vector class in
   * deal.II. It has almost the same member functions, and is often
   * exchangeable. However, Paralution LocalVector only supports int, float,
   * and double.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Vectors
   * @author Bruno Turcksin, 2013
   */
  //TODO: lots of functions are missing
  template <typename Number>
  class Vector : public Subscriptor
  {
  public :
    /**
     * Declare some the standard types used in all containers. These types
     * parallel those in the <tt>C</tt> standard libraries
     * <tt>vector<...></tt> class.
     */
    typedef dealii::types::global_dof_index size_type;
    typedef Number                         *iterator;
    typedef const Number                   *const_iterator;

    /**
     * @name 1: Basic Object-handling
     */
    //@{

    /**
     * Default constructor that generates and empy (zero size) vector. The
     * function <tt>reinit()</tt> will have to give the vector the correct
     * size.
     */
    Vector();

    /**
     * Copy constructor. Sets the dimension to that of the given vector, and
     * copies all the elements. The copied vector stays on the host/device
     * where it is create.
     */
    //TODO: Look to use choose between CopyFrom and Clone. Difference between
    //copy and clone: copy the vector stays on his host/device, with clone
    //the vector goes on the same host/device.
    Vector(const Vector<Number> &v);

    /**
     * Constructor. Set dimensionto @p n and initialize all the elements
     * with zero.
     *
     * The constructor is made explicit to avoid accidents like this:
     * <tt>v=0;</tt>. Presumably, the user want to set every element of the
     * vector to zero, but instead, what happens is this call:
     * <tt>v=Vector@<Number@>(0);</tt>, i.e. the vector is replaced by one
     * length zero.
     */
    explicit Vector(const size_type n);

    /**
     * Initialize the vector with a given range of values pointed to by the
     * iterators. This function is there in anlogy to the @p std::vector
     * class.
     */
    //TODO
    // template <typename InputIterator>
    // Vector (const InputIterator first,
    //         const InputIterator last);

    /**
     * Destructor.
     */
    ~Vector();

    /**
     * Change the dimension of the vector to @p N. The vector is filled with
     * zeros.
     */
    //TODO look to add fast
    void reinit(const size_type N);

    /**
     * Return dimension of the vector.
     */
    std::size_t size() const;

    /**
     * Make the @p Vector class a bit loke the <tt>vector<></tt> class of
     * the C++ standard library by returning iterators to the strat and end
     * of the elements of this vector. The iterator is created on the host
     * or the device and cannot be moved.
     */
    iterator begin();

    /**
     * Return constant iterator to the start of the vectors. The iterator is
     * created on the host of the device and cannot be moved.
     */
    const_iterator begin() const;

    /**
     * Return an iterator pointing to the element past the end of the array.
     * The iterator is created on the host or the device and cannot be
     * moved.
     */
    iterator end();

    /**
     * Return a constant iterator pointing to the element past the end of
     * the array. The iterator is created on the host or the device and
     * cannot be moved.
     */
    const_iterator end() const;
    //@}

    /**
     * @name 2: Data-Acess
     */
    //@{
    /**
     * Access the value of the @p ith component. Works only on the host.
     */
    Number operator() (const size_type i) const;

    /**
     * Access the @p ith component as writeable reference. Works only on the
     * host.
     */
    Number &operator() (const size_type i);

    /**
     * Access the value of the @p ith component. Works only on the host.
     *
     * Exactly the same as operator().
     */
    Number operator[] (const size_type i) const;

    /**
     * Access the @p ith component as a writeable reference.
     *
     * Exactly thte asame as operator().
     */
    Number &operator[] (const size_type i);

    /**
     * Return a constant reference to the underlying Paralution LocalVector
     * class.
     */
    const paralution::LocalVector<Number>& paralution_vector() const;

    /**
     * Return a (modifyable) reference to the underlying Paralution
     * LocalVector class.
     */
    paralution::LocalVector<Number>& paralution_vector();
    //@}

    /**
     * @name 3: Modification of vectors
     */
    //@{
    /**
     * Add the given vector to the present one.
     */
    Vector<Number>& operator+= (const Vector<Number> &v);

    /**
     * Substract the given vector from the present one.
     */
    Vector<Number>& operator-= (const Vector<Number> &v);

    /**
     * Addition of @p s to all components. Note that @p s is a scalar and
     * not a vector.
     */
    void add(const Number s);

    /**
     * Scale each element of the vector by a constant value.
     */
    Vector<Number>& operator*= (const Number factor);

    /**
     * Scale each element of the vector by the inverse of the given value.
     */
    Vector<Number>& operator/= (const Number factor);
    //@}

  private :
    /**
     * Underlying Paralution LocalVector<Number>.
     */
    paralution::LocalVector<Number> local_vector;
  };




// ------------------- inline functions --------------

  template <typename Number>
  inline Vector<Number>::Vector() {}



  template <typename Number>
  inline Vector<Number>::Vector(const Vector <Number> &v)
  {
    local_vector.CopyFrom(v.paralution_vector());
  }



  template <typename Number>
  inline Vector<Number>::Vector(const size_type n)
  {
    local_vector.Allocate("deal_ii_vector",n);
  }



  template <typename Number>
  inline Vector<Number>::~Vector()
  {
    local_vector.Clear();
  }



  template <typename Number>
  void Vector<Number>::reinit(const size_type n)
  {
    local_vector.Clear();
    local_vector.Allocate("deal_ii_vector",n);
  }



  template <typename Number>
  inline std::size_t Vector<Number>::size() const
  {
    return static_cast<size_type>(local_vector.get_size());
  }



  template <typename Number>
  inline typename Vector<Number>::iterator Vector<Number>::begin()
  {
    return &(local_vector[0]);
  }



  template <typename Number>
  inline typename Vector<Number>::const_iterator Vector<Number>::begin() const
  {
    return &(local_vector[0]);
  }



  template <typename Number>
  inline typename Vector<Number>::iterator Vector<Number>::end()
  {
    return &(local_vector[0])+local_vector.get_size();
  }



  template <typename Number>
  inline typename Vector<Number>::const_iterator Vector<Number>::end() const
  {
    return &(local_vector[0])+local_vector.get_size();
  }



  template <typename Number>
  inline Number Vector<Number>::operator() (const size_type i) const
  {
    AssertIndexRange(i,static_cast<size_type>(local_vector.get_size()));

    return local_vector[i];
  }



  template <typename Number>
  inline Number &Vector<Number>::operator() (const size_type i)
  {
    AssertIndexRange(i,static_cast<size_type>(local_vector.get_size()));

    return local_vector[i];
  }



  template <typename Number>
  inline Number Vector<Number>::operator[] (const size_type i) const
  {
    AssertIndexRange(i,static_cast<size_type>(local_vector.get_size()));

    return local_vector[i];
  }



  template <typename Number>
  inline Number &Vector<Number>::operator[] (const size_type i)
  {
    AssertIndexRange(i,static_cast<size_type>(local_vector.get_size()));

    return local_vector[i];
  }



  template <typename Number>
  inline paralution::LocalVector<Number> const &Vector<Number>::paralution_vector() const
  {
    return local_vector;
  }



  template <typename Number>
  inline paralution::LocalVector<Number>& Vector<Number>::paralution_vector()
  {
    return local_vector;
  }



  template <typename Number>
  inline Vector<Number>& Vector<Number>::operator+= (Vector<Number> const &v)
  {
    Assert(size()==v.size(),ExcDimensionMismatch(size(),v.size()));

    local_vector.ScaleAdd(1.,v.paralution_vector());

    return *this;
  }



  template <typename Number>
  inline Vector<Number>& Vector<Number>::operator-= (Vector<Number> const &v)
  {
    Assert(size()==v.size(),ExcDimensionMismatch(size(),v.size()));

    local_vector.ScaleAddScale(1.,v.paralution_vector(),-1.);

    return *this;
  }



  template <typename Number>
  inline void Vector<Number>::add(const Number s)
  {
    Assert(numbers::is_finite(s),ExcNumberNotFinite());

    size_type size = local_vector.get_size();
    for (size_type i=0; i<size; ++i)
      local_vector[i] += s;
  }



  template <typename Number>
  inline Vector<Number>& Vector<Number>::operator*= (const Number factor)
  {
    Assert(numbers::is_finite(factor),ExcNumberNotFinite());

    local_vector.Scale(factor);

    return *this;
  }



  template <typename Number>
  inline Vector<Number>& Vector<Number>::operator/= (const Number factor)
  {
    Assert(numbers::is_finite(factor),ExcNumberNotFinite());

    const Number inv_factor(1./factor);

    Assert(numbers::is_finite(inv_factor),ExcNumberNotFinite());

    local_vector.Scale(inv_factor);

    return *this;
  }
}

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION

/*----------------------------   paralution_vector.h     ---------------------------*/

#endif
/*----------------------------   paralution_vector.h     ---------------------------*/
