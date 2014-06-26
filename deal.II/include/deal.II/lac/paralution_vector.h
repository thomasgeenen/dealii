// ---------------------------------------------------------------------
// $Id: paralution_vector.h 30040 2013-07-18 17:06:48Z maier $
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

#ifndef __deal2__paralution_vector_h
#define __deal2__paralution_vector_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PARALUTION

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/memory_consumption.h>

#include <paralution.hpp>

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
   * @author Bruno Turcksin, 2013, 2014
   */
  template <typename Number>
  class Vector : public Subscriptor
  {
  public :
    /**
     * Declare some of the standard types used in all containers. These types
     * parallel those in the <tt>C</tt> standard libraries
     * <tt>vector<...></tt> class.
     */
    typedef dealii::types::global_dof_index size_type;
    typedef Number                          value_type;
    typedef Number                         *iterator;
    typedef const Number                   *const_iterator;

    /**
     * A variable that indicates whether this vector supports distributed data
     * storage. If true, then this vector also needs an appropriate compress()
     * function that allows communicating recent set or add operations to
     * individual elements to be communicated to other processors.
     *
     * For the current class, the variable equals false, since it does not
     * support parallel data storage.
     */
    static const bool supports_distributed_data = false;

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
     * Copy constructor. Set the dimension to that of the given vector, and
     * copies all the elements. If @p copy_backend is set to false, the copied
     * vector stays on the host/device where it is created. Otherwise the
     * copied vector is moved to the host/device of the given vector.
     */
    Vector(const Vector<Number> &v, bool copy_backend = false);

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
     * Destructor.
     */
    ~Vector();

    /**
     * This function does nothing but is there for compatibility with the
     * @p PETScWrappers::Vector class.
     *
     * For the PETSC vector wrapper class, this function compresses the
     * underlying representation of the PETSc object, i.e. flushes the
     * buffers of the vector obkect if it has any. This function is necessary
     * after writing into a vector element-by-element and before anything else
     * can be done on it.
     *
     * However, for the implementation of this class, it is immaterial and
     * thus an empty function.
     */
    void compress (::dealii::VectorOperation::values operation
                   =::dealii::VectorOperation::unknown) const;

    /**
     * Change the dimension of the vector to @p N. The vector is filled with
     * zeros.
     */
    void reinit(const size_type N);

    /**
     * Change the dimension of the vector to that of the vector v. The vector
     * is filled with zeros.
     */
    void reinit(const Vector<Number> &v);

    /**
     * Free the vector.
     */
    void clear();

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

    /**
     * Set all components of the vector to the given number s.
     */
    Vector<Number>& operator= (const Number s);

    /**
     * Copy the given vector in the present one. Resize if necessary.
     */
    Vector<Number>& operator= (const Vector<Number> &v);

    /**
     * Mean value of the elements of this vector.
     */
    Number mean_value() const;

    /**
     * $l_2$-norm of the vector. The ssquare root of the sum of the squares of
     * the elements.
     */
    Number l2_norm() const;

    /**
     * Return true if the vector contains ghost elements. Since this not a
     * distributed vector the method always returns false.
     */
    bool has_ghost_elements() const;

    /**
     * Returns true if the given global index is in the local range of this
     * processor. Since this is not a distributed vector the method always
     * returns true.
     */
    bool in_local_range (const size_type global_index) const;

    /**
     * Return an index set that describes wich elements of this vector are
     * owned by the current processor. Note that this index set does not
     * include elements this vector may store locally as ghost elements but
     * that are in fact owned by another processor. As a consequence, the
     * index sets returned returned on different processors if this is a
     * distributed vector will form disjoint sets that add up to the the
     * complete index set. Obviously, if a vector is created on only one
     * processor, then the result would satisfy
     * @code
     *   vec.locally_owned_elements() == complete_index_set(vec.size())
     * @endcode
     *
     * Since the current data type does not support parallel data storage
     * across different processors, the returned index set is the complete
     * index set.
     */
    IndexSet locally_owned_elements() const;
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
     * Exactly the same as operator().
     */
    Number &operator[] (const size_type i);

    /**
     * A collective get operation: instead of getting individual elements of a
     * vector, this function allows to get a whole set of elements at once.
     * The indices of the elements to be read are stated in the first
     * argument, the corresponding values are returned in the second.
     */
    void extract_subvector_to (const std::vector<size_type> &indices,
                               std::vector<Number>          &values) const;

    /**
     * Just as the above, but with pointers. Useful in minimizing copying of
     * the data around.
     */
    template <typename ForwardIterator, typename OutputIterator>
    void extract_subvector_to (ForwardIterator       indices_begin,
                               const ForwardIterator indices_end,
                               OutputIterator        values_begin) const;

    /**
     * Return a constant reference to the underlying Paralution LocalVector
     * class.
     */
    const paralution::LocalVector<Number>& paralution_vector() const;

    /**
     * Return a (modifiable) reference to the underlying Paralution
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
     * A collective add operation: This function adds a whole set of values
     * stored in @p values to the vector components specified by @p indices.
     */
    template <typename Number2>
    void add (const std::vector<size_type> &indices,
              const std::vector<Number2>   &values);

    /**
     * This is a second collective add operation. As a difference, this
     * function takes a deal.II vector of values.
     */
    template <typename Number2>
    void add (const std::vector<size_type>    &indices,
              const ::dealii::Vector<Number2> &values);

    /**
     * Take an address where <tt>n_eleements</tt> are stored contiguously and
     * add them into the vector. Handles all case which are not covered by the
     * other two <tt>add()</tt> functions above.
     */
    template <typename Number2>
    void add (const size_type  n_elements,
              const size_type *indices,
              const Number2   *values);

    /**
     * Addition of @p s to all components. Note that @p s is a scalar and
     * not a vector.
     */
    void add (const Number s);

    /**
     * Simple vector addition, equal to the <tt>operator +=</tt>.
     */
    void add (const Vector<Number> &v);

    /**
     * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
     */
    void add (const Number a, const Vector<Number> &V);

    /**
     * Scaling and simple vector addition, i.e. <tt>*this += a*V + v*W</tt>.
     */
    void add (const Number a, const Vector<Number> &V,
              const Number b, const Vector<Number> &W);

    /**
     * Scaling and simple vector addition, i.e. <tt>*this = s*(*this)+V</tt>.
     */
    void sadd (const Number s, const Vector<Number> &V);

    /**
     * Scaling and simple addition, i.e. <tt>*this = s*(*this)+a*V</tt>.
     */
    void sadd (const Number s, const Number a, const Vector<Number> &V);

    /**
     * Scaling and multiple addition.
     * <tt>*this = s*(*this) + a*v + b*W</tt>.
     */
    void sadd (const Number s, const Number a, const Vector<Number> &V,
               const Number b, const Vector<Number> &W);

    /**
     * Scaling and multiple addition.
     * <tt>*this = s*(*this) + a*v + b*W + c*X</tt>.
     */
    void sadd (const Number s, const Number a, const Vector<Number> &V,
               const Number v, const Vector<Number> &W, const Number c,
               const Vector<Number> &X);

    /**
     * Scale each element of the vector by a constant value.
     */
    Vector<Number>& operator*= (const Number factor);

    /**
     * Scale each element of the vector by the inverse of the given value.
     */
    Vector<Number>& operator/= (const Number factor);

    /**
     * Assignement <tt>*this = a*u</tt>.
     */
    void equ (const Number a, const Vector<Number> &u);

    /**
     * Assignement <tt>*this = a*u + b*v</tt>.
     */
    void equ (const Number a, const Vector<Number> &u,
              const Number b, const Vector<Number> &v);

    /**
     * Assignement <tt>*this = a*u + b*v + c*w</tt>.
     */
    void equ (const Number a, const Vector<Number> &u,
              const Number b, const Vector<Number> &v,
              const Number c, const Vector<Number> &w);
    /**
     * Return the scalar product <tt>*this^T *x</tt>.
     */
    Number scalar_product(const Vector<Number> &x) const;

    //@}

    /**
     * @name 4: Mixed stuff
     */
    //@{

    /**
     * Move the Vector to the accelerator.
     */
    void move_to_accelerator();

    /**
     * Move the Vector to the host.
     */
    void move_to_host();

    /**
     * Move the Vector to the accelerator. The function returns immediately and
     * performs the asynchronous transfer in the background.
     */
    void move_to_accelerator_async();

    /**
     * Move the Vector to the host. The function returns immediately and
     * performs the asynchronous transfer in the background.
     */
    void move_to_host_async();

    /**
     * Print to a stream. @p precision denotes the desired precision with
     * which values shall be printed, @p scientific whether scientific
     * notation shall be used. If @p across is @p true then the vector is
     * printed in a line, while if @p false then the elements are printed on a
     * separate line each.
     */
    void print (std::ostream      &out,         
                const unsigned int precision,   
                const bool         scientific,  
                const bool         across) const;

    /**
     * Synchronize the code when move_to_host_async or move_to_accelerator_async
     * is used.
     */
    void sync();

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t memory_consumption () const;

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
  inline void Vector<Number>::compress(::dealii::VectorOperation::values operation) const
  {}



  template <typename Number>
  void Vector<Number>::reinit(const size_type n)
  {
    local_vector.Clear();
    local_vector.Allocate("deal_ii_vector",n);
  }



  template <typename Number>
  inline void Vector<Number>::reinit(const Vector<Number> &v)
  {
    local_vector.Clear();
    local_vector.Allocate("deal_ii_vector",v.size());
  }



  template <typename Number>
  inline void Vector<Number>::clear()
  {
    local_vector.Clear();
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
  inline Vector<Number>& Vector<Number>::operator= (const Number s)
  {
    Assert(numbers::is_finite(s),ExcNumberNotFinite());

    local_vector.SetValues(s);

    return *this;
  }



  template <typename Number>
  inline Vector<Number>& Vector<Number>::operator= (const Vector<Number> &v)
  {
    local_vector.CopyFrom(v.paralution_vector());

    return *this;
  }



  template <typename Number>
  inline Number Vector<Number>::l2_norm() const
  {
    return local_vector.Norm();
  }



  template <typename Number>
  inline bool Vector<Number>::has_ghost_elements() const
  {
    return false;
  }



  template <typename Number>
  inline bool Vector<Number>::in_local_range(const size_type global_index) const
  {
    return true;
  }



  template <typename Number>
  inline IndexSet Vector<Number>::locally_owned_elements() const
  {
    return complete_index_set(size());
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



  // This function is not inlined but it needs to be in the header file because
  // of ForwardIterator and OutputIterator.
  template <typename Number>
  template <typename ForwardIterator, typename OutputIterator>
  void Vector<Number>::extract_subvector_to (ForwardIterator       indices_begin,
                                             const ForwardIterator indices_end,
                                             OutputIterator        values_begin) const
  {
    while (indices_begin != indices_end)
      {
        *values_begin = operator()(*indices_begin);
        ++indices_begin;
        ++values_begin;
      }
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

    local_vector.AddScale(v.paralution_vector(),1.);

    return *this;
  }



  template <typename Number>
  inline Vector<Number>& Vector<Number>::operator-= (Vector<Number> const &v)
  {
    Assert(size()==v.size(),ExcDimensionMismatch(size(),v.size()));

    local_vector.AddScale(v.paralution_vector(),-1.);

    return *this;
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



  template <typename Number>
  inline Number Vector<Number>::scalar_product(const Vector<Number> &x) const
  {
    return local_vector.Dot(x.paralution_vector());
  }



  template <typename Number>
  inline void Vector<Number>::move_to_accelerator()
  {
    local_vector.MoveToAccelerator();
  }



  template <typename Number>
  inline void Vector<Number>::move_to_host()
  {
    local_vector.MoveToHost();
  }



  template <typename Number>
  inline void Vector<Number>::move_to_accelerator_async()
  {
    local_vector.MoveToAcceleratorAsync();
  }


  template <typename Number>
  inline void Vector<Number>::move_to_host_async()
  {
    local_vector.MoveToHostAsync();
  }



  template <typename Number>
  inline void Vector<Number>::sync()
  {
    local_vector.Sync();
  }



  template <typename Number>
  std::size_t Vector<Number>::memory_consumption () const
  {
    //TODO: This is a very poor approzimation.
    return sizeof(*this) + (size() * sizeof(Number));
  }
}

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION

/*----------------------------   paralution_vector.h     ---------------------------*/

#endif
/*----------------------------   paralution_vector.h     ---------------------------*/
