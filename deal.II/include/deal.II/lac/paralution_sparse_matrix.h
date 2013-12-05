// ---------------------------------------------------------------------
// $Id: paralution_sparse_matrix.h 30040 2013-07-18 17:06:48Z maier $
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

#ifndef __deal2__paralution_matrix_h
#define __deal2__paralution_matrix_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PARALUTION

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/sparse_matrix.h>
#include <algorithm>

#include "paralution.hpp"
#include "base/local_matrix.hpp"

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup ParalutionWrappers
 *@{
 */
namespace ParalutionWrappers
{
  /**
   * This class implements a wrapper to use the Paralution sparse matrix class
   * LocalMatrix. This class is designed for use in either serial
   * implementation or as a localized copy on each processors.
   *
   * The interface of this class is modeled after the existing SparseMatrix
   * class in deal.II. It has almost the same member functions, and is often
   * exchangeable. However, Paralution LocalMatix only supports float and double.
   *
   * @ingroup ParalutionWrappers
   * @ingroup Matrix1
   * @author Bruno Turcksin, 2013
   */
  template <typename Number>
  class SparseMatrix : public Subscriptor
  {
  public :
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    /**
     * @name 1: Constructors and initialization
     */
    //@{
    /**
     * Constructor; initializes the matrix to be empty, without any
     * structure, i.e. the matrix is not usable at all. This constructor is
     * therefore only useful for matrices which are members of a class. All
     * other matrices should be created at a point in the data flow where
     * all necessary information is available.
     *
     * You have to initialize the matrix before usage with reinit(const
     * SparsityPattern &sparsity_pattern).
     */
    SparseMatrix();

    /**
     * Copy constructor. This constructor is only allowed to be called if
     * the matrix to be copied is empty. This is for the same reason as for
     * the SparsityPattern, see there fore the details.
     *
     * If you really want to copy a whole matrix, you can do so by using the
     * copy_from() function.
     */
    //TODO
    //SparseMarix(const SparseMatrix &sparse_matrix);

    /**
     * Constructo. Takes the given matrix sparsity structure to represent
     * the sparsity pattern of this matrix. You can change the sparsity
     * pattern later on by calling the reinit(const SparsityPattern&)
     * function.
     *
     * The constructor is marked explicit so as to disallow that someone
     * passes a sparsity pattern in place of a spaese matrix to some
     * function, where an empty matrix would be generated then.
     */
    explicit SparseMatrix(const SparsityPattern &sparsity_pattern);

    /**
     * Destructor. Free all memory, but do not release the memory of the
     * sparsity structure.
     */
    ~SparseMatrix();

    /**
     * Reinitialize the sparse matrix with the given sparsity pattern. The
     * latter tells the matrix how many nonzero elements there need to be
     * reserved.
     *
     * The elements of the matrix are set to zero by this function.
     */
    void reinit(SparsityPattern const &sparsity_pattern);

    /**
     * Release all memory and return to a state just like after having
     * called the default constructor. It also forgets the sparsity pattern
     * it was previously tied to.
     */
    virtual void clear();

    /**
     * This function convert the underlying SparseMatrix to
     * Paralution::LocalMatrix. This function frees the SparseMatrix.
     */
    void convert_to_paralution_csr();
    //@}
    /**
     * @name 2: Information on the matrix
     */
    //@{
    /**
     * Return the dimension of the image space. To remember: the matrix is
     * of dimension $m \times n$. This function works only after
     * convert_to_paralution_csr has been called.
     */
    //TODO make the function work all the time
    size_type m() const;

    /**
     * Return the dimension of the range space. To remember: the matrix is
     * of dimension $m \times n$.
     */
    //TODO make the function work all the time
    size_type n() const;
    //@}
    /**
     * @name 3: Modifying entries
     */
    //@{
    /**
     * Add <tt>value</tt> to the element (<i>i,j</i>). Throws an error if
     * the entry does not exist or if <tt>value</tt> is not a finite number.
     * Still it is allowed to store zerp values in non-existent fields.
     *
     * This function only works on the host.
     */
    void add(const size_type i,
             const size_type j,
             const Number value);

    /**
     * Add all elements given in a FullMatrix into sparse matrix locations
     * given by <tt>indices</tt>. In other words, this function adds the
     * elements in <tt>full_matrix</tt> to the respective entries in calling
     * matrix, using the local-to-global indexing specified by
     * <tt>indices</tt> for both the rows and the columns of the matrix.
     * This function assumes a quadratic sparse matrix and a quadratic
     * full_matrix, the usual situation in FE calculations.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The defaul value is
     * <tt>true</tt>, i.e., zero values won't be added into the matrix.
     *
     * This function only works on the host.
     */
    template <typename Number2>
    void add(const std::vector<size_type> &indices,
             const FullMatrix<Number2> &full_matrix,
             const bool elide_zero_values = true);

    /**
     * Same function as before, but now including the possibility to use
     * rectaangular full_matrices and different local-to-global indexing on
     * rows and columns, respectively.
     *
     * This function only works on the host.
     */
    template <typename Number2>
    void add (const std::vector<size_type> &row_indices,
              const std::vector<size_type> &col_indices,
              const FullMatrix<Number2>    &full_matrix,
              const bool                    elide_zero_values = true);

    /**
     * Set several elements in the specified row of the matrix with column
     * indices as given by <tt>col_indices</tt> to the respective value.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The default value is
     * <tt>true</tt>, i.e., zero values won't be added into the matrix.
     *
     * This function only works on the host.
     */
    template <typename Number2>
    void add (const size_type               row,
              const std::vector<size_type> &col_indices,
              const std::vector<Number2>   &values,
              const bool                    elide_zero_values = true);

    /**
     * Add an array of values given by <tt>values</tt> in the given global
     * matrix row at columns specified by col_indices in the sparse_matrix.
     *
     * The optional parameter <tt>elide_zero_values</tt> can be used to
     * specify whether zero values should be added anyway or these should be
     * filtered away and only non-zero data is added. The default value is
     * <tt>true></tt>, i.e., zero values won't be added into the matrix.
     *
     * This function only works on the host.
     */
    template <typename Number2>
    void add (const size_type  row,
              const size_type  n_cols,
              const size_type *col_indices,
              const Number2   *values,
              const bool       elide_zero_values = true,
              const bool       col_indices_are_sorted = false);

    //@}
    /**
     * @name 4: Access to underlying Paralution data
     */
    /**
     * Return a constant reference to the underlying Paralution LocalMatrix
     * data.
     */
    paralution::LocalMatrix<Number> const &paralution_matrix() const;

    /**
     * Return a reference to the underlying Paralution LocalMatrix data.
     */
    paralution::LocalMatrix<Number>& paralution_matrix();
    //@}

  private :
    /**
     * Underlying Paralution LocalMatrix<Number>.
     */
    paralution::LocalMatrix<Number> local_matrix;

    /**
     * Temporary SparseMatrix used to build the Paralution LocalMatrix.
     */
    ::dealii::SparseMatrix<Number> sparse_matrix;
  };



// ------------------- inline functions --------------

  template <typename Number>
  inline SparseMatrix<Number>::SparseMatrix() {}



  template <typename Number>
  inline SparseMatrix<Number>::SparseMatrix(SparsityPattern const &sparsity_pattern)
  {
    reinit(sparsity_pattern);
  }



  template <typename Number>
  inline SparseMatrix<Number>::~SparseMatrix()
  {
    local_matrix.clear();
  }



  template <typename Number>
  inline void SparseMatrix<Number>::reinit(SparsityPattern const &sparsity_pattern)
  {
    local_matrix.clear();
    sparse_matrix.reinit(sparsity_pattern);
  }



  template <typename Number>
  inline typename SparseMatrix<Number>::size_type SparseMatrix<Number>::m() const
  {
    return local_matrix.get_nrows();
  }



  template <typename Number>
  inline typename SparseMatrix<Number>::size_type SparseMatrix<Number>::n() const
  {
    return local_matrix.get_ncols();
  }



  template <typename Number>
  inline void SparseMatrix<Number>::add(const size_type i,
                                        const size_type j,
                                        const Number    value)
  {
    sparse_matrix.add(i,j,value);
  }



  template <typename Number>
  template <typename Number2>
  inline void SparseMatrix<Number>::add(const std::vector<size_type> &indices,
                                        const FullMatrix<Number2>    &full_matrix,
                                        const bool                    elide_zero_values)
  {
    sparse_matrix.add(indices,full_matrix,elide_zero_values);
  }



  template <typename Number>
  template <typename Number2>
  inline void SparseMatrix<Number>::add(const size_type               row,
                                        const std::vector<size_type> &col_indices,
                                        const std::vector<Number2>   &values,
                                        const bool                    elide_zero_values)
  {
    sparse_matrix.add(row,col_indices,values,elide_zero_values);
  }

  template <typename Number>
  template <typename Number2>
  inline void SparseMatrix<Number>::add(const size_type  row,
                                        const size_type  n_cols,
                                        const size_type *col_indices,
                                        const Number2   *values,
                                        const bool       elide_zero_values,
                                        const bool       col_indices_are_sorted)
  {
    sparse_matrix.add(row,n_cols,col_indices,values,elide_zero_values,col_indices_are_sorted);
  }



  template <typename Number>
  inline paralution::LocalMatrix<Number> const &SparseMatrix<Number>::paralution_matrix() const
  {
    return local_matrix;
  }



  template <typename Number>
  inline paralution::LocalMatrix<Number>& SparseMatrix<Number>::paralution_matrix()
  {
    return local_matrix;
  }
}


/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION

/*----------------------------   paralution_sparse_matrix.h     ---------------------------*/

#endif
/*----------------------------   paralution_sparse_matrix.h     ---------------------------*/
