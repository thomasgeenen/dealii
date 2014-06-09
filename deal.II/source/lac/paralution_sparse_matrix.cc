// ---------------------------------------------------------------------
// $Id: paralution_sparse_matrix.cc 31567 2013-11-06 18:01:36Z turcksin $
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

#include <deal.II/lac/paralution_sparse_matrix.h>

#ifdef DEAL_II_WITH_PARALUTION

DEAL_II_NAMESPACE_OPEN

namespace ParalutionWrappers
{
  template <typename Number>
  SparseMatrix<Number>::SparseMatrix(const SparseMatrix &sparse_matrix, bool copy_backend)
  {
    is_local_matrix = sparse_matrix.is_paralution_matrix();
    // This is costly but copy_from is used to have a consistent behavior between
    // sparse_matrix and local_matrix.
    if (is_local_matrix==false)
      sparse_matrix.copy_from(sparse_matrix.sparse_matrix());
    else
    {
      if (copy_backend==false)
        local_matrix.CopyFrom(sparse_matrix.paralution_matrix());
      else
        local_matrix.CloneFrom(sparse_matrix.paralution_matrix());
    }
  }



  template <typename Number>
  void SparseMatrix<Number>::convert_to_paralution_csr()
  {
    local_matrix.AllocateCSR("deal_ii_local_matrix",sparse_matrix.n_nonzero_elements(),
                             sparse_matrix.m(),sparse_matrix.n());

    //TODO: replace deprecated functions
    std::vector<Number> val;
    val.reserve(sparse_matrix.n_nonzero_elements());
    typename ::dealii::SparseMatrix<Number>::iterator it = sparse_matrix.begin();
    typename ::dealii::SparseMatrix<Number>::iterator end_it = sparse_matrix.end();
    for (; it!=end_it; ++it)
      val.push_back(it->value());

    const size_type *colnums_ptr = sparse_matrix.get_sparsity_pattern().get_column_numbers();
    const std::size_t *rowstart_ptr = sparse_matrix.get_sparsity_pattern().get_rowstart_indices();
    std::vector<int> colnums(colnums_ptr,colnums_ptr+sparse_matrix.n_nonzero_elements());
    std::vector<int> rowstart(rowstart_ptr,rowstart_ptr+sparse_matrix.m()+1);


    // do the copying around of entries so that the diagonal entry is in the
    // right place. note that this is easy to detect: since all entries apart
    // from the diagonal entry are sorted, we know that the diagonal entry is
    // in the wrong place if and only if its column index is larger than the
    // column index of the second entry in a row
    //
    // ignore rows with only one or no entry
    for (size_type row=0; row<sparse_matrix.m(); ++row)
      {
        // we may have to move some elements that are left of the diagonal
        // but presently after the diagonal entry to the left, whereas the
        // diagonal entry has to move to the right. we could first figure out
        // where to move everything to, but for simplicity we just make a
        // series of swaps instead (this is kind of a single run of
        // bubble-sort, which gives us the desired result since the array is
        // already "almost" sorted)
        //
        // in the first loop, the condition in the while-header also checks
        // that the row has at least two entries and that the diagonal entry
        // is really in the wrong place
        int cursor = rowstart[row];
        while ((cursor < rowstart[row+1]-1) &&
               (colnums[cursor] > colnums[cursor+1]))
          {
            std::swap (colnums[cursor], colnums[cursor+1]);
            std::swap (val[cursor], val[cursor+1]);
            ++cursor;
          }
      }

    local_matrix.CopyFromCSR(&rowstart[0],&colnums[0],&val[0]);

    // Free the memory used by sparse_matrix.
    sparse_matrix.clear();
    is_local_matrix = true;
  }



  template <typename Number>
  void SparseMatrix<Number>::transpose()
  {
    if (is_local_matrix==false)
    {
      AssertThrow(false,
          ExcMessage("Transpose cannot be used if the underlying matrix is a dealii::SparseMatrix."));
    }
    else
    {
      local_matrix.Transpose();
    }
  }



  template <typename Number>
  void SparseMatrix<Number>::copy_from(const SparseMatrix &sparse_matrix)
  {
    is_local_matrix = sparse_matrix.is_paralution_matrix();
    if (is_local_matrix==false)
      sparse_matrix.copy_from(sparse_matrix.sparse_matrix());
    else
      local_matrix.CopyFrom(sparse_matrix.paralution_matrix());
  }



  template <typename Number>
  void SparseMatrix<Number>::copy_from_async(const SparseMatrix &sparse_matrix)
  {
    is_local_matrix = sparse_matrix.is_paralution_matrix();
    if (is_local_matrix==false)
      sparse_matrix.copy_from(sparse_matrix.sparse_matrix());
    else
      local_matrix.CopyFromAsync(sparse_matrix.paralution_matrix());
  }

  template <typename Number>
  void SparseMatrix<Number>::convert_format(matrix_format format)
  {
    if (is_local_matrix==false)
    {
      switch (format)
      {
        case DENSE :
          {
            local_matrix.ConvertToDense();
            break;
          }
        case CSR :
          {
            local_matrix.ConvertToCSR();
            break;
          }
        case MCSR :
          {
            local_matrix.ConvertToMCSR();
            break;
          }
        case BCSR :
          {
            local_matrix.ConvertToBCSR();
            break;
          }
        case COO :
          {
            local_matrix.ConvertToDIA();
            break;
          }
        case ELL :
          {
            local_matrix.ConvertToELL();
            break;
          }
        case HYB :
          {
            local_matrix.ConvertToHYB();
            break;
          }
        default :
          {
            AssertThrow(false,ExcMessage("Wrong format of Paralution matrix."));
          }
      }
    }
  }
}

// Explicit instantiations
namespace ParalutionWrappers
{
  template void SparseMatrix<float>::convert_to_paralution_csr();

  template void SparseMatrix<double>::convert_to_paralution_csr();
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION
