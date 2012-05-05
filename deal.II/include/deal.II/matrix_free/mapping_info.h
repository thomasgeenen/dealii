//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__matrix_free_mapping_info_h
#define __deal2__matrix_free_mapping_info_h


#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/matrix_free/helper_functions.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
namespace MatrixFreeFunctions
{
  /**
   * The class that stores all geometry-dependent data related with cell
   * interiors for use in the matrix-free class.
   *
   * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
   */
  template <int dim, typename Number>
  struct MappingInfo
  {
    typedef VectorizedArray<Number>               vector_t;
    typedef Point<dim,vector_t>                             point;
    typedef Tensor<1,dim,vector_t>                          tensor1;
    typedef Tensor<2,dim,vector_t>                          tensor2;
    typedef Tensor<3,dim,vector_t>                          tensor3;
    typedef Tensor<1,(dim>1?dim*(dim-1)/2:1),Tensor<1,dim,vector_t> > tensorUT;
    static const std::size_t n_vectors
      = VectorizedArray<Number>::n_array_elements;

                                /**
                                 * Determines how many bits of an unsigned int
                                 * are used to distinguish the cell types
                                 * (Cartesian, with constant Jacobian, or
                                 * general)
                                 */
    static const std::size_t  n_cell_type_bits = 2;

                                /**
                                 * Determines how many types of different
                                 * cells can be detected at most. Corresponds
                                 * to the number of bits we reserved for it.
                                 */
    static const unsigned int n_cell_types = 1U<<n_cell_type_bits;

                                /**
                                 * Empty constructor.
                                 */
    MappingInfo();

                                /**
                                 * Computes the information in the given
                                 * cells. The cells are specified by the level
                                 * and the index within the level (as given by
                                 * CellIterator::level() and
                                 * CellIterator::index(), in order to allow
                                 * for different kinds of iterators,
                                 * e.g. standard DoFHandler, multigrid, etc.)
                                 * on a fixed Triangulation. In addition, a
                                 * mapping and several quadrature formulas are
                                 * given.
                                 */
    void initialize (const dealii::Triangulation<dim>                &tria,
                     const std::vector<std::pair<unsigned int,unsigned int> > &cells,
                     const std::vector<unsigned int>         &active_fe_index,
                     const Mapping<dim>                      &mapping,
                     const std::vector<dealii::hp::QCollection<1> >  &quad,
                     const UpdateFlags                        update_flags);

                                /**
                                 * Helper function to determine which update
                                 * flags must be set in the internal functions
                                 * to initialize all data as requested by the
                                 * user.
                                 */
    UpdateFlags
    compute_update_flags (const UpdateFlags                        update_flags,
                          const std::vector<dealii::hp::QCollection<1> >  &quad) const;

                                /**
                                 * Returns the type of a given cell as
                                 * detected during initialization.
                                 */
    unsigned int get_cell_type (const unsigned int cell_chunk_no) const
    {
      AssertIndexRange (cell_chunk_no, cell_type.size());
      return cell_type[cell_chunk_no] % n_cell_types;
    };

                                /**
                                 * Returns the type of a given cell as
                                 * detected during initialization.
                                 */
    unsigned int get_cell_data_index (const unsigned int cell_chunk_no) const
    {
      AssertIndexRange (cell_chunk_no, cell_type.size());
      return cell_type[cell_chunk_no] >> n_cell_type_bits;
    };

                                /**
                                 * Clears all data fields in this class.
                                 */
    void clear ();

                                /**
                                 * Returns the memory consumption of this
                                 * class in bytes.
                                 */
    std::size_t memory_consumption() const;

                                /**
                                 * Prints a detailed summary of memory
                                 * consumption in the different structures of
                                 * this class to the given output stream.
                                 */
    template <typename STREAM>
    void print_memory_consumption(STREAM         &out,
                                  const SizeInfo &size_info) const;

                                /**
                                 * Stores whether a cell is Cartesian, has
                                 * constant transform data (Jacobians) or is
                                 * general. cell_type % 4 gives this
                                 * information (0: Cartesian, 1: constant
                                 * Jacobian throughout cell, 2: general cell),
                                 * and cell_type / 4 gives the index in the
                                 * data field of where to find the information
                                 * in the fields Jacobian and JxW values
                                 * (except for quadrature points, for which
                                 * the index runs as usual).
                                 */
    std::vector<unsigned int> cell_type;

                                /**
                                 * The first field stores the inverse Jacobian
                                 * for Cartesian cells: There, it is a
                                 * diagonal rank-2 tensor, so we actually just
                                 * store a rank-1 tensor. It is the same on
                                 * all cells, therefore we only store it once
                                 * per cell, and use similarities from one
                                 * cell to another, too (on structured meshes,
                                 * there are usually many cells with the same
                                 * Jacobian).
                                 *
                                 * The second field stores the Jacobian
                                 * determinant for Cartesian cells (without
                                 * the quadrature weight, which depends on the
                                 * quadrature point, whereas the determinant
                                 * is the same on each quadrature point).
                                 */
    AlignedVector<std::pair<tensor1,vector_t> > cartesian;

                                /**
                                 * The first field stores the Jacobian for
                                 * non-Cartesian cells where all the Jacobians
                                 * on the cell are the same (i.e., constant,
                                 * which comes from a linear transformation
                                 * from unit to real cell). Also use
                                 * similarities from one cell to another (on
                                 * structured meshes, there are usually many
                                 * cells with the same Jacobian).
                                 *
                                 * The second field stores the Jacobian
                                 * determinant for non-Cartesian cells with
                                 * constant Jacobian throughout the cell
                                 * (without the quadrature weight, which
                                 * depends on the quadrature point, whereas
                                 * the determinant is the same on each
                                 * quadrature point).
                                 */
    AlignedVector<std::pair<tensor2,vector_t> > linear;

                                /**
                                 * Definition of a structure that stores data
                                 * that depends on the quadrature formula (if
                                 * we have more than one quadrature formula on
                                 * a given problem, these fields will be
                                 * different)
                                 */
    struct MappingInfoDependent
    {
                                /**
                                 * This field stores the row starts for the
                                 * inverse Jacobian transformations,
                                 * quadrature weights and second derivatives.
                                 */
      std::vector<unsigned int> rowstart_jacobians;

                                /**
                                 * This field stores the inverse Jacobian
                                 * transformation from unit to real cell,
                                 * which is needed for most gradient
                                 * transformations (corresponds to
                                 * FEValues::inverse_jacobian) for general
                                 * cells.
                                 */
      AlignedVector<tensor2> jacobians;

                                /**
                                 * This field stores the Jacobian
                                 * determinant times the quadrature weights
                                 * (JxW in deal.II speak) for general cells.
                                 */
      AlignedVector<vector_t> JxW_values;

                                /**
                                 * Stores the diagonal part of the gradient of
                                 * the inverse Jacobian transformation. The
                                 * first index runs over the derivatives
                                 * $\partial^2/\partial x_i^2$, the second
                                 * over the space coordinate. Needed for
                                 * computing the Laplacian of FE functions on
                                 * the real cell. Uses a separate storage from
                                 * the off-diagonal part $\partial^2/\partial
                                 * x_i \partial x_j, i\neq j$ because that is
                                 * only needed for computing a full Hessian.
                                 */
      AlignedVector<tensor2>  jacobians_grad_diag;

                                /**
                                 * Stores the off-diagonal part of the
                                 * gradient of the inverse Jacobian
                                 * transformation. Because of symmetry, only
                                 * the upper diagonal part is needed. The
                                 * first index runs through the derivatives
                                 * row-wise, i.e., $\partial^2/\partial x_1
                                 * \partial x_2$ first, then
                                 * $\partial^2/\partial x_1 \partial x_3$, and
                                 * so on. The second index is the spatial
                                 * coordinate. Not filled currently.
                                 */
      AlignedVector<tensorUT> jacobians_grad_upper;

                                /**
                                 * Stores the row start for quadrature points
                                 * in real coordinates for both types of
                                 * cells. Note that Cartesian cells will have
                                 * shorter fields (length is @p n_q_points_1d)
                                 * than non-Cartesian cells (length is @p
                                 * n_q_points).
                                 */
      std::vector<unsigned int> rowstart_q_points;

                                /**
                                 * Stores the quadrature points in real
                                 * coordinates for Cartesian cells (does not
                                 * need to store the full data on all points)
                                 */
      AlignedVector<point>    quadrature_points;

                                /**
                                 * The dim-dimensional quadrature formula
                                 * underlying the problem (constructed from a
                                 * 1D tensor product quadrature formula).
                                 */
      dealii::hp::QCollection<dim>    quadrature_formula;

                                /**
                                 * The (dim-1)-dimensional quadrature formula
                                 * underlying the problem (constructed from a
                                 * 1D tensor product quadrature formula).
                                 */
      dealii::hp::QCollection<dim-1>  quadrature_formula_faces;

                                /**
                                 * The number of quadrature points for the
                                 * current quadrature formula.
                                 */
      std::vector<unsigned int> n_q_points;

                                /**
                                 * The number of quadrature points for the
                                 * current quadrature formula when applied to
                                 * a face. Only set if the quadrature formula
                                 * is derived from a tensor product, since it
                                 * is not defined from the full quadrature
                                 * formula otherwise.
                                 */
      std::vector<unsigned int> n_q_points_face;

                                /**
                                 * The quadrature weights (vectorized data
                                 * format) on the unit cell.
                                 */
      std::vector<AlignedVector<vector_t> > quadrature_weights;

                                /**
                                 * This variable stores the number of
                                 * quadrature points for all quadrature
                                 * indices in the underlying element for
                                 * easier access to data in the hp case.
                                 */
      std::vector<unsigned int> quad_index_conversion;

                                /**
                                 * Returns the quadrature index for a given
                                 * number of quadrature points. If not in hp
                                 * mode or if the index is not found, this
                                 * function always returns index 0. Hence,
                                 * this function does not check whether the
                                 * given degree is actually present.
                                 */
      unsigned int
      quad_index_from_n_q_points (const unsigned int n_q_points) const
      {
        for (unsigned int i=0; i<quad_index_conversion.size(); ++i)
          if (n_q_points == quad_index_conversion[i])
            return i;
        return 0;
      }

      
                                /**
                                 * Prints a detailed summary of memory
                                 * consumption in the different structures of
                                 * this class to the given output stream.
                                 */
      template <typename STREAM>
      void print_memory_consumption(STREAM         &out,
                                    const SizeInfo &size_info) const;

                                /**
                                 * Returns the memory consumption in bytes.
                                 */
      std::size_t memory_consumption () const;
    };

                                /**
                                 * Contains all the stuff that depends on the
                                 * quadrature formula

                                 */
    std::vector<MappingInfoDependent> mapping_data_gen;

                                /**
                                 * Stores whether JxW values have been
                                 * initialized
                                 */
    bool JxW_values_initialized;

                                /**
                                 * Stores whether we computed second
                                 * derivatives.
                                 */
    bool second_derivatives_initialized;

                                /**
                                 * Stores whether we computed quadrature
                                 * points.
                                 */
    bool quadrature_points_initialized;

                                /**
                                 * Internal temporary data used for the
                                 * initialization.
                                 */
    struct CellData
    {
      CellData (const double jac_size_in) :
        jac_size (jac_size_in) {}

      void resize (const unsigned int size)
      {
        if (general_jac.size() != size)
          {
            quadrature_points.resize(size);
            general_jac.resize(size);
            general_jac_grad.resize(size);
          }
      }

      AlignedVector<tensor1>  quadrature_points;
      AlignedVector<tensor2>  general_jac;
      AlignedVector<tensor3>  general_jac_grad;
      tensor2                 const_jac;
      const double            jac_size;
    };

                                /**
                                 * Helper function called internally during
                                 * the initialize function.
                                 */
   void evaluate_on_cell (const dealii::Triangulation<dim> &tria,
                          const std::pair<unsigned int,unsigned int> *cells,
                          const unsigned int  cell,
                          const unsigned int  my_q,
                          unsigned int (&cell_t_prev)[n_vectors],
                          unsigned int (&cell_t)[n_vectors],
                          FEValues<dim,dim> &fe_values,
                          CellData          &cell_data) const;
  };

} // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif