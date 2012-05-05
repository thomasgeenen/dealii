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

#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/matrix_free/mapping_info.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
namespace MatrixFreeFunctions
{
  // ----------------- actual MappingInfo functions -------------------------

  template <int dim, typename Number>
  MappingInfo<dim,Number>::MappingInfo()
    :
    JxW_values_initialized (false),
    second_derivatives_initialized (false),
    quadrature_points_initialized (false)
  {}



  template <int dim, typename Number>
  void
  MappingInfo<dim,Number>::clear ()
  {
    JxW_values_initialized = false;
    quadrature_points_initialized = false;
    second_derivatives_initialized = false;
    mapping_data_gen.clear();
    cell_type.clear();
    cartesian.clear();
    linear.clear();
  }



  template <int dim, typename Number>
  UpdateFlags
  MappingInfo<dim,Number>::
  compute_update_flags (const UpdateFlags update_flags,
                        const std::vector<dealii::hp::QCollection<1> > &quad) const
  {
                                // this class is build around the evaluation
                                // this class is build around the evaluation
                                // of inverse gradients, so compute them in
                                // any case
    UpdateFlags new_flags = update_inverse_jacobians;

                                // if the user requested gradients, need
                                // inverse Jacobians
    if (update_flags & update_gradients || update_flags & update_inverse_jacobians)
      new_flags |= update_inverse_jacobians;

                                // for JxW, would only need JxW values.
    if (update_flags & update_JxW_values)
      new_flags |= update_JxW_values;

                                // for Hessian information, need inverse
                                // Jacobians and the derivative of Jacobians
                                // (these two together will give use the
                                // gradients of the inverse Jacobians, which
                                // is what we need)
    if (update_flags & update_hessians || update_flags & update_jacobian_grads)
      new_flags |= update_jacobian_grads;

    if (update_flags & update_quadrature_points)
      new_flags |= update_quadrature_points;

                                // there is one more thing: if we have a
                                // quadrature formula with only one quadrature
                                // point on the first component, but more
                                // points on later components, we need to have
                                // Jacobian gradients anyway in order to
                                // determine whether the Jacobian is constant
                                // throughout a cell
    bool formula_with_one_point = false;
    for (unsigned int i=0; i<quad[0].size(); ++i)
      if (quad[0][i].size() == 1)
        {
          formula_with_one_point = true;
          break;
        }
    if (formula_with_one_point == true)
      for (unsigned int comp=1; comp<quad.size(); ++comp)
        for (unsigned int i=0; i<quad[comp].size(); ++i)
          if (quad[comp][i].size() > 1)
            {
              new_flags |= update_jacobian_grads;
              goto end_set;
            }
  end_set:
    return new_flags;
  }



  namespace internal
  {
    template <int dim>
    double get_jacobian_size (const dealii::Triangulation<dim> &tria)
    {
      if (tria.n_cells() == 0)
        return 1;
      else return tria.begin()->diameter();
    }
  }



  template <int dim, typename Number>
  void
  MappingInfo<dim,Number>::initialize
  (const dealii::Triangulation<dim>                                 &tria,
   const std::vector<std::pair<unsigned int,unsigned int> > &cells,
   const std::vector<unsigned int>                          &active_fe_index,
   const Mapping<dim>                                       &mapping,
   const std::vector<dealii::hp::QCollection<1> >                   &quad,
   const UpdateFlags                                         update_flags_input)
  {
    clear();
    const unsigned int n_quads = quad.size();
    const unsigned int n_cells = cells.size();
    Assert (n_cells%n_vectors == 0, ExcInternalError());
    const unsigned int n_macro_cells = n_cells/n_vectors;
    mapping_data_gen.resize (n_quads);
    cell_type.resize (n_macro_cells);

                                // dummy FE that is used to set up an FEValues
                                // object. Do not need the actual finite
                                // element because we will only evaluate
                                // quantities for the mapping that are
                                // independent of the FE
    FE_Nothing<dim> dummy_fe;
    UpdateFlags update_flags = compute_update_flags (update_flags_input, quad);

    if (update_flags & update_JxW_values)
      JxW_values_initialized = true;
    if (update_flags & update_jacobian_grads)
      second_derivatives_initialized = true;
    if (update_flags & update_quadrature_points)
      quadrature_points_initialized = true;

                                // when we make comparisons about the size of
                                // Jacobians we need to know the approximate
                                // size of typical entries in Jacobians. We
                                // need to fix the Jacobian size once and for
                                // all. We choose the diameter of the first
                                // cell (on level zero, which is the best
                                // accuracy we can hope for, since diameters
                                // on finer levels are computed by differences
                                // of nearby cells). If the mesh extends over
                                // a certain domain, the precision of double
                                // values is essentially limited by this
                                // precision.
    const double jacobian_size = internal::get_jacobian_size(tria);

                                // objects that hold the data for up to
                                // n_vectors cells while we fill them up. Only
                                // after all n_vectors cells have been
                                // processed, we can insert the data into the
                                // data structures of this class
    CellData data (jacobian_size);

    for (unsigned int my_q=0; my_q<n_quads; ++my_q)
      {
        MappingInfoDependent & current_data = mapping_data_gen[my_q];
        const unsigned int n_hp_quads = quad[my_q].size();
        AssertIndexRange (0, n_hp_quads);
        current_data.n_q_points.reserve (n_hp_quads);
        current_data.n_q_points_face.reserve (n_hp_quads);
        current_data.quadrature_weights.resize (n_hp_quads);
        std::vector<unsigned int> n_q_points_1d (n_hp_quads),
          step_size_cartesian (n_hp_quads);
        if (n_hp_quads > 1)
          current_data.quad_index_conversion.resize(n_hp_quads);
        for (unsigned int q=0; q<n_hp_quads; ++q)
          {
            n_q_points_1d[q] = quad[my_q][q].size();
            const unsigned int n_q_points =
              Utilities::fixed_power<dim>(n_q_points_1d[q]);
            current_data.n_q_points.push_back (n_q_points);

            current_data.n_q_points_face.push_back 
              (Utilities::fixed_power<dim-1>(n_q_points_1d[q]));
            current_data.quadrature_formula.push_back
              (Quadrature<dim>(quad[my_q][q]));
            current_data.quadrature_formula_faces.push_back
              (Quadrature<dim-1>(quad[my_q][q]));

                                // set quadrature weights in vectorized form
            current_data.quadrature_weights[q].resize(n_q_points);
            for (unsigned int i=0; i<n_q_points; ++i)
              current_data.quadrature_weights[q][i] =
                current_data.quadrature_formula[q].get_weights()[i];

            if (n_hp_quads > 1)
              current_data.quad_index_conversion[q] = n_q_points;

                                // To walk on the diagonal for lexicographic
                                // ordering, we have to jump one index ahead
                                // in each direction. For direction 0, this is
                                // just the next point, for direction 1, it
                                // means adding n_q_points_1d, and so on.
            step_size_cartesian[q] = 0;
            unsigned int factor = 1;
            for (unsigned int d=0; d<dim; ++d)
              {
                step_size_cartesian[q] += factor;
                factor *= n_q_points_1d[q];
              }
          }

                                // if there are no cells, there is nothing to
                                // do
        if (cells.size() == 0)
          continue;

        tensor3  jac_grad, grad_jac_inv;
        tensor1  tmp;

                                // encodes the cell types of the current
                                // cell. Since several cells must be
                                // considered together, this variable holds
                                // the individual info of the last chunk of
                                // cells
        unsigned int cell_t [n_vectors], cell_t_prev [n_vectors];
        for (unsigned int j=0; j<n_vectors; ++j)
          cell_t_prev[j] = numbers::invalid_unsigned_int;

                                // fe_values object that is used to compute
                                // the mapping data. for the hp case there
                                // might be more than one finite
                                // element. since we manually select the
                                // active FE index and not via a
                                // hp::DoFHandler<dim>::active_cell_iterator,
                                // we need to manually select the correct
                                // finite element, so just hold a vector of
                                // FEValues
        std::vector<std_cxx1x::shared_ptr<FEValues<dim> > >
          fe_values (current_data.quadrature_formula.size());
        UpdateFlags update_flags_feval =
          (update_flags & update_inverse_jacobians ? update_jacobians : update_default) |
          (update_flags & update_jacobian_grads ? update_jacobian_grads : update_default) |
          (update_flags & update_quadrature_points ? update_quadrature_points : update_default);

                                // resize the fields that have fixed size or
                                // for which we know something from an earlier
                                // loop
        current_data.rowstart_q_points.resize (n_macro_cells+1);
        if (my_q > 0)
          {
            const unsigned int n_cells_var =
              mapping_data_gen[0].rowstart_jacobians.size()-1;
            current_data.rowstart_jacobians.reserve (n_cells_var+1);
            const unsigned int reserve_size = n_cells_var *
              current_data.n_q_points[0];
            if (mapping_data_gen[0].jacobians.size() > 0)
              current_data.jacobians.reserve (reserve_size);
            if (mapping_data_gen[0].JxW_values.size() > 0)
              current_data.jacobians.reserve (reserve_size);
            if (mapping_data_gen[0].jacobians_grad_diag.size() > 0)
              current_data.jacobians_grad_diag.reserve (reserve_size);
            if (mapping_data_gen[0].jacobians_grad_upper.size() > 0)
              current_data.jacobians_grad_upper.reserve (reserve_size);
          }

                                // a hash structure that is used to detect
                                // similarities between mapping data from one
                                // cell to the next.
        std::vector<std::pair<unsigned int, int> > hash_collection;
        internal::HashValue hash_value (jacobian_size);

                                // loop over all cells
        for (unsigned int cell=0; cell<n_macro_cells; ++cell)
          {
                                // GENERAL OUTLINE: First generate the data in
                                // format "number" for n_vectors cells, and
                                // then find the most general type of cell for
                                // appropriate vectorized formats. then fill
                                // this data in
            const unsigned int fe_index = active_fe_index.size() > 0 ?
              active_fe_index[cell] : 0;
            const unsigned int n_q_points = current_data.n_q_points[fe_index];
            if (fe_values[fe_index].get() == 0)
              fe_values[fe_index].reset 
                (new FEValues<dim> (mapping, dummy_fe,
                                    current_data.quadrature_formula[fe_index],
                                    update_flags_feval));
            FEValues<dim> &fe_val = *fe_values[fe_index];
            data.resize (n_q_points);

                                // if the fe index has changed from the
                                // previous cell, set the old cell type to
                                // invalid (otherwise, we might detect
                                // similarity due to some cells further ahead)
            if (cell > 0 && active_fe_index.size() > 0 &&
                active_fe_index[cell] != active_fe_index[cell-1])
              cell_t_prev[n_vectors-1] = numbers::invalid_unsigned_int;
            evaluate_on_cell (tria, &cells[cell*n_vectors],
                              cell, my_q, cell_t_prev, cell_t, fe_val, data);

                                // now reorder the data into vectorized
                                // types. if we are here for the first time,
                                // we need to find out whether the Jacobian
                                // allows for some simplification (Cartesian,
                                // linear) taking n_vectors cell together and
                                // we have to insert that data into the
                                // respective fields. Also, we have to
                                // compress different cell indicators into one
                                // structure.

            if (my_q == 0)
              {
                                // find the most general cell type (most
                                // general type is 2 (general cell))
                unsigned int most_general_type = 0;
                for (unsigned int j=0; j<n_vectors; ++j)
                  if (cell_t[j] > most_general_type)
                    most_general_type = cell_t[j];
                AssertIndexRange (most_general_type, 3);

                                // Cartesian cell with diagonal Jacobian: only
                                // insert the diagonal of the inverse and the
                                // Jacobian determinant
                unsigned int insert_position = numbers::invalid_unsigned_int;
                typedef std::vector<std::pair<unsigned int,int> >::iterator iter;
                if (most_general_type == 0)
                  {
                    std::pair<tensor1,vector_t> new_entry;
                    for (unsigned int d=0; d<dim; ++d)
                      new_entry.first[d] = data.const_jac[d][d];
                    insert_position = cartesian.size();

                                // check whether everything is the same as on
                                // another cell before. find an insertion point
                                // in the list of keys that we have
                                // collected. put negative position so that
                                // the insertion position from lower_bound is
                                // that with the same key
                    const unsigned int hash = hash_value.template operator()<dim,Number> (data.const_jac, true);
                    std::pair<unsigned int,int> insertion (hash, -insert_position);
                    iter pos = std::lower_bound (hash_collection.begin(),
                                                 hash_collection.end(),
                                                 insertion);

                                // ok, found a data field with the same
                                // key. check whether we really hit a
                                // duplicate, i.e., whether the hash really
                                // was effective
                    bool duplicate = true;
                    if (pos != hash_collection.end() &&
                        pos->first == hash)
                      {
                        for (unsigned int d=0; d<dim; ++d)
                          for (unsigned int j=0; j<n_vectors; ++j)
                            if (std::fabs(data.const_jac[d][d][j]-
                                          cartesian[-pos->second].first[d][j])>
                                hash_value.scaling)
                              duplicate = false;
                      }
                    else
                      duplicate = false;

                                // if no duplicate, insert the data
                    if (duplicate == false)
                      {
                        hash_collection.insert (pos, insertion);
                        cartesian.push_back (new_entry);
                      }
                                // else, remember the position
                    else
                      insert_position = -pos->second;
                  }

                                // Constant Jacobian case. same strategy as
                                // before, but with other data fields
                else if (most_general_type == 1)
                  {
                    insert_position = linear.size();

                                // check whether everything is the same as on
                                // the previous cell
                    const unsigned int hash = hash_value.template operator()<dim,Number> (data.const_jac, false);
                    std::pair<unsigned int,int> insertion (hash, -insert_position);
                    iter pos = std::lower_bound (hash_collection.begin(),
                                                 hash_collection.end(),
                                                 insertion);

                                // ok, found a data field with the same
                                // key. check whether we really hit a
                                // duplicate
                    bool duplicate = true;
                    if (pos != hash_collection.end() &&
                        pos->first == hash)
                      {
                        for (unsigned int d=0; d<dim; ++d)
                          for (unsigned int e=0; e<dim; ++e)
                            for (unsigned int j=0; j<n_vectors; ++j)
                              if (std::fabs(data.const_jac[d][e][j]-
                                            linear[-pos->second].first[d][e][j])>
                                  hash_value.scaling)
                                duplicate = false;
                      }
                    else
                      duplicate = false;

                    if (duplicate == false)
                      {
                        hash_collection.insert (pos, insertion);
                        linear.push_back (std::pair<tensor2,vector_t>(data.const_jac,
                                                                      make_vectorized_array (Number(0.))));
                      }
                    else
                      insert_position = -pos->second;
                  }

                                // general cell case: first resize the data
                                // field to fit the new data. if we are here
                                // the first time, assume that there are many
                                // general cells to come, so reserve some
                                // memory in order to not have too many
                                // reallocations and memcpy's. The scheme used
                                // here involves at most one reallocation.
                else
                  {
                    Assert (most_general_type == 2, ExcInternalError());
                    insert_position = current_data.rowstart_jacobians.size();
                    if (current_data.rowstart_jacobians.size() == 0)
                      {
                        unsigned int reserve_size = (n_macro_cells-cell+1)/2;
                        current_data.rowstart_jacobians.reserve
                          (reserve_size);
                        reserve_size *= n_q_points;
                        current_data.jacobians.reserve (reserve_size);
                        if (update_flags & update_JxW_values)
                          current_data.JxW_values.reserve (reserve_size);
                        if (update_flags & update_jacobian_grads)
                          current_data.jacobians_grad_diag.reserve (reserve_size);
                        if (update_flags & update_jacobian_grads)
                          current_data.jacobians_grad_upper.reserve (reserve_size);
                      }
                  }

                cell_type[cell] = ((insert_position << n_cell_type_bits) +
                                   most_general_type);

              } // end if (my_q == 0)

                                // general cell case: now go through all
                                // quadrature points and collect the
                                // data. done for all different quadrature
                                // formulas, so do it outside the above loop.
            if (get_cell_type(cell) == 2)
              {
                const unsigned int previous_size =
                  current_data.jacobians.size();
                current_data.rowstart_jacobians.push_back (previous_size);
                if (update_flags & update_JxW_values)
                  {
                    AssertDimension (previous_size,
                                     current_data.JxW_values.size());
                  }
                if (update_flags & update_jacobian_grads)
                  {
                    AssertDimension (previous_size,
                                     current_data.jacobians_grad_diag.size());
                    AssertDimension (previous_size,
                                     current_data.jacobians_grad_upper.size());
                  }
                for (unsigned int q=0; q<n_q_points; ++q)
                  {
                    tensor2  &jac = data.general_jac[q];
                    tensor3  &jacobian_grad = data.general_jac_grad[q];
                    for (unsigned int j=0; j<n_vectors; ++j)
                      if (cell_t[j] < 2)
                        {
                          for (unsigned int d=0; d<dim; ++d)
                            for (unsigned int e=0; e<dim; ++e)
                              {
                                jac[d][e][j] = data.const_jac[d][e][j];
                                for (unsigned int f=0; f<dim; ++f)
                                  jacobian_grad[d][e][f][j] = 0.;
                              }
                        }

                    const vector_t det = determinant (jac);
                    current_data.jacobians.push_back (transpose(invert(jac)));
                    const tensor2 &inv_jac = current_data.jacobians.back();

                                // TODO: deal.II does not use abs on
                                // determinants. Is there an assumption
                                // somewhere that negative determinants are
                                // not allowed?
                    if (update_flags & update_JxW_values)
                      current_data.JxW_values.push_back
                        (std::abs(det) * current_data.quadrature_weights[fe_index][q]);

                    if (update_flags & update_jacobian_grads)
                      {
                                // for second derivatives on the real cell,
                                // need also the gradient of the inverse
                                // Jacobian J. This involves some calculus and
                                // is done vectorized. This is very cheap
                                // compared to what fe_values does (in early
                                // 2011). If L is the gradient of the jacobian
                                // on the unit cell, the gradient of the
                                // inverse is given by (multidimensional
                                // calculus) - J * (J * L) * J (the third J is
                                // because we need to transform the gradient L
                                // from the unit to the real cell, and then
                                // apply the inverse Jacobian). Compare this
                                // with 1D with j(x) = 1/k(phi(x)), where j =
                                // phi' is the inverse of the jacobian and k
                                // is the derivative of the jacobian on the
                                // unit cell. Then j' = phi' k'/k^2 = j k'
                                // j^2.

                                // compute: jac_grad = J*grad_unit(J^-1)
                        for (unsigned int d=0; d<dim; ++d)
                          for (unsigned int e=0; e<dim; ++e)
                            for (unsigned int f=0; f<dim; ++f)
                              {
                                jac_grad[f][e][d] = (inv_jac[f][0] *
                                                     jacobian_grad[d][e][0]);
                                for (unsigned int g=1; g<dim; ++g)
                                  jac_grad[f][e][d] += (inv_jac[f][g] *
                                                        jacobian_grad[d][e][g]);
                              }

                                // compute: transpose (-jac * jac_grad[d] * jac)
                        for (unsigned int d=0; d<dim; ++d)
                          for (unsigned int e=0; e<dim; ++e)
                            {
                              for (unsigned int f=0; f<dim; ++f)
                                {
                                  tmp[f] = vector_t();
                                  for (unsigned int g=0; g<dim; ++g)
                                    tmp[f] -= jac_grad[d][f][g] * inv_jac[g][e];
                                }

                                // needed for non-diagonal part of Jacobian
                                // grad
                              for (unsigned int f=0; f<dim; ++f)
                                {
                                  grad_jac_inv[f][d][e] = inv_jac[f][0] * tmp[0];
                                  for (unsigned int g=1; g<dim; ++g)
                                    grad_jac_inv[f][d][e] += inv_jac[f][g] * tmp[g];
                                }
                            }

                        {
                          vector_t grad_diag[dim][dim];
                          for (unsigned int d=0; d<dim; ++d)
                            for (unsigned int e=0; e<dim; ++e)
                              grad_diag[d][e] = grad_jac_inv[d][d][e];
                          current_data.jacobians_grad_diag.push_back
                            (Tensor<2,dim,vector_t>(grad_diag));
                        }

                                // sets upper-diagonal part of Jacobian
                        tensorUT grad_upper;
                        for (unsigned int d=0, count=0; d<dim; ++d)
                          for (unsigned int e=d+1; e<dim; ++e, ++count)
                            for (unsigned int f=0; f<dim; ++f)
                              grad_upper[count][f] = grad_jac_inv[d][e][f];
                        current_data.jacobians_grad_upper.push_back(grad_upper);
                      }
                  }
              }

            if (update_flags & update_quadrature_points)
              {
                                 // eventually we turn to the quadrature points
                                // that we can compress in case we have
                                // Cartesian cells. we also need to reorder
                                // them into arrays of vectorized data types.
                                // first go through the cells and find out how
                                // much memory we need to allocate for the
                                // quadrature points. We store n_q_points for
                                // all cells but Cartesian cells. For
                                // Cartesian cells, only need to store the
                                // values on a diagonal through the cell
                                // (n_q_points_1d). This will give (with some
                                // little indexing) the location of all
                                // quadrature points.
                const unsigned int old_size =
                  current_data.quadrature_points.size();
                current_data.rowstart_q_points[cell] = old_size;

                tensor1 quad_point;

                if (get_cell_type(cell) == 0)
                  {
                    current_data.quadrature_points.resize (old_size+
                                                           n_q_points_1d[fe_index]);
                    for (unsigned int q=0; q<n_q_points_1d[fe_index]; ++q)
                      for (unsigned int d=0; d<dim; ++d)
                        current_data.quadrature_points[old_size+q][d] =
                          data.quadrature_points[q*step_size_cartesian[fe_index]][d];
                  }
                else
                  {
                    current_data.quadrature_points.resize (old_size + n_q_points);
                    for (unsigned int q=0; q<n_q_points; ++q)
                      for (unsigned int d=0; d<dim; ++d)
                        current_data.quadrature_points[old_size+q][d] =
                          data.quadrature_points[q][d];
                  }
              }
          } // end for ( cell < n_macro_cells )
        current_data.rowstart_jacobians.push_back 
          (current_data.jacobians.size());
        current_data.rowstart_q_points[n_macro_cells] =
          current_data.quadrature_points.size();

                                // finally, need to invert and transpose the
                                // Jacobians in the cartesian and linear
                                // fields and compute the JxW value.
        if (my_q == 0)
          {
            for (unsigned int i=0; i<cartesian.size(); ++i)
              {
                vector_t det = cartesian[i].first[0];
                for (unsigned int d=1; d<dim; ++d)
                  det *= cartesian[i].first[d];
                for (unsigned int d=0; d<dim; ++d)
                  cartesian[i].first[d] = 1./cartesian[i].first[d];
                cartesian[i].second = std::abs(det);
              }
            for (unsigned int i=0; i<linear.size(); ++i)
              {
                vector_t det = determinant(linear[i].first);
                linear[i].first = transpose(invert(linear[i].first));
                linear[i].second = std::abs(det);
              }
          }
      }
  }



  template<int dim, typename Number>
  void
  MappingInfo<dim,Number>::evaluate_on_cell (const dealii::Triangulation<dim> &tria,
                                             const std::pair<unsigned int,unsigned int> *cells,
                                             const unsigned int  cell,
                                             const unsigned int  my_q,
                                             unsigned int (&cell_t_prev)[n_vectors],
                                             unsigned int (&cell_t)[n_vectors],
                                             FEValues<dim,dim> &fe_val,
                                             CellData          &data) const
  {
    const unsigned int n_q_points = fe_val.n_quadrature_points;
    const UpdateFlags update_flags = fe_val.get_update_flags();

                                // this should be the same value as
                                // used in HashValue::scaling (but we
                                // not have that field here)
    const double zero_tolerance_double = data.jac_size *
      std::numeric_limits<double>::epsilon() * 1024.;
    for (unsigned int j=0; j<n_vectors; ++j)
      {
        typename dealii::Triangulation<dim>::cell_iterator
          cell_it (&tria, cells[j].first, cells[j].second);
        fe_val.reinit(cell_it);
        cell_t[j] = numbers::invalid_unsigned_int;

                                // extract quadrature points and store them
                                // temporarily. if we have Cartesian cells, we
                                // can compress the indices
        if (update_flags & update_quadrature_points)
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const Point<dim> &point = fe_val.quadrature_point(q);
              for (unsigned int d=0; d<dim; ++d)
                data.quadrature_points[q][d][j] = point[d];
            }

                                // if this is not the first quadrature formula
                                // and we already have determined that this
                                // cell is either Cartesian or with constant
                                // Jacobian, we have nothing more to do.
        if (my_q > 0 && get_cell_type(cell) < 2)
          continue;

                                // first round: if the transformation is
                                // detected to be the same as on the old cell,
                                // we only need to copy over the data.
        if (fe_val.get_cell_similarity() == CellSimilarity::translation
            &&
            my_q == 0)
          {
            if (j==0)
              {
                Assert (cell>0, ExcInternalError());
                cell_t[j] = cell_t_prev[n_vectors-1];
              }
            else
              cell_t[j] = cell_t[j-1];
          }

        const DerivativeForm<1,dim,dim> & jac_0 = fe_val.jacobian(0);

        if (my_q == 0)
          {
                                // check whether the Jacobian is constant on
                                // this cell the first time we come around
                                // here
            if (cell_t[j] == numbers::invalid_unsigned_int)
              {
                bool jacobian_constant = true;
                for (unsigned int q=1; q<n_q_points; ++q)
                  {
                    const DerivativeForm<1,dim,dim> &jac = fe_val.jacobian(q);
                    for (unsigned int d=0; d<dim; ++d)
                      for (unsigned int e=0; e<dim; ++e)
                        if (std::fabs(jac_0[d][e]-jac[d][e]) >
                            zero_tolerance_double)
                          jacobian_constant = false;
                    if (jacobian_constant == false)
                      break;
                  }

                                // check whether the Jacobian is diagonal to
                                // machine accuracy
                bool cell_cartesian = jacobian_constant;
                for (unsigned int d=0; d<dim; ++d)
                  for (unsigned int e=0; e<dim; ++e)
                    if (d!=e)
                      if (std::fabs(jac_0[d][e]) >
                          zero_tolerance_double)
                        {
                          cell_cartesian=false;
                          break;
                        }

                                // in case we have only one quadrature point,
                                // we can have non-constant Jacobians, but we
                                // cannot detect it by comparison from one
                                // quadrature point to the next: in that case,
                                // need to look at second derivatives and see
                                // whether there are some non-zero entries
                                // (this is necessary since we determine the
                                // constness of the Jacobian for the first
                                // quadrature formula and might not look at
                                // them any more for the second, third
                                // quadrature formula
                if (cell_cartesian == false && n_q_points == 1)
                  {
                    Assert (update_flags & update_jacobian_grads,
                            ExcInternalError());
                    const DerivativeForm<1,dim,dim> &jac = fe_val.jacobian(0);
                    const DerivativeForm<2,dim,dim> &jacobian_grad =
                      fe_val.jacobian_grad(0);
                    for (unsigned int d=0; d<dim; ++d)
                      for (unsigned int e=0; e<dim; ++e)
                        for (unsigned int f=0; f<dim; ++f)
                          {
                            double jac_grad_comp = (jac[f][0] *
                                                    jacobian_grad[d][e][0]);
                            for (unsigned int g=1; g<dim; ++g)
                              jac_grad_comp += (jac[f][g] *
                                                jacobian_grad[d][e][g]);
                            if (std::fabs(jac_grad_comp) >
                                zero_tolerance_double)
                              jacobian_constant = false;
                          }
                  }
                                // set cell type
                if (cell_cartesian == true)
                  cell_t[j] = 0;
                else if (jacobian_constant == true)
                  cell_t[j] = 1;
                else
                  cell_t[j] = 2;
              }

                                // Cartesian cell
            if (cell_t[j] == 0)
              {
                                // set Jacobian into diagonal and clear
                                // off-diagonal part
                for (unsigned int d=0; d<dim; ++d)
                  {
                    data.const_jac[d][d][j] = jac_0[d][d];
                    for (unsigned int e=d+1; e<dim; ++e)
                      {
                        data.const_jac[d][e][j] = 0.;
                        data.const_jac[e][d][j] = 0.;
                      }
                  }
                continue;
              }

                                // cell with linear mapping
            else if (cell_t[j] == 1)
              {
                                // compress out very small values
                for (unsigned int d=0; d<dim; ++d)
                  for (unsigned int e=0; e<dim; ++e)
                    data.const_jac[d][e][j] =
                      std::fabs(jac_0[d][e]) < zero_tolerance_double ?
                      0 : jac_0[d][e];
                continue;
              }
          }

                                // general cell case

                                // go through all quadrature points and fill
                                // in the data into the temporary data
                                // structures with slots for the vectorized
                                // data types
        for (unsigned int q=0; q<n_q_points; ++q)
          {

                                // compress out very small numbers which are
                                // only noise. Then it is cleaner to use zero
                                // straight away (though it does not save any
                                // memory)
            const DerivativeForm<1,dim,dim> &jac = fe_val.jacobian(q);
            for (unsigned int d=0; d<dim; ++d)
              for (unsigned int e=0; e<dim; ++e)
                data.general_jac[q][d][e][j] =
                  std::fabs(jac[d][e]) < zero_tolerance_double ? 0. : jac[d][e];

                                // need to do some calculus based on the
                                // gradient of the Jacobian, in order to find
                                // the gradient of the inverse Jacobian which
                                // is needed in user code. however, we would
                                // like to perform that on vectorized data
                                // types instead of doubles or floats. to this
                                // end, copy the gradients first
            if (update_flags & update_jacobian_grads)
              {
                const DerivativeForm<2,dim,dim> &jacobian_grad = fe_val.jacobian_grad(q);
                for (unsigned int d=0; d<dim; ++d)
                  for (unsigned int e=0; e<dim; ++e)
                    for (unsigned int f=0; f<dim; ++f)
                      data.general_jac_grad[q][d][e][f][j] = jacobian_grad[d][e][f];
              }
          }
      } // end loop over all entries in vectorization (n_vectors cells)

                                // set information for next cell
    for (unsigned int j=0; j<n_vectors; ++j)
      cell_t_prev[j] = cell_t[j];
  }




  template <int dim, typename Number>
  std::size_t MappingInfo<dim,Number>::MappingInfoDependent::memory_consumption() const
  {
    std::size_t
      memory = MemoryConsumption::memory_consumption (jacobians);
    memory += MemoryConsumption::memory_consumption (JxW_values);
    memory += MemoryConsumption::memory_consumption (jacobians_grad_diag);
    memory += MemoryConsumption::memory_consumption (jacobians_grad_upper);
    memory += MemoryConsumption::memory_consumption (rowstart_q_points);
    memory += MemoryConsumption::memory_consumption (quadrature_points);
    memory += MemoryConsumption::memory_consumption (quadrature_formula);
    memory += MemoryConsumption::memory_consumption (quadrature_formula_faces);
    memory += MemoryConsumption::memory_consumption (quadrature_weights);
    memory += MemoryConsumption::memory_consumption (n_q_points);
    memory += MemoryConsumption::memory_consumption (n_q_points_face);
    memory += MemoryConsumption::memory_consumption (quad_index_conversion);
    return memory;
  }



  template <int dim, typename Number>
  std::size_t MappingInfo<dim,Number>::memory_consumption() const
  {
    std::size_t
      memory= MemoryConsumption::memory_consumption (mapping_data_gen);
    memory += MemoryConsumption::memory_consumption (linear);
    memory += MemoryConsumption::memory_consumption (cartesian);
    memory += MemoryConsumption::memory_consumption (cell_type);
    memory += sizeof (this);
    return memory;
  }



  template <int dim, typename Number>
  template <typename STREAM>
  void MappingInfo<dim,Number>::MappingInfoDependent::print_memory_consumption
  (STREAM         &out,
   const SizeInfo &size_info) const
  {
                                // print_mem involves global communication, so
                                // we can disable the check here only if no
                                // processor has any such data
#if DEAL_II_COMPILER_SUPPORTS_MPI
    unsigned int general_size_glob = 0, general_size_loc = jacobians.size();
    MPI_Allreduce (&general_size_loc, &general_size_glob, 1, MPI_UNSIGNED,
                   MPI_MAX, size_info.communicator);
#else
    unsigned int general_size_glob = jacobians.size();
#endif
    if (general_size_glob > 0)
      {
        out << "      Memory Jacobian data:          ";
        size_info.print_mem (out,
                             MemoryConsumption::memory_consumption (jacobians)
                             +
                             MemoryConsumption::memory_consumption (JxW_values));
        out << "      Memory second derivative data: ";
        size_info.print_mem (out,
                             MemoryConsumption::memory_consumption (jacobians_grad_diag)
                             +
                             MemoryConsumption::memory_consumption (jacobians_grad_upper));
      }

#if DEAL_II_COMPILER_SUPPORTS_MPI
    unsigned int quad_size_glob = 0, quad_size_loc = quadrature_points.size();
    MPI_Allreduce (&quad_size_loc, &quad_size_glob, 1, MPI_UNSIGNED,
                   MPI_MAX, size_info.communicator);
#else
    unsigned int quad_size_glob = quadrature_points.size();
#endif
    if (quad_size_glob > 0)
      {
        out << "      Memory quadrature points:      ";
        size_info.print_mem (out,
                             MemoryConsumption::memory_consumption (rowstart_q_points)
                             +
                             MemoryConsumption::memory_consumption (quadrature_points));
      }
  }



  template <int dim, typename Number>
  template <typename STREAM>
  void MappingInfo<dim,Number>::print_memory_consumption(STREAM &out,
                                                         const SizeInfo &size_info) const
  {
    out << "    Cell types:                      ";
    size_info.print_mem (out, MemoryConsumption::memory_consumption (cell_type));
    out << "    Memory transformations compr:    ";
    size_info.print_mem (out, MemoryConsumption::memory_consumption (linear) +
                         MemoryConsumption::memory_consumption (cartesian));
    for (unsigned int j=0; j<mapping_data_gen.size(); ++j)
      {
        out << "    Data component " << j << std::endl;
        mapping_data_gen[j].print_memory_consumption(out, size_info);
      }
  }

} // end of namespace MatrixFreeFunctions
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE