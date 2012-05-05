//------------------  matrix_vector_05.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------  matrix_vector_05.cc  ------------------------


// this function tests the correctness of the implementation of matrix free
// matrix-vector products by comparing with the result of deal.II sparse
// matrix. The mesh uses a parallelogram mesh with hanging nodes (only cell
// type: 1 = linear).

#include "../tests.h"

std::ofstream logfile("matrix_vector_05/output");

#include "matrix_vector_common.h"


template <int dim, int fe_degree>
void test ()
{
  if (dim == 3)
    return;
  Triangulation<dim> tria;
  Tensor<2,dim> points;
  points[0][0] = 0.25;
  points[0][1] = 0.123;
  points[1][0] = 0.09983712334;
  points[1][1] = 0.314159265358979;
  GridGenerator::parallelogram (tria, points);
  typename Triangulation<dim>::active_cell_iterator
    cell = tria.begin_active (),
    endc = tria.end();
  for (; cell!=endc; ++cell)
    if (cell->center().norm()<1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active ();
  for (; cell!=endc; ++cell)
    if (cell->center().norm()<0.2)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(2);
  tria.begin(tria.n_levels()-1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active ();
  for (unsigned int i=0; i<5; ++i)
    {
      cell = tria.begin_active ();
      unsigned int counter = 0;
      for (; cell!=endc; ++cell, ++counter)
	if (counter % (7-i) == 0)
	  cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  do_test<dim, fe_degree, double> (dof, constraints);
}