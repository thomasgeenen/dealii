//------------------  quadrature_points.cc  ------------------------
//    $Id$
//    Version: $Name$
//
//------------------  quadrature_points.cc  ------------------------


// this function tests the correctness of cached quadrature points

#include "../tests.h"

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vectors.h>

#include <fstream>
#include <iostream>

std::ofstream logfile("quadrature_points/output");


template <int dim, int fe_degree>
void test ()
{
  typedef double number;
  Triangulation<dim> tria;
  GridGenerator::hyper_ball (tria);
  static const HyperBallBoundary<dim> boundary;
  tria.set_boundary (0, boundary);
  tria.refine_global(5-dim);

  MappingQ<dim> mapping (4);
  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);
  deallog << "Testing " << fe.get_name() << std::endl;
  //std::cout << "Number of cells: " << tria.n_active_cells() << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  MatrixFree<dim,number> mf_data;
  {
    const QGauss<1> quad (fe_degree+1);
    typename MatrixFree<dim,number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim,number>::AdditionalData::none;
    data.mapping_update_flags = update_quadrature_points;
    mf_data.reinit (mapping, dof, constraints, quad, data);
  }

  double error_points = 0, abs_points = 0;
  const unsigned int n_cells = mf_data.get_size_info().n_macro_cells;
  FEEvaluation<dim,fe_degree+1> fe_eval (mf_data);
  FEValues<dim> fe_values (mapping, fe, mf_data.get_quad(),
			   update_quadrature_points);

  typedef VectorizedArray<double> vector_t;
  for (unsigned int cell=0; cell<n_cells; ++cell)
    {
      fe_eval.reinit(cell);
      for (unsigned int j=0; j<mf_data.n_components_filled(cell); ++j)
	{
	  fe_values.reinit (mf_data.get_cell_iterator(cell,j));
	  for (unsigned int q=0; q<fe_eval.n_q_points; ++q)
	    {
	      abs_points += fe_values.quadrature_point(q).norm();
	      for (unsigned int d=0; d<dim; ++d)
		error_points += std::fabs(fe_values.quadrature_point(q)[d]-
					  fe_eval.quadrature_point(q)[d][j]);
	    }
	}
    }

  deallog << "Norm of difference: " << error_points/abs_points
	  << std::endl << std::endl;
}


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog << std::setprecision (3);

  {
    deallog.threshold_double(1.e-12);
    deallog.push("2d");
    test<2,1>();
    test<2,2>();
    test<2,4>();
    deallog.pop();
    deallog.push("3d");
    test<3,1>();
    test<3,3>();
    deallog.pop();
  }
}
