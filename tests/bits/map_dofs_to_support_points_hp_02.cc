//----------------------------------------------------------------------
//    $Id: map_dofs_to_support_points_hp_02.cc 24468 2011-09-29 15:44:06Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/dofs/dof_tools.h>

// check
//   DoFTools::
//   map_dofs_to_support_points(...);
// for the hp case with different finite elements
// on different cells.


using namespace std;

template <int dim>
void test ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  //define DoFhandler and FEs
  FE_Q<dim> fe1(1);
  FE_Q<dim> fe2(2);
  
  MappingQ<dim> mapping(1);
  hp::MappingCollection<dim> mapping_collection(mapping);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe1);
  fe_collection.push_back(fe2);

  hp::DoFHandler<dim> hp_dof_handler(triangulation);

  //distribute dofs
  hp_dof_handler.begin_active()->set_active_fe_index(1);
  hp_dof_handler.distribute_dofs(fe_collection);

  
  
  //now map the dofs to the support points and show them on the screen
  std::vector<Point<dim> > hp_map(hp_dof_handler.n_dofs());

  DoFTools::map_dofs_to_support_points(mapping_collection, hp_dof_handler, hp_map);

  // output the elements
  for (unsigned int i=0; i<hp_map.size(); i++)
  {
    deallog<< " Location of " << i<<" th DoF: "<<hp_map[i] << " | ";
  }
  deallog<<std::endl;

}

int main ()
{
  std::ofstream logfile("map_dofs_to_support_points_hp_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  return 0;
}