//----------------------------  q_collection.h  ----------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  q_collection.h  ----------------------------
#ifndef __deal2__q_collection_h
#define __deal2__q_collection_h

#include <base/config.h>
#include <base/subscriptor.h>
#include <base/quadrature.h>
#include <base/smartpointer.h>
#include <fe/fe.h>

#include <vector>
#include <boost/shared_ptr.hpp>

/**
 * This class implements a collection of quadrature objects used in the
 * hpDoFHandler(). It is thus to a hpDoFHandler() what a Quadrature is to a
 * DoFHandler. This collection stores pointers to the original quadrature
 * elements added to it. Hence the handling of the memory consumed by
 * the quadrature objects must be done by the user of the class.
 *
 * Although it is strongly recommended to supply an appropriate quadrature
 * for each finite element type used in a hp-computation, the QCollection
 * class implements a conversion constructor from a single quadrature.
 * Therefore it is possible to offer only a single quadrature to the
 * hpFEValues class instead of a QCollection. The reason for this
 * mechanism lies in the structure of the deal.II library. At many places
 * throughout the library pseudo quadrature rules are constructed to obtain
 * the values of the finite element function at certain points. An example
 * is the DataOut class. Due to this conversion constructor, these lines
 * of code continue to work without modifications.
 * 
 * @author Oliver Kayser-Herold, 2005
 */
template <int dim>
class QCollection : public Subscriptor
{
  public:
                                     /**
                                      * Default constructor. Initialises
				      * this QCollection.
				      */
    QCollection ();

                                     /**
                                      * Conversion constructor. This
				      * constructor creates a QCollection
				      * from a single quadrature rule. In
				      * the newly created QCollection, this
				      * quadrature is used for all active_fe
				      * indices.
				      */
    QCollection (const Quadrature<dim> &quadrature);

                                     /**
                                      * Returns the number of quadrature
				      * pointers stored in this object.
				      */
    unsigned int n_quadratures () const;

                                     /**
                                      * Returns a reference to the quadrature rule which
				      * was specified by the user for the
				      * active_fe_index which is provided
				      * as a parameter to this method.
				      */
    const Quadrature<dim> &
    get_quadrature (const unsigned int active_fe_index) const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

                                     /**
                                      * Adds a new quadrature rule to the QCollection.
				      * The quadrature rules have to be added in the order
				      * of the active_fe_indices. Thus the reference
				      * to the quadrature rule for active_fe_index 0
				      * has to be added first, followed by the
				      * quadrature rule for active_fe_index 1.
				      */
    unsigned int add_quadrature (const Quadrature<dim> &new_quadrature);

    DeclException0 (ExcNoQuadrature);
    
  private:
				     /**
				      * Upon construction of a <tt>QCollection</tt>
				      * the later functionality of the class is specified.
				      * Either it is a real collection, which provides
				      * different quadrature rules for each active_fe_index
				      * or its a "unreal" collection, which returns
				      * a the same quadrature rule for all active_fe_indices.
				      * This boolean remembers which type this object
				      * is.
				      */
    bool single_quadrature;

				     /**
				      * The real container, which stores pointers
				      * to the different quadrature objects.
				      */
    std::vector<SmartPointer<const Quadrature<dim> > > quadratures;
};



/* --------------- inline functions ------------------- */

template <int dim>
inline
unsigned int
QCollection<dim>::n_quadratures () const 
{
    return quadratures.size();
}


#endif
