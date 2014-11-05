// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/trilinos_vector_base.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <cmath>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace Tpetra
  {

    namespace internal
    {
      VectorReference::operator TrilinosScalar () const
      {
        Assert (index < vector.size(),
                ExcIndexRange (index, 0, vector.size()));

        local_dof_index local_index = vector.vector->getMap().getLocalElement(index);

        Assert (local_index!=Teuchos::OrdinalTraits<local_dof_index>::invalid(),
                VectorBase::ExcAccessToNonLocalElement (index, vector.local_size(),
                                                        vector.vector->Map().MinMyGID(),
                                                        vector.vector->Map().MaxMyGID()));


        // Create the view
        Teuchos::ArrayRCP<TrilinosScalar> vector_view = vector.vector->get1dView();

        return vector_view[local_index];
      }
    }



    VectorBase::VectorBase ()
      :
      ghost_elements_readable  (false),
#ifdef DEAL_II_WITH_MPI
      vector(new Tpetra::MultiVector<TrilinosScalar,local_dof_index,global_dof_index,node>(
               Teuchos::RCP<Tpetra::Map<local_dof_index,global_dof_index,node> >(new
                   Tpetra::Map<local_dof_index,global_dof_index,node> (0,0,
                       Teuchos::RCP<Teuchos::MpiComm<int> >(MPI_COMM_SELF)))))
#else
      vector(new Tpetra::MultiVector<TrilinosScalar,local_dof_index,global_dof_index,node>(
               Teuchos::RCP<Tpetra::Map<local_dof_index,global_dof_index,node> >(new
                   Tpetra::Map<local_dof_index,global_dof_index,node> (0,0,
                       Teuchos::RCP<Teuchos::SerialComm<int> >()))))
#endif
    {}



    VectorBase::VectorBase (const VectorBase &v)
      :
      Subscriptor(),
      ghost_elements_readable (false),
      vector(new Tpetra::MultiVector<TrilinosScalar,local_dof_index,global_dof_index,
             node>(*v.vector)),
      nonlocal_vector(new Tpetra::MultiVector<TrilinosScalar,local_dof_index,global_dof_index,
                      node>(v.nonlocal_vector.getMap(),1))
    {}



    VectorBase::~VectorBase ()
    {}



    void
    VectorBase::clear ()
    {
      // When we clear the vector, reset the pointer and generate an empty
      // vector.
#ifdef DEAL_II_WITH_MPI
      Teuchos::<Tpetra::Map<local_dof_index,global_dof_index,node> > map(0,0,
          Teuchos::RCP<Teuchos::MpiComm<int> > (MPI_COMM_SELF));
#else
      Teuchos::<Tpetra::Map<local_dof_index,global_dof_index,node> > map(0,0,
          Teuchos::RCP<Teuchos::SerialComm<int> > ());
#endif

      ghost_elements_readable = false;
      vector.reset (new Tpetra::MultiVector<TrilinosScalar,local_dof_index,
                    global_dof_index,node>(map,1));

      nonlocal_vector.reset (new Tpetra::MultiVector<TrilinosScalar,local_dof_index,
                             global_dof_index,node>(map,1));
    }



    VectorBase &
    VectorBase::operator = (const VectorBase &v)
    {
      Assert (vector.get() != 0,
              ExcMessage("Vector is not constructed properly."));
      Assert (nonlocal_vector.get() != 0,
              ExcMessage("Vector is not constructed properly."));

      if (local_range() != v.local_range())
        {
          vector.reset (new Tpetra::MultiVector<local_dof_index,global_dof_index,
                        node>(*v.vector));
          nonlocal_vector.reset (new Tpetra::MultiVector<local_dof_index,global_dof_index,
                                 node>(*v.nonlocal_vector));
          ghost_elements_readable = v.ghost_elements_readable;
        }
      else
        {
          Assert (vector->getMap().isSameAs(v.vector->getMap()) == true,
                  ExcMessage ("The Tpetra maps in the assignment operator ="
                              " do not match, even though the local_range "
                              " seems to be the same. Check vector setup!"));
          Assert (nonlocal_vector->getMap().isSameAs(v.nonlocal_vector->getMap()) == true,
                  ExcMessage ("The Tpetra maps in the assignment operator ="
                              " do not match, even though the local_range "
                              " seems to be the same. Check vector setup!"));

          vector->update(1.0, *v.vector, 0.0);
          nonlocal_vector->update(1.0, *v.nonlocal_vector, 0.0);
        }

      return *this;
    }



    template <typename number>
    VectorBase &
    VectorBase::operator = (const ::dealii::Vector<number> &v)
    {
      Assert (size() == v.size(),
              ExcDimensionMismatch(size(), v.size()));
      // Create the view
      Teuchos::ArrayRCP<TrilinosScalar> vector_view = vector->get1dViewNonConst();

      // this is probably not very efficient
      // but works. in particular, we could do
      // better if we know that
      // number==TrilinosScalar because then we
      // could elide the copying of elements
      //
      // let's hope this isn't a
      // particularly frequent operation
      std::pair<size_type, size_type>
      local_range = this->local_range ();
      for (size_type i=local_range.first; i<local_range.second; ++i)
        vector_view[i-local_range.first] = v(i);

      zero_out_ghosts;

      return *this;
    }



    void
    VectorBase::compress (::dealii::VectorOperation::values action)
    {
      Tpetra::CombineMode mode;
      if (given_last_action==::dealii::VectorOperation::add)
        mode = Add;
      else if (given_last_action==::dealii::VectorOperation::insert)
        mode = Insert;

      Tpetra::Export<local_dof_index,global_dof_index,node> exporter(nonlocal_vector->getMap(),
          vector->getMap());
      vector->doExport(*nonlocal_vector, exporter, mode);
      zero_out_ghosts();
    }



    TrilinosScalar
    VectorBase::operator () (const size_type index) const
    {
      // Extract local indices in the vector.
      local_dof_index local_i = vector->getMap().getLocalElement(index);
      TrilinosScalar value = 0.;

      Assert (local_i != Teuchos::OrdinalTraits<local_dof_index>::invalid() ||
              ghost_elements_readable == true,
              ExcMessage("You tried to read a ghost element of this vector, "
                         "but it has not imported its ghost values."));
      if  (local_i != Teuchos::OrdinalTraits<local_dof_index>::invalid())
        {
          // Create the view
          Teuchos::ArrayRCP<TrilinosScalar> vector_view = vector->get1dView();

          value = vector_view[local_i];
        }
      else
        {
          local_dof_index nonlocal_i = nonlocal_vector->getMap().getLocalElement(index);

          Assert (local_i != Teuchos::OrdinalTraits<local_dof_index>::invalid(),
                  ExcMessage("You tried to access element " << index
                             << " of a distributed vector, but this element is not stored "
                             << "on the current processor."));

          // Create the view
          Teuchos::ArrayRCP<TrilinosScalar> nonlocal_vector_view = nonlocal_vector->get1dView();

          value = nonlocal_vector_view[local_i];
        }

      return value;
    }



    void
    VectorBase::add (const VectorBase &v,
                     const bool        allow_different_maps)
    {
      if (allow_different_maps == false)
        *this += v;
      else
        {
          AssertThrow (size() == v.size(),
                       ExcDimensionMismatch (size(), v.size()));

          Tpetra::Import<local_dof_index,global_dof_index,node> data_exchange (vector->getMap(),
              v.vector->getMap());
          vector->doImport(*v.vector, data_exchange, Tpetra::CombineMode::ADD);
          //TODO what to do with the nonlocal vector
        }
    }



    bool
    VectorBase::operator == (const VectorBase &v) const
    {
      Assert (size() == v.size(),
              ExcDimensionMismatch(size(), v.size()));
      if (local_size() != v.local_size())
        return false;

      // Create the views
      Teuchos::ArrayRCP<TrilinosScalar> vector_view = vector->get1dView();
      Teuchos::ArrayRCP<TrilinosScalar> v_view = v->get1dView();
      size_type i;
      for (i=0; i<local_size(); i++)
        if (v_view[i]!=vector_view[i])
          return false;

      return true;
    }



    bool
    VectorBase::operator != (const VectorBase &v) const
    {
      Assert (size() == v.size(),
              ExcDimensionMismatch(size(), v.size()));

      return (!(*this==v));
    }



    bool
    VectorBase::all_zero () const
    {
      // Create the view
      Teuchos::ArrayRCP<TrilinosScalar> vector_view = vector->get1dView();
      // get a representation of the vector and
      // loop over all the elements
      TrilinosScalar *start_ptr = vector_view[0];
      const TrilinosScalar *ptr  = start_ptr,
                            *eptr = start_ptr + local_size();
      unsigned int flag = 0;
      while (ptr != eptr)
        {
          if (*ptr != 0.)
            {
              flag = 1;
              break;
            }
          ++ptr;
        }

#ifdef DEAL_II_WITH_MPI
      // in parallel, check that the vector
      // is zero on _all_ processors.
      const Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = vector->getMap().getComm();
      unsigned int num_nonzero = Utilities::MPI::sum(flag, mpi_comm->getRawComm());
      return num_nonzero == 0;
#else
      return flag == 0;
#endif

    }



    bool
    VectorBase::is_non_negative () const
    {
#ifdef DEAL_II_WITH_MPI
      // if this vector is a parallel one, then
      // we need to communicate to determine
      // the answer to the current
      // function. this still has to be
      // implemented
      AssertThrow(local_size() == size(), ExcNotImplemented());
#endif
      // get a representation of the vector and
      // loop over all the elements
      Teuchos::ArrayRCP<TrilinosScalar> vector_view = vector.vector->get1dView();
      TrilinosScalar *start_ptr = &vector_view[0];

      // TODO: This
      // won't work in parallel like
      // this. Find out a better way to
      // this in that case.
      const TrilinosScalar *ptr  = start_ptr,
                            *eptr = start_ptr + size();
      bool flag = true;
      while (ptr != eptr)
        {
          if (*ptr < 0.0)
            {
              flag = false;
              break;
            }
          ++ptr;
        }

      return flag;
    }



    // TODO: up to now only local
    // data printed out! Find a
    // way to neatly output
    // distributed data...
    void
    VectorBase::print (const char *format) const
    {
      Assert (vector->getGlobalLength()!=0, ExcEmptyObject());

      Teuchos::ArrayRCP<TrilinosScalar> vector_view = vector.vector->get1dView();

      for (size_type j=0; j<size(); ++j)
        {
          double t = vector_view[j];

          if (format != 0)
            std::printf (format, t);
          else
            std::printf (" %5.2f", double(t));
        }
      std::printf ("\n");
    }



    void
    VectorBase::print (std::ostream      &out,
                       const unsigned int precision,
                       const bool         scientific,
                       const bool         across) const
    {
      AssertThrow (out, ExcIO());

      // get a representation of the
      // vector and loop over all
      // the elements TODO: up to
      // now only local data printed
      // out! Find a way to neatly
      // output distributed data...
      Teuchos::ArrayRCP<TrilinosScalar> vector_view = vector.vector->get1dView();
      TrilinosScalar *val = &vector_view[0];

      AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      out.precision (precision);
      if (scientific)
        out.setf (std::ios::scientific, std::ios::floatfield);
      else
        out.setf (std::ios::fixed, std::ios::floatfield);

      if (across)
        for (size_type i=0; i<size(); ++i)
          out << val[i] << ' ';
      else
        for (size_type i=0; i<size(); ++i)
          out << val[i] << std::endl;
      out << std::endl;

      // restore the representation
      // of the vector
      AssertThrow (out, ExcIO());
    }



    void
    VectorBase::swap (VectorBase &v)
    {
      std::swap(last_action, v.last_action);
      std::swap(compressed, v.compressed);
      std::swap(vector, v.vector);
    }



    std::size_t
    VectorBase::memory_consumption () const
    {
      //TODO[TH]: No accurate memory
      //consumption for Trilinos vectors
      //yet. This is a rough approximation with
      //one index and the value per local
      //entry.
      return sizeof(*this)
             + this->local_size()*( sizeof(double)+
                                    sizeof(TrilinosWrappers::types::int_type) );
    }

  }
} /* end of namespace TrilinosWrappers */


namespace TrilinosWrappers
{
#include "trilinos_vector_base.inst"
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
