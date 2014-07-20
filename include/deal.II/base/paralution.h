// ---------------------------------------------------------------------
// $Id: paralution.h 32926 2014-05-16 17:12:16Z turcksin $
//
// Copyright (C) 2014 by the deal.II authors
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

#ifndef __deal2__paralution_h
#define __deal2__paralution_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PARALUTION

#include <deal.II/base/multithread_info.h>

#include <paralution.hpp>

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  /**
   * A namespace for utility functions that abstract certain operations using
   * the Paralution library.
   */
  namespace Paralution
  {
    /**
     * A class that is used to initialize and to stop Paralution. If a program
     * uses Paralution one would typically just create an object of this type at
     * the beginninh of <code>main()</code>. The constructor of this class then
     * runs <code>paralution::init_paralution()</code>. At the end of the
     * program, the compiler will invoke the destructor of this object which in
     * turns calls <code>paralution::stop_paralution()</code>.
     */
    class Paralution_InitFinalize
    {
    public :
      /**
       * Constructor. Initialize Paralution by calling <tt>paralution::init_paralution</tt>
       * and set the number of threads used by Paralution to the given parameters.
       */
      Paralution_InitFinalize(const unsigned int max_num_threads);

      /**
       * Destructor. Calls <tt>paralution::stop_paralution</tt>.
       */
      ~Paralution_InitFinalize();

      /**
       * Set a device.
       */
      void set_device(const unsigned int device) const;

      /**
       * Set a specific gpu.
       */
      void set_gpu_cuda(const unsigned int gpu) const;

      /**
       * Set OpenCL compute units.
       */
      void set_ocl_compute_units(const unsigned int compute_unit) const;

      /**
       * Set a specific OpenCL platform.
       */
      void set_ocl_platform(const unsigned int platform) const;

      /**
       * Set a specific OpenCL platform and device.
       */
      void set_ocl(const unsigned int platform, const unsigned int device) const;

      /**
       * Set the number of OpenMP thread.
       */
      void set_omp_threads(const unsigned int n_threads) const;

      /**
       * Print information about the platform.
       */
      void info() const;
    };



// ------------------- inline functions --------------
    inline Paralution_InitFinalize::Paralution_InitFinalize(const unsigned int max_num_threads)
    {
      // Initialize paralution
      paralution::init_paralution();

      // Set the number of OpenMP threads
      set_omp_threads(max_num_threads);

      // Set the number of TBB threads
      multithread_info.set_thread_limit(max_num_threads);
    }


    inline Paralution_InitFinalize::~Paralution_InitFinalize()
    {
      paralution::stop_paralution();
    }



    inline void Paralution_InitFinalize::set_device(const unsigned int device) const
    {
      paralution::set_device_paralution(device);
    }



    inline void Paralution_InitFinalize::set_gpu_cuda(const unsigned int gpu) const
    {
      paralution::set_gpu_cuda_paralution(gpu);
    }



    inline void Paralution_InitFinalize::set_ocl_compute_units(const unsigned int compute_unit) const
    {
      paralution::set_ocl_compute_units_paralution(compute_unit);
    }



    inline void Paralution_InitFinalize::set_ocl_platform(const unsigned int platform) const
    {
      paralution::set_ocl_platform_paralution(platform);
    }



    inline void Paralution_InitFinalize::set_ocl(const unsigned int platform,
                                                 const unsigned int device) const
    {
      paralution::set_ocl_paralution(platform, device);
    }



    inline void Paralution_InitFinalize::set_omp_threads(const unsigned int n_threads) const
    {
      paralution::set_omp_threads_paralution(n_threads);
    }



    inline void Paralution_InitFinalize::info() const
    {
      paralution::info_paralution();
    }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PARALUTION

#endif
