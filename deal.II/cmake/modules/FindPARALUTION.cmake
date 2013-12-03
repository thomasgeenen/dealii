## ---------------------------------------------------------------------
## $Id: FindPARALUTION.cmake 31527 2013-11-03 09:58:45Z maier $
##
## Copyright (C) 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Try to find the PARALUTION library
#
# This module exports
#
#   PARALUTION_INCLUDE_DIRS
#   PARALUTION_LIBRARIES
#   PARALUTION_LINKER_FLAGS
#

SET_IF_EMPTY(PARALUTION_DIR "$ENV{PARALUTION_DIR}")

INCLUDE(FindPackageHandleStandardArgs)


FIND_PATH(PARALUTION_INCLUDE_DIR paralution.hpp
  HINTS
  ${PARALUTION_DIR}
  PATH_SUFFIXES
    inc
  )

FIND_LIBRARY(PARALUTION_LIBRARY
  NAMES paralution
  HINTS
  ${PARALUTION_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

SET(_output ${PARALUTION_LIBRARY})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PARALUTION DEFAULT_MSG
  _output # Cosmetic: Gives nice output
  PARALUTION_LIBRARY
  PARALUTION_INCLUDE_DIR
  )

MARK_AS_ADVANCED(
  PARALUTION_LIBRARY
  PARALUTION_INCLUDE_DIR
  )

IF(PARALUTION_FOUND)
  SET(PARALUTION_INCLUDE_DIRS
    ${PARALUTION_INCLUDE_DIR}
    )
  SET(PARALUTION_LIBRARIES
    ${PARALUTION_LIBRARY}
    )
  SET(PARALUTION_LINKER_FLAGS
    ${PARALUTION_LINKER_FLAGS}
    )

  MARK_AS_ADVANCED(PARALUTION_DIR)
ELSE()
  SET(PARALUTION_DIR "" CACHE PATH
    "An optional hint to a paralution directory"
    )
ENDIF()
