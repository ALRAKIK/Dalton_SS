/* gen1int: compute one-electron integrals using rotational London atomic-orbitals
   Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen

   gen1int is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   gen1int is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with gen1int. If not, see <http://www.gnu.org/licenses/>.

   This file is the header file of partial geometric derivatives.

   2012-05-24, Bin Gao:
   * first version
*/

#if !defined(GEOM_PART_H)
#define GEOM_PART_H

/* QMatrix library for matrix operations */
#include "qmatrix.h"

/* information of partial geometric derivatives */
typedef struct {
  QInt num_cent;         /* number of differentiated centers */

  QInt *idx_atoms;       /* indices of atoms to be used as differentiated centers */
  QInt order_geo;        /* order of total geometric derivatives */
  QInt max_num_cent;     /* maximum number of differentiated centers */
  QInt num_paths;        /* total number of different allowed paths in N-ary tree */
  QInt num_unique_geo;   /* number of unique total geometric derivatives in the N-ary tree */
  QInt idx_path;         /* index of current path */
  Qint visit_height;     /* height of atom to visit */
  QInt *idx_node;        /* indices of the selected atom nodes */
  QInt *wt_node;         /* weights of the selected atom nodes */
  QInt *idx_cent;        /* indices of the generated differentiated centers */
  QInt *order_cent;      /* orders of geometric derivatives of the differentiated centers */
  QInt path_offset;      /* offset of unique derivatives of current path */
  QInt path_num_unique;  /* number of unique derivatives of current path */
  QInt path_num_redunt;  /* number of redundant derivatives of current path */
} PartGeom;

/* functions of partial geometric derivatives */
extern QErrorCode PartGeom
extern QErrorCode PartGeom
extern QErrorCode PartGeom
extern QErrorCode PartGeom

#endif
