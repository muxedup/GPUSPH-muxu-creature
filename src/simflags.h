/*  Copyright 2014 Alexis Herault, Giuseppe Bilotta, Robert A. Dalrymple, Eugenio Rustico, Ciro Del Negro

    Istituto Nazionale di Geofisica e Vulcanologia
        Sezione di Catania, Catania, Italy

    Università di Catania, Catania, Italy

    Johns Hopkins University, Baltimore, MD

    This file is part of GPUSPH.

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Set of boolean aspects of the simulation, to determine if
 * any of the features is enabled (XSPH, adaptive timestep, moving
 * boundaries, inlet/outlet, DEM, Ferrari correction, etc)
 */

#ifndef _SIMFLAGS_H
#define _SIMFLAGS_H

#include "common_types.h"

// no options
#define ENABLE_NONE			0UL

// adaptive timestepping
#define ENABLE_DTADAPT			1UL

// XSPH
#define ENABLE_XSPH				(ENABLE_DTADAPT << 1)

// DEM
#define ENABLE_DEM				(ENABLE_XSPH << 1)

// moving boundaries
#define ENABLE_MOVING_BODIES	(ENABLE_DEM << 1)

// inlet/outlet
#define ENABLE_INLET_OUTLET		(ENABLE_MOVING_BODIES << 1)

// water depth computation
#define ENABLE_WATER_DEPTH		(ENABLE_INLET_OUTLET << 1)

// TODO testpoints, rigid bodies, ferrari, ...

#endif