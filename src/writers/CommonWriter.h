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
#ifndef H_COMMONWRITER_H
#define H_COMMONWRITER_H

/* The CommonWriter is a Writer dedicated to saving time series of common
 * aggregated data related to the particle system, such as its total energy,
 * the wavegages, etc.
 * It is independent of all other writers, active by default, and by default
 * writes a datapoint every time any other writer also writes (i.e. it doesn't
 * have a frequency of its own). It is of course possible to override this
 * behavior by setting its frequency to a finite value (e.g. 0 to write at
 * every timestep, a negative value to disable it, or a positive value to
 * write with a fixed frequency) in the problem.
 * The files created by the CommonWriter should all follow the following format:
 *
 * 0. tab separated
 * 1. a header on the first row (column titles)
 * 2. the simulation time should always be the first column
 * 3. one row per time point
 *
 * This makes it easy to import these files into common tools, and trivial
 * to plot them with e.g. gnuplot
 */

#include "Writer.h"

using std::ofstream;
using std::ostream;

class CommonWriter : public Writer
{
public:
	CommonWriter(const GlobalData *_gdata);
	~CommonWriter();

	void write(uint numParts, BufferList const& buffers, uint node_offset, double t, const bool testpoints);

	void write_energy(double t, float4 *energy);
	void write_WaveGage(double t, GageList const& gage);
	void write_objects(double t);
	void write_objectforces(double t, uint numobjects,
		const float3* computedforces, const float3* computedtorques,
		const float3* appliedforces, const float3* appliedtorques);
	void write_flux(double t, float *fluxes);

	bool need_write(double t) const;

private:
	/* Save a summary of phys_params, sim_params and options */
	void write_simparams(ostream &out);
	void write_physparams(ostream &out);
	void write_options(ostream &out);
	void write_summary();

	ofstream		m_energyfile;
	ofstream		m_WaveGagefile;
	ofstream		m_objectfile;
	ofstream		m_objectforcesfile;
	ofstream		m_fluxfile;

};
#endif

