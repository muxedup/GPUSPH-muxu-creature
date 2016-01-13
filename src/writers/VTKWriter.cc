/*  Copyright 2011-2013 Alexis Herault, Giuseppe Bilotta, Robert A. Dalrymple, Eugenio Rustico, Ciro Del Negro

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

#include <sstream>
#include <fstream>
#include <stdexcept>

#include "VTKWriter.h"
// GlobalData is required for writing the device index. With some order
// of inclusions, a forward declaration might be required
#include "GlobalData.h"

#include "vector_print.h"

using namespace std;

// TODO for the time being, we assume no more than 256 devices
// upgrade to UInt16 / ushort if it's ever needed

typedef unsigned char dev_idx_t;
static const char dev_idx_str[] = "UInt8";

VTKWriter::VTKWriter(const GlobalData *_gdata)
  : Writer(_gdata),
	m_planes_fname(),
	m_blockidx(-1)
{
	m_fname_sfx = ".vtu";

	string time_fname = open_data_file(m_timefile, "VTUinp", "", ".pvd");

	// Writing header of VTUinp.pvd file
	if (m_timefile) {
		m_timefile << "<?xml version='1.0'?>\n";
		m_timefile << "<VTKFile type='Collection' version='0.1'>\n";
		m_timefile << " <Collection>\n";
	}
}


VTKWriter::~VTKWriter()
{
	mark_timefile();
	m_timefile.close();
}

void VTKWriter::add_block(std::string const& blockname, std::string const& fname, double t)
{
	++m_blockidx;
	m_timefile << "  <DataSet timestep='" << t << "' group='" << m_blockidx <<
		"' name='" << blockname << "' file='" << fname << "'/>" << endl;
}

void VTKWriter::start_writing(double t)
{
	Writer::start_writing(t);

	m_blockidx = -1;

	const bool has_planes = gdata->s_hPlanes.size() > 0;

	if (has_planes) {
		if (m_planes_fname.size() == 0) {
			save_planes();
		}
		add_block("Planes", m_planes_fname, t);
	}
}

void VTKWriter::mark_written(double t)
{
	mark_timefile();

	Writer::mark_written(t);
}

/* Endianness check: (char*)&endian_int reads the first byte of the int,
 * which is 0 on big-endian machines, and 1 in little-endian machines */
static int endian_int=1;
static const char* endianness[2] = { "BigEndian", "LittleEndian" };

static float zeroes[4];

/* auxiliary functions to write data array entrypoints */
inline void
scalar_array(ofstream &out, const char *type, const char *name, size_t offset)
{
	out << "	<DataArray type='" << type << "' Name='" << name
		<< "' format='appended' offset='" << offset << "'/>" << endl;
}

inline void
vector_array(ofstream &out, const char *type, const char *name, uint dim, size_t offset)
{
	out << "	<DataArray type='" << type << "' Name='" << name
		<< "' NumberOfComponents='" << dim
		<< "' format='appended' offset='" << offset << "'/>" << endl;
}

inline void
vector_array(ofstream &out, const char *type, uint dim, size_t offset)
{
	out << "	<DataArray type='" << type
		<< "' NumberOfComponents='" << dim
		<< "' format='appended' offset='" << offset << "'/>" << endl;
}

// Binary dump a single variable of a given type
template<typename T>
inline void
write_var(ofstream &out, T const& var)
{
	out.write(reinterpret_cast<const char *>(&var), sizeof(T));
}

// Binary dump an array of variables of given type and size
template<typename T>
inline void
write_arr(ofstream &out, T const *var, size_t len)
{
	out.write(reinterpret_cast<const char *>(var), sizeof(T)*len);
}


void
VTKWriter::write(uint numParts, BufferList const& buffers, uint node_offset, double t, const bool testpoints)
{
	const double4 *pos = buffers.getData<BUFFER_POS_GLOBAL>();
	const hashKey *particleHash = buffers.getData<BUFFER_HASH>();
	const float4 *vel = buffers.getData<BUFFER_VEL>();
	const float4 *vol = buffers.getData<BUFFER_VOLUME>();
	const float *sigma = buffers.getData<BUFFER_SIGMA>();
	const particleinfo *info = buffers.getData<BUFFER_INFO>();
	const float3 *vort = buffers.getData<BUFFER_VORTICITY>();
	const float4 *normals = buffers.getData<BUFFER_NORMALS>();
	const float4 *gradGamma = buffers.getData<BUFFER_GRADGAMMA>();
	const float *tke = buffers.getData<BUFFER_TKE>();
	const float *eps = buffers.getData<BUFFER_EPSILON>();
	const float *turbvisc = buffers.getData<BUFFER_TURBVISC>();
	const float *spsturbvisc = buffers.getData<BUFFER_SPS_TURBVISC>();
	const float4 *eulervel = buffers.getData<BUFFER_EULERVEL>();
	const float *priv = buffers.getData<BUFFER_PRIVATE>();
	const vertexinfo *vertices = buffers.getData<BUFFER_VERTICES>();

	string filename;

	ofstream fid;
	filename = open_data_file(fid, "PART", current_filenum());

	// Header
	//====================================================================================
	fid << "<?xml version='1.0'?>" << endl;
	fid << "<VTKFile type='UnstructuredGrid'  version='0.1'  byte_order='" <<
		endianness[*(char*)&endian_int & 1] << "'>" << endl;
	fid << " <UnstructuredGrid>" << endl;
	fid << "  <Piece NumberOfPoints='" << numParts << "' NumberOfCells='" << numParts << "'>" << endl;
	fid << "   <PointData Scalars='Pressure' Vectors='Velocity'>" << endl;

	size_t offset = 0;

	// pressure
	scalar_array(fid, "Float32", "Pressure", offset);
	offset += sizeof(float)*numParts+sizeof(int);

	// density
	scalar_array(fid, "Float32", "Density", offset);
	offset += sizeof(float)*numParts+sizeof(int);

	// mass
	scalar_array(fid, "Float32", "Mass", offset);
	offset += sizeof(float)*numParts+sizeof(int);

	// gamma
	if (gradGamma) {
		scalar_array(fid, "Float32", "Gamma", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	// turbulent kinetic energy
	if (tke) {
		scalar_array(fid, "Float32", "TKE", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	// turbulent epsilon
	if (eps) {
		scalar_array(fid, "Float32", "Epsilon", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	// eddy viscosity
	if (turbvisc) {
		scalar_array(fid, "Float32", "Eddy viscosity", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	// SPS eddy viscosity
	if (spsturbvisc) {
		scalar_array(fid, "Float32", "SPS turbulent viscosity", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	/* Fluid number is only included if there are more than 1 */
	const bool write_fluid_num = (gdata->problem->physparams()->numFluids() > 1);

	/* Object number is only included if there are any */
	// TODO a better way would be for GPUSPH to expose the highest
	// object number ever associated with any particle, so that we
	// could check that
	const bool write_part_obj = (gdata->problem->simparams()->numbodies > 0);

	// particle info
	if (info) {
		scalar_array(fid, "UInt16", "Part type+flags", offset);
		offset += sizeof(ushort)*numParts+sizeof(int);
		// fluid number
		if (write_fluid_num) {
			// Limit to 256 fluids
			scalar_array(fid, "UInt8", "Fluid number", offset);
			offset += sizeof(uchar)*numParts+sizeof(int);
		}
		// object number
		if (write_part_obj) {
			// TODO UInt16 or UInt8 based on number of objects
			scalar_array(fid, "UInt16", "Part object", offset);
			offset += sizeof(ushort)*numParts+sizeof(int);
		}
		scalar_array(fid, "UInt32", "Part id", offset);
		offset += sizeof(uint)*numParts+sizeof(int);
	}

	if (vertices) {
		vector_array(fid, "UInt32", "Vertices", 4, offset);
		offset += sizeof(uint)*4*numParts+sizeof(int);
	}

	// device index
	if (MULTI_DEVICE) {
		scalar_array(fid, dev_idx_str, "DeviceIndex", offset);
		offset += sizeof(dev_idx_t)*numParts+sizeof(int);
	}

	// cell index
	scalar_array(fid, "UInt32", "CellIndex", offset);
	offset += sizeof(uint)*numParts+sizeof(int);

	// velocity
	vector_array(fid, "Float32", "Velocity", 3, offset);
	offset += sizeof(float)*3*numParts+sizeof(int);

	if (eulervel) {
		// Eulerian velocity
		vector_array(fid, "Float32", "Eulerian velocity", 3, offset);
		offset += sizeof(float)*3*numParts+sizeof(int);

		// Eulerian density
		scalar_array(fid, "Float32", "Eulerian density", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	// gradient gamma
	if (gradGamma) {
		vector_array(fid, "Float32", "Gradient Gamma", 3, offset);
		offset += sizeof(float)*3*numParts+sizeof(int);
	}

	// vorticity
	if (vort) {
		vector_array(fid, "Float32", "Vorticity", 3, offset);
		offset += sizeof(float)*3*numParts+sizeof(int);
	}

	// normals
	if (normals) {
		vector_array(fid, "Float32", "Normals", 3, offset);
		offset += sizeof(float)*3*numParts+sizeof(int);

		scalar_array(fid, "Float32", "Criteria", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	// private
	if (priv) {
		scalar_array(fid, "Float32", "Private", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	// volume
	if (vol) {
		vector_array(fid, "Float32", "Volume", 4, offset);
		offset += sizeof(float)*4*numParts+sizeof(int);
	}

	// sigma
	if (sigma) {
		scalar_array(fid, "Float32", "Sigma", offset);
		offset += sizeof(float)*numParts+sizeof(int);
	}

	fid << "   </PointData>" << endl;

	// position
	fid << "   <Points>" << endl;
	vector_array(fid, "Float64", 3, offset);
	offset += sizeof(double)*3*numParts+sizeof(int);
	fid << "   </Points>" << endl;

	// Cells data
	fid << "   <Cells>" << endl;
	scalar_array(fid, "Int32", "connectivity", offset);
	offset += sizeof(uint)*numParts+sizeof(int);
	scalar_array(fid, "Int32", "offsets", offset);
	offset += sizeof(uint)*numParts+sizeof(int);
	scalar_array(fid, "UInt8", "types", offset);
	offset += sizeof(uchar)*numParts+sizeof(int);
	fid << "   </Cells>" << endl;
	fid << "  </Piece>" << endl;

	fid << " </UnstructuredGrid>" << endl;
	fid << " <AppendedData encoding='raw'>\n_";
	//====================================================================================

	int numbytes=sizeof(float)*numParts;

	// pressure
	write_var(fid, numbytes);
	for (uint i=node_offset; i < node_offset + numParts; i++) {
		float value = 0.0;
		if (TESTPOINT(info[i]))
			value = vel[i].w;
		else
			value = m_problem->pressure(vel[i].w, fluid_num(info[i]));
		write_var(fid, value);
	}

	// density
	write_var(fid, numbytes);
	for (uint i=node_offset; i < node_offset + numParts; i++) {
		float value = 0.0;
		if (TESTPOINT(info[i]))
			// TODO FIXME: Testpoints compute pressure only
			// In the future we would like to have a density here
			// but this needs to be done correctly for multifluids
			value = NAN;
		else
			value = vel[i].w;
		write_var(fid, value);
	}

	// mass
	write_var(fid, numbytes);
	for (uint i=node_offset; i < node_offset + numParts; i++) {
		float value = pos[i].w;
		write_var(fid, value);
	}

	// gamma
	if (gradGamma) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float value = gradGamma[i].w;
			write_var(fid, value);
		}
	}

	// turbulent kinetic energy
	if (tke) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float value = tke[i];
			write_var(fid, value);
		}
	}

	// turbulent epsilon
	if (eps) {
		write_var(fid, numbytes);
		for (uint i=0; i < numParts; i++) {
			float value = eps[i];
			write_var(fid, value);
		}
	}

	// eddy viscosity
	if (turbvisc) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float value = turbvisc[i];
			write_var(fid, value);
		}
	}

	// SPS turbulent viscosity
	if (spsturbvisc) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float value = spsturbvisc[i];
			write_var(fid, value);
		}
	}

	// particle info
	if (info) {
		// type + flags
		numbytes=sizeof(ushort)*numParts;
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			ushort value = type(info[i]);
			write_var(fid, value);
		}

		// fluid number
		if (write_fluid_num) {
			numbytes=sizeof(uchar)*numParts;
			write_var(fid, numbytes);
			for (uint i=node_offset; i < node_offset + numParts; i++) {
				uchar value = fluid_num(info[i]);
				write_var(fid, value);
			}
		}

		if (write_part_obj) {
			numbytes=sizeof(ushort)*numParts;
			write_var(fid, numbytes);
			for (uint i=node_offset; i < node_offset + numParts; i++) {
				ushort value = object(info[i]);
				write_var(fid, value);
			}
		}

		// id
		numbytes=sizeof(uint)*numParts;
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			uint value = id(info[i]);
			write_var(fid, value);
		}
	}

	// vertices
	if (vertices) {
		numbytes = sizeof(uint)*4*numParts;
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			uint *value = (uint*)(vertices + i);
			write_arr(fid, value, 4);
		}
	}

	// device index
	if (MULTI_DEVICE) {
		numbytes = sizeof(dev_idx_t)*numParts;
		write_var(fid, numbytes);
		// The previous way was to compute the theoretical containing cell solely according on the particle position. This, however,
		// was inconsistent with the actual particle distribution among the devices, since one particle can be physically out of the
		// containing cell until next calchash/reorder.
		// The current policy is: just list the particles according to how the global array is partitioned. In other words, we rely
		// on the particle index to understad which device downloaded the particle data.
		for (uint d = 0; d < gdata->devices; d++) {
			// compute the global device ID for each device
			dev_idx_t value = gdata->GLOBAL_DEVICE_ID(gdata->mpi_rank, d);
			// write one for each particle (no need for the "absolute" particle index)
			for (uint p = 0; p < gdata->s_hPartsPerDevice[d]; p++)
				write_var(fid, value);
		}
		// There two alternate policies: 1. use particle hash or 2. compute belonging device.
		// To use the particle hash, instead of just relying on the particle index, use the following code:
		/*
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			uint value = gdata->s_hDeviceMap[ cellHashFromParticleHash(particleHash[i]) ];
			write_var(fid, value);
		}
		*/
		// This should be equivalent to the current "listing" approach. If for any reason (e.g. debug) one needs to write the
		// device index according to the current spatial position, it is enough to compute the particle hash from its position
		// instead of reading it from the particlehash array. Please note that this would reflect the spatial split but not the
		// actual assignments: until the next calchash is performed, one particle remains in the containing device even if it
		// it is slightly outside the domain.
	}

	// linearized cell index (NOTE: particles might be slightly off the belonging cell)
	numbytes = sizeof(uint)*numParts;
	write_var(fid, numbytes);
	for (uint i=node_offset; i < node_offset + numParts; i++) {
		uint value = cellHashFromParticleHash( particleHash[i] );
		write_var(fid, value);
	}

	numbytes=sizeof(float)*3*numParts;

	// velocity
	write_var(fid, numbytes);
	for (uint i=node_offset; i < node_offset + numParts; i++) {
		float *value = zeroes;
		//if (FLUID(info[i]) || TESTPOINTS(info[i]))
			value = (float*)(vel + i);
		write_arr(fid, value, 3);
	}

	if (eulervel) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float *value = zeroes;
			value = (float*)(eulervel + i);
			write_arr(fid, value, 3);
		}

		numbytes=sizeof(float)*numParts;

		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float value = eulervel[i].w;
			write_var(fid, value);
		}

		numbytes=sizeof(float)*3*numParts;
	}

	// gradient gamma
	if (gradGamma) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float *value = zeroes;
			value = (float*)(gradGamma + i);
			write_arr(fid, value, 3);
		}
	}

	// vorticity
	if (vort) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float *value = zeroes;
			if (FLUID(info[i])) {
				value = (float*)(vort + i);
			}
			write_arr(fid, value, 3);
		}
	}

	// normals
	if (normals) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float *value = zeroes;
			if (FLUID(info[i])) {
				value = (float*)(normals + i);
			}
			write_arr(fid, value, 3);
		}

		numbytes=sizeof(float)*numParts;
		// criteria
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float value = 0;
			if (FLUID(info[i]))
				value = normals[i].w;
			write_var(fid, value);
		}
	}

	numbytes=sizeof(float)*numParts;

	// private
	if (priv) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float value = priv[i];
			write_var(fid, value);
		}
	}

	numbytes=sizeof(float)*numParts*4;

	// volume
	if (vol) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float *value = (float*)(vol + i);
			write_arr(fid, value, 4);
		}
	}

	numbytes=sizeof(float)*numParts;

	// sigma
	if (sigma) {
		write_var(fid, numbytes);
		for (uint i=node_offset; i < node_offset + numParts; i++) {
			float value = sigma[i];
			write_var(fid, value);
		}
	}

	numbytes=sizeof(double)*3*numParts;

	// position
	write_var(fid, numbytes);
	for (uint i=node_offset; i < node_offset + numParts; i++) {
		double *value = (double*)(pos + i);
		write_arr(fid, value, 3);
	}

	numbytes=sizeof(int)*numParts;
	// connectivity
	write_var(fid, numbytes);
	for (uint i=0; i < numParts; i++) {
		uint value = i;
		write_var(fid, value);
	}
	// offsets
	write_var(fid, numbytes);
	for (uint i=0; i < numParts; i++) {
		uint value = i+1;
		write_var(fid, value);
	}

	// types (currently all cells type=1, single vertex, the particle)
	numbytes=sizeof(uchar)*numParts;
	write_var(fid, numbytes);
	for (uint i=0; i < numParts; i++) {
		uchar value = 1;
		write_var(fid, value);
	}

	fid << " </AppendedData>" << endl;
	fid << "</VTKFile>" << endl;

	add_block("Particles", filename, t);

}

void
VTKWriter::write_WaveGage(double t, GageList const& gage)
{
	ofstream fp;
	string filename = open_data_file(fp, "WaveGage", current_filenum());

	size_t num = gage.size();

	// For gages without points, z will be NaN, and we'll set
	// it to match the lowest world coordinate
	const double worldBottom = gdata->worldOrigin.z;

	// Header
	fp << "<?xml version='1.0'?>" << endl;
	fp << "<VTKFile type='UnstructuredGrid'  version='0.1'  byte_order='" <<
		endianness[*(char*)&endian_int & 1] << "'>" << endl;
	fp << " <UnstructuredGrid>" << endl;
	fp << "  <Piece NumberOfPoints='" << num << "' NumberOfCells='" << num << "'>" << endl;

	//Writing Position
	fp << "   <Points>" << endl;
	fp << "	<DataArray type='Float32' NumberOfComponents='3' format='ascii'>" << endl;
	for (size_t i=0; i <  num; i++)
		fp << gage[i].x << "\t" << gage[i].y << "\t" <<
			(isfinite(gage[i].z) ? gage[i].z : worldBottom) << "\t";
	fp << endl;
	fp << "	</DataArray>" << endl;
	fp << "   </Points>" << endl;

	// Cells data
	fp << "   <Cells>" << endl;
	fp << "	<DataArray type='Int32' Name='connectivity' format='ascii'>" << endl;
	for (size_t i = 0; i < num; i++)
		fp << i << "\t" ;
	fp << endl;
	fp << "	</DataArray>" << endl;
	fp << "" << endl;

	fp << "	<DataArray type='Int32' Name='offsets' format='ascii'>" << endl;
	for (size_t i = 0; i < num; i++)
		fp << (i+1) << "\t" ;
	fp << endl;
	fp << "	</DataArray>" << endl;

	fp << "" << endl;
	fp << "	<DataArray type='Int32' Name='types' format='ascii'>" << endl;
	for (size_t i = 0; i < num; i++)
		fp << 1 << "\t" ;
	fp << endl;
	fp << "	</DataArray>" << endl;

	fp << "   </Cells>" << endl;

	fp << "  </Piece>" << endl;
	fp << " </UnstructuredGrid>" << endl;
	fp << "</VTKFile>" <<endl;

	fp.close();

	add_block("WaveGages", filename, t);
}

static inline void chomp(double3 &pt, double eps=FLT_EPSILON)
{
		if (fabs(pt.x) < eps)
			pt.x = 0;
		if (fabs(pt.y) < eps)
			pt.y = 0;
		if (fabs(pt.z) < eps)
			pt.z = 0;
}

// check that pt is between inf and sup, with FLT_EPSILON relative tolerange
static inline bool bound(float pt, float inf, float sup)
{
	// when inf or sup is zero, the tolerance must be absolute, not relative
	// Also note the use of absolue value to ensure the limits are expanded
	// in the right direction
	const float lower = inf ? inf - FLT_EPSILON*fabs(inf) : -FLT_EPSILON;
	const float upper = sup ? sup + FLT_EPSILON*fabs(sup) : FLT_EPSILON;
	return (pt > lower) && (pt < upper);
}

void
VTKWriter::save_planes()
{
	ofstream fp;
	m_planes_fname = open_data_file(fp, "PLANES");

	fp << set_vector_fmt(" ");

	PlaneList const& planes = gdata->s_hPlanes;
	const double3 wo = gdata->problem->get_worldorigin();
	const double3 ow = wo + gdata->problem->get_worldsize();

	typedef std::vector<std::pair<double4, int> > CheckList;
	typedef std::vector<double3> CoordList;

	// We want to find the intersection of the planes defined in the boundary
	// with the bounding box of the plane (wo to ow). We do this by finding the intersection
	// with each pair of planes of the bounding box. The CheckList is composed of such pairs,
	// ordered such that the intersections are returned in sequence (otherwise the resulting
	// planes in the VTK would come out butterfly-shaped.
	// The number associated with each pair of planes is the index of the coordinate that must
	// be found by the intersection.
	CheckList checks;

	checks.push_back(make_pair(
			make_double4(wo.x, wo.y, 0, 1), 2));
	checks.push_back(make_pair(
			make_double4(wo.x, 0, wo.z, 1), 1));
	checks.push_back(make_pair(
			make_double4(wo.x, ow.y, 0, 1), 2));
	checks.push_back(make_pair(
			make_double4(wo.x, 0, ow.z, 1), 1));

	checks.push_back(make_pair(
			make_double4(ow.x, ow.y, 0, 1), 2));
	checks.push_back(make_pair(
			make_double4(ow.x, 0, ow.z, 1), 1));

	checks.push_back(make_pair(
			make_double4(ow.x, wo.y, 0, 1), 2));
	checks.push_back(make_pair(
			make_double4(0, wo.y, wo.z, 1), 0));
	checks.push_back(make_pair(
			make_double4(0, wo.y, ow.z, 1), 0));

	checks.push_back(make_pair(
			make_double4(0, ow.y, ow.z, 1), 0));

	checks.push_back(make_pair(
			make_double4(ow.x, 0, wo.z, 1), 1));
	checks.push_back(make_pair(
			make_double4(0, ow.y, wo.z, 1), 0));

	CoordList centers;
	CoordList normals;
	std::vector< CoordList > all_intersections;

	// we will store one point per plane (center)
	// followed by the intersections for each plane with the domain bounding box
	size_t npoints = planes.size();

	// find the intersection of each plane with the domain bounding box
	PlaneList::const_iterator plane(planes.begin());
	for (; plane != planes.end(); ++plane) {
		centers.push_back(gdata->calcGlobalPosOffset(plane->gridPos, plane->pos) + wo);
		double3 &cpos = centers.back();
		chomp(cpos);

		normals.push_back(make_double3(plane->normal));
		chomp(normals.back());
		double3 const& normal = normals.back();

		double4 implicit = make_double4(normal, -dot(cpos, normal));

#if DEBUG_VTK_PLANES
		cout << "plane through " << cpos << " normal " << normal << endl;
		cout << "\timplicit " << implicit << endl;
#endif

		all_intersections.push_back( std::vector<double3>() );

		std::vector<double3> & intersections = all_intersections.back();

		CheckList::const_iterator check(checks.begin());
		for (; check != checks.end(); ++check) {
			const double4 &ref = check->first;
			const int coord = check->second;
			double3 pt = make_double3(ref);
			switch (coord) {
			case 0:
				if (!normal.x) continue;
				pt.x = -dot(implicit, ref)/normal.x;
				if (!bound(pt.x, wo.x, ow.x)) continue;
				break;
			case 1:
				if (!normal.y) continue;
				pt.y = -dot(implicit, ref)/normal.y;
				if (!bound(pt.y, wo.y, ow.y)) continue;
				break;
			case 2:
				if (!normal.z) continue;
				pt.z = -dot(implicit, ref)/normal.z;
				if (!bound(pt.z, wo.z, ow.z)) continue;
				break;
			}
			chomp(pt);
			intersections.push_back(pt);
#if DEBUG_VTK_PLANES
			cout << "\t(" << (check-checks.begin()) << ")" << endl;
			cout << "\tcheck " << ref << " from " << coord << endl;
			cout << "\t\tpoint " << intersections.back() << endl;
#endif
		}
		npoints += intersections.size();
	}

	size_t offset = 0;

	fp << "<?xml version='1.0'?>" << endl;
	fp << "<VTKFile type='UnstructuredGrid'  version='0.1'  byte_order='" <<
		endianness[*(char*)&endian_int & 1] << "'>" << endl;
	fp << " <UnstructuredGrid>" << endl;
	fp << "  <Piece NumberOfPoints='" << npoints
		<< "' NumberOfCells='" << planes.size() << " '>" << endl;

	fp << "   <Points>" << endl;

	fp << "<DataArray type='Float64' NumberOfComponents='3'>" << endl;

	// intersection points
	for (std::vector<CoordList>::const_iterator pl(all_intersections.begin());
		pl < all_intersections.end(); ++pl) {
		CoordList const& pts = *pl;
		for (CoordList::const_iterator pt(pts.begin()); pt != pts.end(); ++pt)
			fp << *pt << endl;
	}

	// center points
	for (CoordList::const_iterator pt(centers.begin()); pt != centers.end(); ++pt)
		fp << *pt << endl;

	fp << "</DataArray>" << endl;

	fp << "   </Points>" << endl;

	fp << "   <Cells>" << endl;
	fp << "<DataArray type='Int32' Name='connectivity'>" << endl;
	// intersection points
	offset = 0;
	for (std::vector<CoordList>::const_iterator pl(all_intersections.begin());
		pl < all_intersections.end(); ++pl) {
		CoordList const& pts = *pl;
		for (int i = 0; i < pts.size(); ++i) {
			fp << " " << offset + i;
		}
		offset += pts.size();
		fp << endl;
	}
	fp << "</DataArray>" << endl;
	fp << "<DataArray type='Int32' Name='offsets'>" << endl;
	offset = 0;
	for (int i = 0; i < planes.size(); ++i) {
		offset += all_intersections[i].size();
		fp << offset << endl;
	}
	fp << "</DataArray>" << endl;
	fp << "<DataArray type='Int32' Name='types'>" << endl;
	for (int i = 0; i < planes.size(); ++i) {
		fp << 7 << " "; // POLYGON
	}
	fp << endl;
	fp << "</DataArray>" << endl;
	fp << "   </Cells>" << endl;

	fp << "   <PointData />" << endl;

	fp << "   <CellData Normals='Normals'>" << endl;
	fp << "<DataArray type='Float64' Name='Normals' NumberOfComponents='3'>" << endl;
	for (CoordList::const_iterator pt(normals.begin()); pt != normals.end(); ++pt)
		fp << *pt << endl;
	fp << "</DataArray>" << endl;
	fp << "   </CellData>" << endl;

	fp << "  </Piece>" << endl;
	fp << " </UnstructuredGrid>" << endl;
	fp << "</VTKFile>" <<endl;

	fp.close();
}

void
VTKWriter::mark_timefile()
{
	if (!m_timefile)
		return;
	// Mark the current position, close the XML, go back
	// to the marked position
	ofstream::pos_type mark = m_timefile.tellp();
	m_timefile << " </Collection>\n";
	m_timefile << "</VTKFile>" << endl;
	m_timefile.seekp(mark);
}
