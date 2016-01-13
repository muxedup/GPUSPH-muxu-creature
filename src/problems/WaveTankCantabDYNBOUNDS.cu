#include <cmath>
#include <iostream>
#include <stdexcept>

#include "WaveTankCantabDYNBOUNDS.h"
#include "particledefine.h"
#include "GlobalData.h"
#include "cudasimframework.cu"

#define OFFSET_X (-lx/2)
#define OFFSET_Y (-ly/2)
#define OFFSET_Z (-lz/2)

using namespace std;

WaveTankCantabDYNBOUNDS::WaveTankCantabDYNBOUNDS(GlobalData *_gdata) : Problem(_gdata)
{
	SETUP_FRAMEWORK(
		viscosity<SPSVISC>,
		boundary<DYN_BOUNDARY>,
		flags<ENABLE_DTADAPT>
	);
	
	// Size and origin of the simulation domain
	lx = 22.;
	ly = 32.0;
	lz = 3.5;
	
	// Data for problem setup
	slope_length = 10.2;  // horizontal extent of tank covered by sloping beach
	h_length = 10;
	height = 2.5;
	beta = .1974;
	
	H = 1;
	
	m_size = make_double3(lx, ly, lz);
	m_origin = make_double3(OFFSET_X, OFFSET_Y, OFFSET_Z);
	
	m_usePlanes = false;
	use_bottom_plane = false;
	
	set_deltap(0.15);
	cout << "deltap = " << m_deltap << endl;
	
	if (m_simparams->boundarytype == DYN_BOUNDARY && !m_usePlanes) {
		// number of layers
		dyn_layers = ceil(m_simparams->kernelradius*m_simparams->sfactor);
		// extra layers are one less (since other boundary types still have
		// one layer)
		double3 extra_offset = make_double3((dyn_layers-1)*m_deltap);
		m_origin -= extra_offset;
		m_size += 2*extra_offset;
	}
	
	// Wave maker parameters
	int numPeriods = 36;
	float paddlePeriod = 2.4;
	float incidentAngle = 0;
	m_simparams->tend = numPeriods * paddlePeriod;
	
	// SPH Stuff
	m_simparams->dt = 2e-5;
	m_simparams->dtadaptfactor = 0.2;
	m_simparams->buildneibsfreq = 10;
	
	// Physical parameters
	m_physparams->gravity = make_float3(0, 0, -9.81);
	const float g = length(m_physparams->gravity);
	const float r0 = m_deltap;
	m_physparams->r0 = r0;
	
	add_fluid(1000.f);
	set_kinematic_visc(0,1.0e-6);
	
    m_physparams->artvisccoeff= 0.2f;
	m_physparams->smagfactor = 0.12*0.12*m_deltap*m_deltap;
	m_physparams->kspsfactor = (2.0/3.0)*0.0066*m_deltap*m_deltap;
	m_physparams->epsartvisc = 0.01*m_simparams->slength*m_simparams->slength;
	
	add_writer(VTKWRITER, 0.1);
	m_name = "WaveTankCantabDYNBOUNDS";
	
}

WaveTankCantabDYNBOUNDS::~WaveTankCantabDYNBOUNDS(void) {
	release_memory();
}

void WaveTankCantabDYNBOUNDS::release_memory(void) {
	parts.clear();
	boundary_parts.clear();
}

int WaveTankCantabDYNBOUNDS::fill_parts() {
	cout << "filling particles" << endl;
	const float r0 = m_physparams->r0;
	
	boundary_parts.reserve(100);
	experiment_box = Cube(m_origin, m_size.x, m_size.y, m_size.z);
	
	experiment_box.FillIn(boundary_parts, r0, dyn_layers, false);
	return parts.size() + boundary_parts.size();
}

uint WaveTankCantabDYNBOUNDS::fill_planes()
{
	return (m_usePlanes ? 5 : 0);
}

void WaveTankCantabDYNBOUNDS::copy_planes(double4* planes)
{
	if (!m_usePlanes) return;

	planes[0] = make_double4(0, 0, 1.0, -m_origin.z);
	planes[1] = make_double4(0, 1.0, 0, -m_origin.x);
	planes[2] = make_double4(0, -1.0, 0, m_origin.x + ly);
	planes[3] = make_double4(1.0, 0, 0, -m_origin.y);
	planes[4] = make_double4(-1.0, 0, 0, m_origin.y + lz);
}

void WaveTankCantabDYNBOUNDS::copy_to_array(BufferList &buffers) {
	float4 *pos = buffers.getData<BUFFER_POS>();
	hashKey *hash = buffers.getData<BUFFER_HASH>();
	float4 *vel = buffers.getData<BUFFER_VEL>();
	particleinfo *info = buffers.getData<BUFFER_INFO>();
	vertexinfo *vertices = buffers.getData<BUFFER_VERTICES>();
	float4 *boundelm = buffers.getData<BUFFER_BOUNDELEMENTS>();

	cout << "Boundary parts: " << boundary_parts.size() << "\n";
	for (uint i = 0; i < boundary_parts.size(); i++) {
#if 1
		double water_column = H - boundary_parts[i](2);
		if (water_column < 0)
			water_column = 0;
		float rho = density(water_column, 0);
#else
		float rho = m_physparams->rho0[0];
#endif
		vel[i] = make_float4(0, 0, 0, rho);
		info[i] = make_particleinfo(PT_BOUNDARY, 0, i);
		calc_localpos_and_hash(boundary_parts[i], info[i], pos[i], hash[i]);
	}
	int j = boundary_parts.size();
	cout << "Boundary part mass: " << pos[j-1].w << "\n";
	cout << "Everything Uploaded" << endl;
}