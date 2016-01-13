/*  Copyright 2013 Alexis Herault, Giuseppe Bilotta, Robert A.
 	Dalrymple, Eugenio Rustico, Ciro Del Negro

	Conservatoire National des Arts et Metiers, Paris, France

	Istituto Nazionale di Geofisica e Vulcanologia,
    Sezione di Catania, Catania, Italy

    Universita di Catania, Catania, Italy

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
#include <stdexcept>
#include <string>
//#include <cmath>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Problem.h"
#include "vector_math.h"
#include "vector_print.h"
#include "utils.h"

// here we need the complete definition of the GlobalData struct
#include "GlobalData.h"

// COORD1, COORD2, COORD3
#include "linearization.h"

using namespace std;

Problem::Problem(GlobalData *_gdata) :
	m_problem_dir(_gdata->clOptions->dir),
	m_dem(NULL),
	m_size(make_double3(NAN, NAN, NAN)),
	m_origin(make_double3(NAN, NAN, NAN)),
	m_deltap(NAN),
	gdata(_gdata),
	m_options(_gdata->clOptions),
	m_physparams(new PhysParams()),
	m_simframework(NULL),
	m_bodies_storage(NULL)
{
}

bool
Problem::initialize()
{
	if (simparams()->gage.size() > 0 && !m_simframework->hasPostProcessEngine(SURFACE_DETECTION)) {
		printf("Wave gages present: force-enabling surface detection\n");
		m_simframework->addPostProcessEngine(SURFACE_DETECTION);
	}
	// run post-construction functions
	check_dt();
	check_maxneibsnum();
	calculateFerrariCoefficient();
	create_problem_dir();

	printf("Problem calling set grid params\n");
	set_grid_params();

	return true;
}

Problem::~Problem(void)
{
	delete [] m_bodies_storage;
	delete m_simframework;
	delete m_physparams;

}


/// Allocate storage required for the integration of the kinematic data
/// of moving bodies.
void
Problem::allocate_bodies_storage()
{
	const uint nbodies = simparams()->numbodies;

	if (nbodies) {
		// TODO: this should depend on the integration scheme
		m_bodies_storage = new KinematicData[nbodies];
	}
}

void
Problem::add_moving_body(Object* object, const MovingBodyType mbtype)
{
	// Moving bodies are put at the end of the bodies vector,
	// ODE bodies and moving bodies for which we want a force feedback
	// are put at the beginning of the bodies vector (ODE bodies first).
	// The reason behind this ordering is the way the forces on bodies
	// are reduced by a parallel prefix sum: all the bodies that require
	// force computing must have consecutive ids.
	const uint index = m_bodies.size();
	if (index >= MAX_BODIES)
		throw runtime_error ("Number of moving bodies superior to MAX_BODIES. Increase MAXBODIES\n");
	MovingBodyData *mbdata = new MovingBodyData;
	mbdata->index = index;
	mbdata->type = mbtype;
	mbdata->object = object;
	mbdata->kdata.crot = object->GetCenterOfGravity();
	mbdata->kdata.lvel = make_double3(0.0f);
	mbdata->kdata.avel = make_double3(0.0f);
	mbdata->kdata.orientation = object->GetOrientation();
	switch (mbdata->type) {
		case MB_ODE : {
			const dBodyID bodyid = object->m_ODEBody;
			mbdata->kdata.crot = make_double3(dBodyGetPosition(bodyid));
			mbdata->kdata.lvel = make_double3(dBodyGetLinearVel(bodyid));
			mbdata->kdata.avel = make_double3(dBodyGetAngularVel(bodyid));
			m_bodies.insert(m_bodies.begin() + simparams()->numODEbodies, mbdata);
			simparams()->numODEbodies++;
			simparams()->numforcesbodies++;
			break;
		}

		case MB_FORCES_MOVING:
			m_bodies.insert(m_bodies.begin() + simparams()->numforcesbodies, mbdata);
			simparams()->numforcesbodies++;
			break;

		case MB_MOVING:
			m_bodies.push_back(mbdata);
			break;
	}

	mbdata->initial_kdata = mbdata->kdata;

	simparams()->numbodies = m_bodies.size();
}


MovingBodyData *
Problem::get_mbdata(const uint index)
{
	if (index >= m_bodies.size()) {
		stringstream ss;
		ss << "get_body: body number " << index << " >= numbodies";
		throw runtime_error(ss.str());
	}
	for (vector<MovingBodyData *>::iterator it = m_bodies.begin() ; it != m_bodies.end(); ++it) {
		if ((*it)->index == index)
			return *it;
	}
	return NULL;
}


MovingBodyData *
Problem::get_mbdata(const Object* object)
{
	for (vector<MovingBodyData *>::iterator it = m_bodies.begin() ; it != m_bodies.end(); ++it) {
		if ((*it)->object == object)
			return *it;
	}
	throw runtime_error("get_body: invalid object\n");
	return NULL;
}

size_t
Problem::get_bodies_numparts(void)
{
	size_t total_parts = 0;
	for (vector<MovingBodyData *>::iterator it = m_bodies.begin() ; it != m_bodies.end(); ++it) {
		total_parts += (*it)->object->GetNumParts();
	}

	return total_parts;
}


size_t
Problem::get_forces_bodies_numparts(void)
{
	size_t total_parts = 0;
	for (vector<MovingBodyData *>::iterator it = m_bodies.begin() ; it != m_bodies.end(); ++it) {
		if ((*it)->type == MB_ODE || (*it)->type == MB_FORCES_MOVING)
			total_parts += (*it)->object->GetNumParts();
	}
	return total_parts;
}


size_t
Problem::get_body_numparts(const int index)
{
	return m_bodies[index]->object->GetNumParts();
}


size_t
Problem::get_body_numparts(const Object* object)
{
	return get_mbdata(object)->object->GetNumParts();
}

/*void
Problem::restore_ODE_body(const uint i, const float *gravity_center, const float *quaternion,
	const float *linvel, const float *angvel)
{
	Object *obj = m_ODE_bodies[i];
	dBodyID odeid = obj->m_ODEBody;

	// re-set the position, rotation and velocities in ODE
	dBodySetAngularVel(odeid, angvel[0], angvel[1], angvel[2]);
	dBodySetLinearVel(odeid, linvel[0], linvel[1], linvel[2]);
	dBodySetPosition(odeid, gravity_center[0], gravity_center[1], gravity_center[2]);

	dBodySetQuaternion(odeid, quaternion);

	// After setting the quaternion, ODE does a forced renormalization
	// that will slightly change the value of the quaternion (except in some
	// trivial cases). While the final result is within machine precision to
	// the set value, the (small) difference will propagate through the
	// simulation, resulting in differences. The following code can be used to
	// check the amount of absolute and relative error in the set quaternion:
#if 0
	dQuaternion rec;
	dQuaternion abs_err, rel_err;
	dBodyCopyQuaternion(odeid, rec);
	for (int i = 0; i < 4; ++i) {
		abs_err[i] = fabs(rec[i] - quaternion[i]);
		float normfactor = fabs(rec[i]+quaternion[i])/2;
		rel_err[i] = normfactor == 0 ? abs_err[i] : abs_err[i]/normfactor;
	}

	printf("object %u quaternion: recovered (%g, %g, %g, %g), was (%g, %g, %g, %g),\n"
		"\tdelta (%g, %g, %g, %g), rel err (%g, %g, %g, %g)\n",
		i, rec[0], rec[1], rec[2], rec[3],
		quaternion[0], quaternion[1], quaternion[2], quaternion[3],
		abs_err[0], abs_err[1], abs_err[2], abs_err[3],
		rel_err[0], rel_err[1], rel_err[2], rel_err[3]);
#endif
}*/


void
Problem::calc_grid_and_local_pos(double3 const& globalPos, int3 *gridPos, float3 *localPos) const
{
	int3 _gridPos = calc_grid_pos(globalPos);
	*gridPos = _gridPos;
	*localPos = make_float3(globalPos - m_origin -
		(make_double3(_gridPos) + 0.5)*m_cellsize);
}

void
Problem::get_bodies_cg(void)
{
	for (uint i = 0; i < simparams()->numbodies; i++) {
		calc_grid_and_local_pos(m_bodies[i]->kdata.crot,
			gdata->s_hRbCgGridPos + i,
			gdata->s_hRbCgPos + i);
	}
}


void
Problem::set_body_cg(const double3& crot, MovingBodyData* mbdata) {
	mbdata->kdata.crot = crot;

}


void
Problem::set_body_cg(const uint index, const double3& crot) {
	set_body_cg(crot, m_bodies[index]);
}


void
Problem::set_body_cg(const Object *object, const double3& crot) {
	set_body_cg(crot, get_mbdata(object));
}


void
Problem::set_body_linearvel(const double3& lvel, MovingBodyData* mbdata) {
	mbdata->kdata.lvel = lvel;

}


void
Problem::set_body_linearvel(const uint index, const double3& lvel) {
	set_body_linearvel(lvel, m_bodies[index]);
}


void
Problem::set_body_linearvel(const Object *object, const double3& lvel)
{
	set_body_linearvel(lvel, get_mbdata(object));
}


void
Problem::set_body_angularvel(const double3& avel, MovingBodyData* mbdata) {
	mbdata->kdata.avel = avel;

}

void
Problem::set_body_angularvel(const uint index, const double3& avel) {
	set_body_angularvel(avel, m_bodies[index]);
}


void
Problem::set_body_angularvel(const Object *object, const double3& avel)
{
	set_body_angularvel(avel, get_mbdata(object));
}


void
Problem::bodies_forces_callback(const double t0, const double t1, const uint step, float3 *forces, float3 *torques)
{ /* default does nothing */ }


void
Problem::post_timestep_callback(const double t)
{ /* default does nothing */ }


void
Problem::moving_bodies_callback(const uint index, Object* object, const double t0, const double t1,
		const float3& force, const float3& torque, const KinematicData& initial_kdata,
		KinematicData& kdata, double3& dx, EulerParameters& dr)
{ /* default does nothing */ }

// input: force, torque, step number, dt
// output: cg, trans, steprot (can be input uninitialized)
void
Problem::bodies_timestep(const float3 *forces, const float3 *torques, const int step,
		const double dt, const double t,
		int3 * & cgGridPos, float3 * & cgPos, float3 * & trans, float * & steprot,
		float3 * & linearvel, float3 * & angularvel)
{
	// Compute time step and time according to the integration scheme
	// TODO: must be done according to the integration scheme
	double dt1 = dt;
	if (step == 1)
		dt1 /= 2.0;
	double t0 = t;
	double t1 = t + dt1;

	//#define _DEBUG_OBJ_FORCES_
	bool ode_bodies = false;
	// For ODE bodies apply forces and torques
	for (int i = 0; i < m_bodies.size(); i++) {
		// Shortcut to body data
		MovingBodyData* mbdata = m_bodies[i];
		// Store kinematic data at the beginning of the time step
		if (step == 1)
			m_bodies_storage[i] = mbdata->kdata;
		// Restore kinematic data from the value stored at the beginning of the time step
		if (step == 2)
			mbdata->kdata = m_bodies_storage[i];

		if (mbdata->type == MB_ODE) {
			ode_bodies = true;
			const dBodyID bodyid = mbdata->object->m_ODEBody;
			// For step 2 restore cg, lvel and avel to the value at the beginning of
			// the timestep
			if (step == 2) {
				dBodySetPosition(bodyid, (dReal) mbdata->kdata.crot.x, (dReal) mbdata->kdata.crot.y,
								(dReal) mbdata->kdata.crot.z);
				dBodySetLinearVel(bodyid, (dReal) mbdata->kdata.lvel.x, (dReal) mbdata->kdata.lvel.y,
								(dReal) mbdata->kdata.lvel.z);
				dBodySetAngularVel(bodyid, (dReal) mbdata->kdata.avel.x, (dReal) mbdata->kdata.avel.y,
								(dReal) mbdata->kdata.avel.z);
				dQuaternion quat;
				mbdata->kdata.orientation.ToODEQuaternion(quat);
				dBodySetQuaternion(bodyid, quat);
			}
			dBodyAddForce(bodyid, forces[i].x, forces[i].y, forces[i].z);
			dBodyAddTorque(bodyid, torques[i].x, torques[i].y, torques[i].z);

			#ifdef _DEBUG_OBJ_FORCES_
			cout << "Before dWorldStep, object " << i << "\tt = " << t << "\tdt = " << dt <<"\n";
			//mbdata->object->ODEPrintInformation(false);
			printf("   F:	%e\t%e\t%e\n", forces[i].x, forces[i].y, forces[i].z);
			printf("   T:	%e\t%e\t%e\n", torques[i].x, torques[i].y, torques[i].z);
			#endif
		}
	}

	// Call ODE solver for ODE bodies
	if (ode_bodies) {
		dSpaceCollide(m_ODESpace, (void *) this, &ODE_near_callback_wrapper);
		dWorldStep(m_ODEWorld, dt1);
		if (m_ODEJointGroup)
			dJointGroupEmpty(m_ODEJointGroup);
	}

	// Walk trough all moving bodies :
	// updates bodies center of rotation, linear and angular velocity and orientation
	for (int i = 0; i < m_bodies.size(); i++) {
		// Shortcut to MovingBodyData
		MovingBodyData* mbdata = m_bodies[i];
		// New center of rotation, linear and angular velocity and orientation
		double3 new_trans = make_double3(0.0);
		EulerParameters new_orientation, dr;
		// In case of an ODE body, new center of rotation position, linear and angular velocity
		// and new orientation have been computed by ODE
		if (mbdata->type == MB_ODE) {
			const dBodyID bodyid = mbdata->object->m_ODEBody;
			const double3 new_crot = make_double3(dBodyGetPosition(bodyid));
			new_trans = new_crot - mbdata->kdata.crot;
			mbdata->kdata.crot = new_crot;
			mbdata->kdata.lvel = make_double3(dBodyGetLinearVel(bodyid));
			mbdata->kdata.avel = make_double3(dBodyGetAngularVel(bodyid));
			const EulerParameters new_orientation = EulerParameters(dBodyGetQuaternion(bodyid));
			dr = new_orientation*mbdata->kdata.orientation.Inverse();
			mbdata->kdata.orientation = new_orientation;
		}
		// Otherwise the user is providing linear and angular velocity trough a call back
		// function
		else {
			const uint index = mbdata->index;
			// Get linear and angular velocities at t + dt/2.O for step 1 or t + dt for step 2
			float3 force = make_float3(0.0f);
			float3 torque = make_float3(0.0f);
			if (mbdata->type == MB_FORCES_MOVING) {
				force = forces[i];
				torque = torques[i];
			}

			moving_bodies_callback(index, mbdata->object, t0, t1, force, torque, mbdata->initial_kdata,
					mbdata->kdata, new_trans, dr);
		}

		calc_grid_and_local_pos(mbdata->kdata.crot, cgGridPos + i, cgPos + i);
		trans[i] = make_float3(new_trans);
		linearvel[i] = make_float3(mbdata->kdata.lvel);
		angularvel[i] = make_float3(mbdata->kdata.avel);

		// Compute and relative rotation respect to the beginning of time step
		float *base_addr = steprot + 9*i;
		dr.ComputeRot();
		dr.GetRotation(base_addr);

		#ifdef _DEBUG_OBJ_FORCES_
		if (i == 1 && trans[i].x != 0.0) {
		cout << "After dWorldStep, object "  << i << "\tt = " << t << "\tdt = " << dt <<"\n";
		mbdata->object->ODEPrintInformation(false);
		printf("   lvel: %e\t%e\t%e\n", linearvel[i].x, linearvel[i].y, linearvel[i].z);
		printf("   avel: %e\t%e\t%e\n", angularvel[i].x, angularvel[i].y, angularvel[i].z);
		printf("    pos: %g\t%g\t%g\n", mbdata->kdata.crot.x, mbdata->kdata.crot.y, mbdata->kdata.crot.z);
		printf("   gpos: %d\t%d\t%d\n", cgGridPos[i].x, cgGridPos[i].y, cgGridPos[i].z);
		printf("   lpos: %e\t%e\t%e\n", cgPos[i].x, cgPos[i].y, cgPos[i].z);
		printf("   trans:%e\t%e\t%e\n", trans[i].x, trans[i].y, trans[i].z);
		printf("   n_ep: %e\t%e\t%e\t%e\n", mbdata->kdata.orientation(0), mbdata->kdata.orientation(1),
				mbdata->kdata.orientation(2), mbdata->kdata.orientation(3));
		printf("   dr: %e\t%e\t%e\t%e\n", dr(0), dr(1),dr(2), dr(3));
		printf("   SR:   %e\t%e\t%e\n", base_addr[0], base_addr[1], base_addr[2]);
		printf("         %e\t%e\t%e\n", base_addr[3], base_addr[4], base_addr[5]);
		printf("         %e\t%e\t%e\n\n", base_addr[6], base_addr[7], base_addr[8]);
		}
		#endif
	}
}

// Copy planes for upload
void
Problem::copy_planes(PlaneList& planes)
{
	return;
}

void
Problem::check_dt(void)
{
	float dt_from_sspeed = INFINITY;
	for (uint f = 0 ; f < physparams()->numFluids(); ++f) {
		float sspeed = physparams()->sscoeff[f];
		dt_from_sspeed = fmin(dt_from_sspeed, simparams()->slength/sspeed);
	}
	dt_from_sspeed *= simparams()->dtadaptfactor;

	float dt_from_gravity = sqrt(simparams()->slength/length(physparams()->gravity));
	dt_from_gravity *= simparams()->dtadaptfactor;

	float dt_from_visc = NAN;
	if (simparams()->visctype != ARTVISC) {
		for (uint f = 0; f < physparams()->numFluids(); ++f)
			dt_from_visc = fminf(dt_from_visc, simparams()->slength*simparams()->slength/physparams()->kinematicvisc[f]);
		dt_from_visc *= 0.125f; // TODO this should be configurable
	}

	float cfl_dt = fminf(dt_from_sspeed, fminf(dt_from_gravity, dt_from_visc));

	if (simparams()->dt > cfl_dt) {
		fprintf(stderr, "WARNING: dt %g bigger than %g imposed by CFL conditions (sspeed: %g, gravity: %g, viscosity: %g)\n",
			simparams()->dt, cfl_dt,
			dt_from_sspeed, dt_from_gravity, dt_from_visc);
	} else if (!simparams()->dt) { // dt wasn't set
			simparams()->dt = cfl_dt;
			printf("setting dt = %g from CFL conditions (soundspeed: %g, gravity: %g, viscosity: %g)\n",
				simparams()->dt,
				dt_from_sspeed, dt_from_gravity, dt_from_visc);
	} else {
			printf("dt = %g (CFL conditions from soundspeed: %g, from gravity %g, from viscosity %g)\n",
				simparams()->dt,
				dt_from_sspeed, dt_from_gravity, dt_from_visc);
	}

}

void
Problem::check_maxneibsnum(void)
{
	// kernel radius times smoothing factor, rounded to the next integer
	double r = simparams()->sfactor*simparams()->kernelradius;
	r = ceil(r);

	// volumes are computed using a coefficient which is sligthly more than π
#define PI_PLUS_EPS 3.2
	double vol = 4*PI_PLUS_EPS*r*r*r/3;
	// and rounded up
	vol = ceil(vol);

	// maxneibsnum is obtained rounding up the volume to the next
	// multiple of 32
	uint maxneibsnum = round_up((uint)vol, 32U);

	// with semi-analytical boundaries, boundary particles
	// are doubled, so we expand by a factor of 1.5,
	// again rounding up
	if (simparams()->boundarytype == SA_BOUNDARY)
		maxneibsnum = round_up(3*maxneibsnum/2, 32U);

	// more in general, it's possible to have different particle densities for the
	// boundaries even with other boundary conditions. we do not have a universal
	// parameter that marks the inter-particle distance for boundary particles,
	// although we know that r0 is normally used for this too.
	// TODO FIXME when the double meaning of r0 as inter-particle distance for
	// boundary particles and as fluid-boundary distance is split into separate
	// variables, the inter-particle distance should be used in the next formula

	// The formula we use is based on the following:
	// 1. a half-sphere has (3/2) pi r^3 particle
	// 2. a full circle has pi (r/q)^2 particles, if q is the ratio beween
	//   the inter-particle distance on the full circle and the inter-particle
	//   distance used in the fluid
	// * the number of neighbors that are seen by a particle which is near
	//   a boundary plane with q*dp interparticle-distance is augmented the number
	//   in 2. over the number in 1., giving (3/2) (1/q)^2 (1/r)
	// * of course this does not affect the entire neighborhood, but only the part
	//   which is close to a boundary, which we estimate to be at most 2/3rds of
	//   the neighborhood, which cancels with the (3/2) factor
	//   TODO check if we should assume 7/8ths instead (particle near vertex
	//   only has 1/8th of a sphere in the fluid, the rest is all boundaries).
	double qq = m_deltap/physparams()->r0; // 1/q
	// double ratio = fmax((21*qq*qq)/(16*r), 1.0); // if we assume 7/8
	double ratio = fmax((qq*qq)/r, 1.0); // only use this if it gives us _more_ particles
	// increase maxneibsnum as appropriate
	maxneibsnum = (uint)ceil(ratio*maxneibsnum);
	// round up to multiple of 32
	maxneibsnum = round_up(maxneibsnum, 32U);

	// if the maxneibsnum was user-set, check against computed minimum
	if (simparams()->maxneibsnum) {
		if (simparams()->maxneibsnum < maxneibsnum) {
			fprintf(stderr, "WARNING: problem-set max neibs num too low! %u < %u\n",
				simparams()->maxneibsnum, maxneibsnum);
		} else {
			printf("Using problem-set max neibs num %u (safe computed value was %u)\n",
				simparams()->maxneibsnum, maxneibsnum);
		}
	} else {
		printf("Using computed max neibs num %u\n", maxneibsnum);
		simparams()->maxneibsnum = maxneibsnum;
	}
}


float
Problem::density(float h, int i) const
{
	float density = physparams()->rho0[i];

	if (h > 0) {
		//float g = length(physparams()->gravity);
		float g = abs(physparams()->gravity.z);
		// TODO g*rho0*h/B could be simplified to g*h*gamma/(c0*c0)
		density = physparams()->rho0[i]*pow(g*physparams()->rho0[i]*h/physparams()->bcoeff[i] + 1,
				1/physparams()->gammacoeff[i]);
		}
	return density;
}

// density to achieve a specific pressure
float
Problem::density_for_pressure(float P, int i) const
{
	return  physparams()->rho0[i]*pow(P/physparams()->bcoeff[i] + 1,
				1/physparams()->gammacoeff[i]);
}


float
Problem::soundspeed(float rho, int i) const
{
	return physparams()->sscoeff[i]*pow(rho/physparams()->rho0[i], physparams()->sspowercoeff[i]);
}


float
Problem::pressure(float rho, int i) const
{
	return physparams()->bcoeff[i]*(pow(rho/physparams()->rho0[i], physparams()->gammacoeff[i]) - 1);
}

void
Problem::add_gage(double3 const& pt)
{
	simparams()->gage.push_back(make_double4(pt.x, pt.y, 0., pt.z));
}

plane_t
Problem::implicit_plane(double4 const& p)
{
	const double4 midPoint = make_double4(m_origin + m_size/2, 1.0);

	plane_t plane;
	const double norm = length3(p);
	const double3 normal = as_double3(p)/norm;
	plane.normal = make_float3(normal);

	/* For the plane point, we pick the one closest to the center of the domain
	 * TODO find a better logic ? */

	const double midDist = dot(midPoint, p)/norm;
	double3 planePoint = as_double3(midPoint) - midDist*normal;

	calc_grid_and_local_pos(planePoint, &plane.gridPos, &plane.pos);

	return plane;
}

plane_t
Problem::make_plane(Point const& pt, Vector const& normal)
{
	plane_t plane;

	plane.normal = make_float3(normal);
	calc_grid_and_local_pos(make_double3(pt), &plane.gridPos, &plane.pos);

	return plane;
}

std::string const&
Problem::create_problem_dir(void)
{
	// if no data save directory was specified, default to a name
	// composed of problem name followed by date and time
	if (m_problem_dir.empty()) {
		time_t  rawtime;
		char	time_str[18];

		time(&rawtime);
		strftime(time_str, 18, "_%Y-%m-%dT%Hh%M", localtime(&rawtime));
		time_str[17] = '\0';
		// if "./tests/" doesn't exist yet...
		mkdir("./tests/", S_IRWXU | S_IRWXG | S_IRWXO);
		m_problem_dir = "./tests/" + m_name + std::string(time_str);
	}

	// TODO it should be possible to specify a directory with %-like
	// replaceable strings, such as %{problem} => problem name,
	// %{time} => launch time, etc.

	mkdir(m_problem_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

	return m_problem_dir;
}

// timer tick, for compatibility with old timer-tick writer frequency API
// remove when the old API is obsoleted
static double deprecated_timer_tick;

void
Problem::set_timer_tick(double t)
{
	fputs("WARNING: set_timer_tick() is deprecated\n", stderr);
	fputs("\tPlease use the floating-point version of add_writer() instead\n", stderr);
	deprecated_timer_tick = t;
}

void
Problem::add_writer(WriterType wt, int freq)
{
	fputs("WARNING: add_writer(WriterType, int) is deprecated\n", stderr);
	fputs("\tPlease use the floating-point version of add_writer() instead\n", stderr);
	add_writer(wt, freq*deprecated_timer_tick);
}

void
Problem::add_writer(WriterType wt, double freq)
{
	m_writers.push_back(make_pair(wt, freq));
}


// override in problems where you want to save
// at specific times regardless of standard conditions
bool
Problem::need_write(double t) const
{
	return false;
}

// overridden in subclasses if they want to write custom stuff
// using the CALLBACKWRITER
void
Problem::writer_callback(CallbackWriter *,
	uint numParts, BufferList const&, uint node_offset, double t,
	const bool testpoints) const
{
	fprintf(stderr, "WARNING: CallbackWriter is being used, but writer_callback wasn't implemented\n");
}


// is the simulation finished at the given time?
bool
Problem::finished(double t) const
{
	double tend(simparams()->tend);
	return tend && (t > tend);
}


float3
Problem::g_callback(const double t)
{
	/* If this was not overridden, it's likely that the caller overridden the deprecated
	 * float version, passthrough */
	static bool reminder_shown = false;
	if (!reminder_shown) {
		fprintf(stderr, "WARNING: g_callback(float) is deprecated, please switch to g_callback(double)\n");
		reminder_shown = true;
	}
	IGNORE_WARNINGS(deprecated-declarations)
	return g_callback(float(t));
	RESTORE_WARNINGS
}

float3
Problem::g_callback(const float t)
{
	static bool reminder_shown = false;
	if (!reminder_shown) {
		fprintf(stderr, "WARNING: gravity callback enabled but not overridden\n");
		reminder_shown = true;
	}
	return make_float3(0.0);
}



// Fill the device map with "devnums" (*global* device ids) in range [0..numDevices[.
// Default algorithm: split along the longest axis
void Problem::fillDeviceMap()
{
	fillDeviceMapByAxis(LONGEST_AXIS);
}

// partition by splitting the cells according to their linearized hash.
void Problem::fillDeviceMapByCellHash()
{
	uint cells_per_device = gdata->nGridCells / gdata->totDevices;
	for (uint i=0; i < gdata->nGridCells; i++)
		// guaranteed to fit in a devcount_t due to how it's computed
		gdata->s_hDeviceMap[i] = devcount_t(min( i/cells_per_device, gdata->totDevices-1));
}

// partition by splitting along the specified axis
void Problem::fillDeviceMapByAxis(SplitAxis preferred_split_axis)
{
	// select the longest axis
	if (preferred_split_axis == LONGEST_AXIS) {
		if (	gdata->worldSize.x >= gdata->worldSize.y &&
				gdata->worldSize.x >= gdata->worldSize.z)
			preferred_split_axis = X_AXIS;
		else
		if (	gdata->worldSize.y >= gdata->worldSize.z)
			preferred_split_axis = Y_AXIS;
		else
			preferred_split_axis = Z_AXIS;
	}
	uint cells_per_split_axis = 0;
	switch (preferred_split_axis) {
		case X_AXIS:
			cells_per_split_axis = gdata->gridSize.x;
			break;
		case Y_AXIS:
			cells_per_split_axis = gdata->gridSize.y;
			break;
		case Z_AXIS:
			cells_per_split_axis = gdata->gridSize.z;
			break;
	}

	// Check that we have enough cells along the split axis. This check should
	// be performed in all split algorithms
	if (cells_per_split_axis / (double) gdata->totDevices < 3.0)
		throw runtime_error ("FATAL: not enough cells along the split axis. Aborting.\n");

	uint cells_per_device_per_split_axis = (uint)round(cells_per_split_axis / (double)gdata->totDevices);

	/*
	printf("Splitting domain along axis %s, %u cells per part\n",
		(preferred_split_axis == X_AXIS ? "X" : (preferred_split_axis == Y_AXIS ? "Y" : "Z") ), cells_per_device_per_split_axis);
	*/
	for (uint cx = 0; cx < gdata->gridSize.x; cx++)
		for (uint cy = 0; cy < gdata->gridSize.y; cy++)
			for (uint cz = 0; cz < gdata->gridSize.z; cz++) {
				uint axis_coordinate;
				switch (preferred_split_axis) {
					case X_AXIS: axis_coordinate = cx; break;
					case Y_AXIS: axis_coordinate = cy; break;
					case Z_AXIS: axis_coordinate = cz; break;
				}
				// everything is just a preparation for the following line
				devcount_t dstDevice = devcount_t(axis_coordinate / cells_per_device_per_split_axis);
				// handle the case when cells_per_split_axis multiplies cells_per_split_axis
				dstDevice = (devcount_t)min(dstDevice, gdata->totDevices - 1);
				// compute cell address
				uint cellLinearHash = gdata->calcGridHashHost(cx, cy, cz);
				// assign it
				gdata->s_hDeviceMap[cellLinearHash] = dstDevice;
			}
}

// Like fillDeviceMapByAxis(), but splits are proportional to the contained fluid particles
void Problem::fillDeviceMapByAxisBalanced(SplitAxis preferred_split_axis)
{
	// Select the longest axis
	if (preferred_split_axis == LONGEST_AXIS) {
		if (	gdata->worldSize.x >= gdata->worldSize.y &&
				gdata->worldSize.x >= gdata->worldSize.z)
			preferred_split_axis = X_AXIS;
		else
		if (	gdata->worldSize.y >= gdata->worldSize.z)
			preferred_split_axis = Y_AXIS;
		else
			preferred_split_axis = Z_AXIS;
	}

	// Set some aux variables - axis 1 is the split axis
	uint cells_per_axis1 = 0;
	uint cells_per_axis2 = 0;
	uint cells_per_axis3 = 0;
	uint cx = 0, cy = 0, cz = 0; // cell coordinates
	uint *c1, *c2, *c3; // abstract from cell coordinates
	uint *axisParticleCounter = NULL;
	switch (preferred_split_axis) {
		case X_AXIS:
			cells_per_axis1 = gdata->gridSize.x;
			cells_per_axis2 = gdata->gridSize.y;
			cells_per_axis3 = gdata->gridSize.z;
			c1 = &cx;
			c2 = &cy;
			c3 = &cz;
			axisParticleCounter = gdata->s_hPartsPerSliceAlongX;
			break;
		case Y_AXIS:
			cells_per_axis1 = gdata->gridSize.y;
			cells_per_axis2 = gdata->gridSize.x;
			cells_per_axis3 = gdata->gridSize.z;
			c1 = &cy;
			c2 = &cx;
			c3 = &cz;
			axisParticleCounter = gdata->s_hPartsPerSliceAlongY;
			break;
		case Z_AXIS:
			cells_per_axis1 = gdata->gridSize.z;
			cells_per_axis2 = gdata->gridSize.x;
			cells_per_axis3 = gdata->gridSize.y;
			c1 = &cz;
			c2 = &cx;
			c3 = &cy;
			axisParticleCounter = gdata->s_hPartsPerSliceAlongZ;
			break;
	}

	// Check that we have enough cells along the split axis. This check should
	// be performed in all split algorithms
	if (cells_per_axis1 / (double) gdata->totDevices < 3.0)
		throw runtime_error ("FATAL: not enough cells along the split axis. Aborting.\n");

	// Compute ideal split values
	const uint particles_per_device = gdata->totParticles / gdata->totDevices;
	const uint particles_per_slice = gdata->totParticles / cells_per_axis1;

	// If a device has "almost" particles_per_device particles, next slice will assign particles_per_slice more
	// particles and make it "overflow" the ideal number; we will instead stop before, at this threshold
	const uint particles_per_device_threshold = particles_per_device - (particles_per_slice / 2);

	// printf("Splitting domain along axis %s, ~%u particles per device\n",
	//	(preferred_split_axis == X_AXIS ? "X" : (preferred_split_axis == Y_AXIS ? "Y" : "Z") ), (uint)particles_per_device);

	// We need at least 3 cells per device, regardless the distribution of fluid; so we track the
	// remaining cells which need to be "reserved" for device numbers yet to be analyzed, excluding
	// the first.
	uint reserved_cells =  3 * (gdata->totDevices - 1);

	// We will iterate on the cells and increase the current device number
	uint currentDevice = 0;
	uint currentDeviceParticles = 0;

	// NOTE: not using "*cx++" since post-increment has precedence over deference
	for (*c1 = 0; *c1 < cells_per_axis1; (*c1)++) {
		// We must increase the current device only if:
		// 1. This is not the last device (i.e. should always be currentDevice < totDevices), and
		// 2a. we got enough particles in previous iteration, or
		// 2b. we reached the reserved cells (thus, we must leave them for next devices)
		if ( (currentDevice < gdata->totDevices - 1) &&
			(currentDeviceParticles >= particles_per_device_threshold ||
			*c1 >= cells_per_axis1 - reserved_cells - 1) ) {
			// switch to next device: reset counter,
			currentDeviceParticles = 0;
			// increase device,
			currentDevice++;
			// update reserved_cells (minus mine)
			reserved_cells -= 3;
		}

		// add particles in current slice
		currentDeviceParticles += axisParticleCounter[ *c1 ];

		// assign all the cells of the current slice to current device
		for (*c2 = 0; *c2 < cells_per_axis2; (*c2)++)
			for (*c3 = 0; *c3 < cells_per_axis3; (*c3)++) {
				// we are actually using c1, c2, c3 in a proper order
				const uint cellLinearHash = gdata->calcGridHashHost(cx, cy, cz);
				// assign it
				gdata->s_hDeviceMap[cellLinearHash] = currentDevice;
			}
	} // iterate on split axis
}

void Problem::fillDeviceMapByEquation()
{
	// 1st equation: diagonal plane. (x+y+z)=coeff
	//uint longest_grid_size = max ( max( gdata->gridSize.x, gdata->gridSize.y), gdata->gridSize.z );
	uint coeff = (gdata->gridSize.x + gdata->gridSize.y + gdata->gridSize.z) / gdata->totDevices;
	// 2nd equation: sphere. Sqrt(cx²+cy²+cz²)=radius
	uint diagonal = (uint) sqrt(	gdata->gridSize.x * gdata->gridSize.x +
									gdata->gridSize.y * gdata->gridSize.y +
									gdata->gridSize.z * gdata->gridSize.z) / 2;
	uint radius_part = diagonal /  gdata->totDevices;
	for (uint cx = 0; cx < gdata->gridSize.x; cx++)
		for (uint cy = 0; cy < gdata->gridSize.y; cy++)
			for (uint cz = 0; cz < gdata->gridSize.z; cz++) {
				uint dstDevice;
				// 1st equation: rough oblique plane split --
				dstDevice = (cx + cy + cz) / coeff;
				// -- end of 1st eq.
				// 2nd equation: spheres --
				//uint distance_from_origin = (uint) sqrt( cx * cx + cy * cy + cz * cz);
				// comparing directly the square would be more efficient but could require long uints
				//dstDevice = distance_from_origin / radius_part;
				// -- end of 2nd eq.
				// handle special cases at the edge
				dstDevice = min(dstDevice, gdata->totDevices - 1);
				// compute cell address
				uint cellLinearHash = gdata->calcGridHashHost(cx, cy, cz);
				// assign it
				gdata->s_hDeviceMap[cellLinearHash] = (uchar)dstDevice;
			}
}

// Partition by performing the splitting the domain in the specified number of slices for each axis.
// Values must be > 0. The number of devices will be the product of the input values.
// This is not meant to be called directly by a problem since the number of splits (and thus the devices)
// would be hardocded. A wrapper method (like fillDeviceMapByRegularGrid) can provide an algorithm to
// properly factorize a given number of GPUs in 2 or 3 values.
void Problem::fillDeviceMapByAxesSplits(uint Xslices, uint Yslices, uint Zslices)
{
	// is any of these zero?
	if (Xslices * Yslices * Zslices == 0)
		printf("WARNING: fillDeviceMapByAxesSplits() called with zero values, using 1 instead");

	if (Xslices == 0) Xslices = 1;
	if (Yslices == 0) Yslices = 1;
	if (Zslices == 0) Zslices = 1;

	// divide and round
	uint devSizeCellsX = (gdata->gridSize.x + Xslices - 1) / Xslices ;
	uint devSizeCellsY = (gdata->gridSize.y + Yslices - 1) / Yslices ;
	uint devSizeCellsZ = (gdata->gridSize.z + Zslices - 1) / Zslices ;

	// iterate on all cells
	for (uint cx = 0; cx < gdata->gridSize.x; cx++)
			for (uint cy = 0; cy < gdata->gridSize.y; cy++)
				for (uint cz = 0; cz < gdata->gridSize.z; cz++) {

				// where are we in the 3D grid of devices?
				uint whichDevCoordX = (cx / devSizeCellsX);
				uint whichDevCoordY = (cy / devSizeCellsY);
				uint whichDevCoordZ = (cz / devSizeCellsZ);

				// round if needed
				if (whichDevCoordX == Xslices) whichDevCoordX--;
				if (whichDevCoordY == Yslices) whichDevCoordY--;
				if (whichDevCoordZ == Zslices) whichDevCoordZ--;

				// compute dest device
				uint dstDevice = whichDevCoordZ * Yslices * Xslices + whichDevCoordY * Xslices + whichDevCoordX;
				// compute cell address
				uint cellLinearHash = gdata->calcGridHashHost(cx, cy, cz);
				// assign it
				gdata->s_hDeviceMap[cellLinearHash] = (uchar)dstDevice;
			}
}

// Wrapper for fillDeviceMapByAxesSplits() computing the number of cuts along each axis.
// WARNING: assumes the total number of devices is divided by a combination of 2, 3 and 5
void Problem::fillDeviceMapByRegularGrid()
{
	float Xsize = gdata->worldSize.x;
	float Ysize = gdata->worldSize.y;
	float Zsize = gdata->worldSize.z;
	devcount_t cutsX = 1;
	devcount_t cutsY = 1;
	devcount_t cutsZ = 1;
	devcount_t remaining_factors = gdata->totDevices;

	// define the product of non-zero cuts to keep track of current number of parallelepipeds
//#define NZ_PRODUCT	((cutsX > 0? cutsX : 1) * (cutsY > 0? cutsY : 1) * (cutsZ > 0? cutsZ : 1))

	while (cutsX * cutsY * cutsZ < gdata->totDevices) {
		devcount_t factor = 1;
		// choose the highest factor among 2, 3 and 5 which divides remaining_factors
		if (remaining_factors % 5 == 0) factor = 5; else
		if (remaining_factors % 3 == 0) factor = 3; else
		if (remaining_factors % 2 == 0) factor = 2; else {
			factor = remaining_factors;
			printf("WARNING: splitting by regular grid but %u is not divided by 2,3,5!\n", remaining_factors);
		}
		// choose the longest axis to split along
		if (Xsize >= Ysize && Xsize >= Zsize) {
			Xsize /= factor;
			cutsX *= factor;
		} else
		if (Ysize >= Xsize && Ysize >= Zsize) {
			Ysize /= factor;
			cutsY *= factor;
		} else {
			Zsize /= factor;
			cutsZ *= factor;
		}
	}

	// should always hold, but double check for bugs
	if (cutsX * cutsY * cutsZ != gdata->totDevices)
		printf("WARNING: splitting by regular grid but final distribution (%u, %u, %u) does not produce %u parallelepipeds!\n",
			cutsX, cutsY, cutsZ, gdata->totDevices);

	fillDeviceMapByAxesSplits(cutsX, cutsY, cutsZ);
}



uint
Problem::max_parts(uint numParts)
{
	if (!(simparams()->simflags & ENABLE_INLET_OUTLET))
		return numParts;

	// we assume that we can't have more particles than by filling the whole domain:
	// if the user knows how many particles there are going to be he should implement
	// his own version of this function
	double3 range = get_worldsize();
	range /= m_deltap; // regular fill
	uint wparts = max(range.x,1)*max(range.y,1)*max(range.z,1);
	printf("  estimating %u particles to fill the world\n", wparts);

	return wparts;
}

// This function computes the Ferrari coefficient based on a length-scale. The formula for the coefficient
// is L/(1000 * deltap), see Mayrhofer et al., 2013. If the length scale is not set then the ferrari coefficient
// will be taken as it is, regardless of whether it is set or not (default value = 0)
void
Problem::calculateFerrariCoefficient()
{
	if (isnan(simparams()->ferrari)) {
		if (isnan(simparams()->ferrariLengthScale)) {
			simparams()->ferrari = 0.0f;
			printf("Ferrari coefficient: %e (default value, disabled)\n", simparams()->ferrari);
			if (simparams()->simflags & ENABLE_FERRARI)
				fprintf(stderr, "WARNING: Ferrari correction enabled, but no coefficient or length scale given!\n");
			return;
		}
		else {
			simparams()->ferrari = simparams()->ferrariLengthScale*1e-3f/m_deltap;
			printf("Ferrari coefficient: %e (computed from length scale: %e)\n", simparams()->ferrari, simparams()->ferrariLengthScale);
			return;
		}
	}
	printf("Ferrari coefficient: %e\n", simparams()->ferrari);
	return;
}


/*! Compute grid and cell size from the kernel influence radius
 * The number of cell is obtained as the ratio between the domain size and the
 * influence radius, rounded down to the closest integer.
 * The reason for rounding down is that we want the cell size to be no smaller
 * than the influence radius, to guarantee that all neighbors of a particle are
 * found at most one cell away in each direction.
 */
void
Problem::set_grid_params(void)
{
	/* When using periodicity, it's important that the world size in the periodic
	 * direction is an exact multiple of the deltap: if this is not the case,
	 * fluid filling might use an effective inter-particle distance which is
	 * “significantly” different from deltap, which would lead particles near
	 * periodic boundaries to have distance _exactly_ deltap across the boundary,
	 * but “significantly” different on the same side. While this in general would not
	 * be extremely important, it can have a noticeable effect at the beginning of the
	 * simulation, when particles are distributed quite regularly and the difference
	 * between effective (inner) distance and cross-particle distance can create
	 * a rather annoying discontinuity.
	 * So warn if m_size.{x,y,z} is not a multiple of deltap in case of periodicity.
	 * TODO FIXME this would not be needed if filling was made taking into account
	 * periodicity and spaced particles accordingly.
	 */
	if (simparams()->periodicbound & PERIODIC_X && !is_multiple(m_size.x, m_deltap))
		fprintf(stderr, "WARNING: problem is periodic in X, but X world size %.9g is not a multiple of deltap (%.g)\n",
			m_size.x, m_deltap);
	if (simparams()->periodicbound & PERIODIC_Y && !is_multiple(m_size.y, m_deltap))
		fprintf(stderr, "WARNING: problem is periodic in Y, but Y world size %.9g is not a multiple of deltap (%.g)\n",
			m_size.y, m_deltap);
	if (simparams()->periodicbound & PERIODIC_Z && !is_multiple(m_size.z, m_deltap))
		fprintf(stderr, "WARNING: problem is periodic in X, but Z world size %.9g is not a multiple of deltap (%.g)\n",
			m_size.z, m_deltap);

	double influenceRadius = simparams()->kernelradius*simparams()->slength;
	// with semi-analytical boundaries, we want a cell size which is
	// deltap/2 + the usual influence radius
	double cellSide = influenceRadius;
	if (simparams()->boundarytype == SA_BOUNDARY)
		cellSide += m_deltap/2.0f;

	m_gridsize.x = (uint)floor(m_size.x / cellSide);
	m_gridsize.y = (uint)floor(m_size.y / cellSide);
	m_gridsize.z = (uint)floor(m_size.z / cellSide);

	// While trying to run a simulation at very low resolution, the user might
	// set a deltap so large that cellSide is bigger than m_size.{x,y,z}, resulting
	// in a corresponding gridsize of 0. Check for this case (by checking if any
	// of the gridsize components are zero) and throw.

	if (!m_gridsize.x || !m_gridsize.y || !m_gridsize.z) {
		stringstream ss;
		ss << "resolution " << simparams()->slength << " is too low! Resulting grid size would be "
			<< m_gridsize;
		throw runtime_error(ss.str());
	}

	m_cellsize.x = m_size.x / m_gridsize.x;
	m_cellsize.y = m_size.y / m_gridsize.y;
	m_cellsize.z = m_size.z / m_gridsize.z;

	/*
	printf("set_grid_params->t:\n");
	printf("Domain size\t: (%f, %f, %f)\n", m_size.x, m_size.y, m_size.z);
	*/
	printf("Influence radius / expected cell side\t: %g, %g\n", influenceRadius, cellSide);
	/*
	printf("Grid   size\t: (%d, %d, %d)\n", m_gridsize.x, m_gridsize.y, m_gridsize.z);
	printf("Cell   size\t: (%f, %f, %f)\n", m_cellsize.x, m_cellsize.y, m_cellsize.z);
	printf("       delta\t: (%.2f%%, %.2f%%, %.2f%%)\n",
		(m_cellsize.x - cellSide)*100/cellSide,
		(m_cellsize.y - cellSide)*100/cellSide,
		(m_cellsize.z - cellSide)*100/cellSide);
	*/
}


// Compute position in uniform grid (clamping to edges)
int3
Problem::calc_grid_pos(const Point& pos) const
{
	int3 gridPos;
	gridPos.x = (int)floor((pos(0) - m_origin.x) / m_cellsize.x);
	gridPos.y = (int)floor((pos(1) - m_origin.y) / m_cellsize.y);
	gridPos.z = (int)floor((pos(2) - m_origin.z) / m_cellsize.z);
	gridPos.x = min(max(0, gridPos.x), m_gridsize.x-1);
	gridPos.y = min(max(0, gridPos.y), m_gridsize.y-1);
	gridPos.z = min(max(0, gridPos.z), m_gridsize.z-1);

	return gridPos;
}

/// Compute the uniform grid components of a vector
int3
Problem::calc_grid_offset(double3 const& vec) const
{
	int3 gridOff;
	gridOff = make_int3(floor(vec/m_cellsize));

	return gridOff;
}

/// Compute the local (fractional grid cell) components of a vector,
/// given the vector and its grid offset
double3
Problem::calc_local_offset(double3 const& vec, int3 const& gridOff) const
{
	return vec - (make_double3(gridOff) + 0.5)*m_cellsize;
}


// Compute address in grid from position
uint
Problem::calc_grid_hash(int3 gridPos) const
{
	return gridPos.COORD3 * m_gridsize.COORD2 * m_gridsize.COORD1 + gridPos.COORD2 * m_gridsize.COORD1 + gridPos.COORD1;
}


void
Problem::calc_localpos_and_hash(const Point& pos, const particleinfo& info, float4& localpos, hashKey& hash) const
{
	int3 gridPos = calc_grid_pos(pos);

	// automatically choose between long hash (cellHash + particleId) and short hash (cellHash)
	hash = calc_grid_hash(gridPos);

	localpos.x = float(pos(0) - m_origin.x - (gridPos.x + 0.5)*m_cellsize.x);
	localpos.y = float(pos(1) - m_origin.y - (gridPos.y + 0.5)*m_cellsize.y);
	localpos.z = float(pos(2) - m_origin.z - (gridPos.z + 0.5)*m_cellsize.z);
	localpos.w = float(pos(3));
}

void
Problem::init_keps(float* k, float* e, uint numpart, particleinfo* info, float4* pos, hashKey* hash)
{
	const float Lm = fmax(2*m_deltap, 1e-5f);
	const float k0 = pow(0.002f*physparams()->sscoeff[0], 2);
	const float e0 = 0.16f*pow(k0, 1.5f)/Lm;

	for (uint i = 0; i < numpart; i++) {
		k[i] = k0;
		e[i] = e0;
	}
}

/* Initialize the particle volumes from their masses and densities. */
void
Problem::init_volume(BufferList &buffers, uint numParticles)
{
	float4 *pos = buffers.getData<BUFFER_POS>();
	float4 *vel = buffers.getData<BUFFER_VEL>();
	float4 *vol = buffers.getData<BUFFER_VOLUME>();

	for (uint i = 0; i < numParticles; ++i) {
		float4 pvol;
		// .x: initial volume, .w current volume.
		// at the beginning they are both equal to mass/density
		pvol.x = pvol.w = pos[i].w/vel[i].w;
		// .y is the log of current/initial
		pvol.y = 0;
		// .z is unused, set to zero
		pvol.z = 0;

		vol[i] = pvol;
	}
}

void
Problem::imposeBoundaryConditionHost(
			MultiBufferList::iterator		bufwrite,
			MultiBufferList::const_iterator	bufread,
					uint*			IOwaterdepth,
			const	float			t,
			const	uint			numParticles,
			const	uint			numOpenBoundaries,
			const	uint			particleRangeEnd)
{
	fprintf(stderr, "WARNING: open boundaries are present, but imposeBoundaryCondtionHost was not implemented\n");
	return;
}

void Problem::imposeForcedMovingObjects(
			float3	&gravityCenters,
			float3	&translations,
			float*	rotationMatrices,
	const	uint	ob,
	const	double	t,
	const	float	dt)
{
	// not implemented
	return;
}




void Problem::PlaneCut(PointVect& points, const double a, const double b,
			const double c, const double d)
{
	PointVect new_points;
	new_points.reserve(points.size());
	//const double norm_factor = sqrt(a*a + b*b + c*c);
	for (uint i = 0; i < points.size(); i++) {
		const Point & p = points[i];
		const double dist = a*p(0) + b*p(1) + c*p(2) + d;

		if (dist >= 0)
			new_points.push_back(p);
	}

	points.clear();
	points = new_points;
}
