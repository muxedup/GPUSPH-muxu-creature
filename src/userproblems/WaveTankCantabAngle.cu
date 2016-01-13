/*  Copyright 2011 Alexis Herault, Giuseppe Bilotta, Robert A. Dalrymple, Eugenio Rustico, Ciro Del Negro

	Istituto de Nazionale di Geofisica e Vulcanologia
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

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "WaveTank.h"
#include "particledefine.h"
#include "GlobalData.h"
#include "cudasimframework.cu"


#define MK_par 2

WaveTank::WaveTank(GlobalData *_gdata) : Problem(_gdata)
{
	// Size and origin of the simulation domain
	lx = 22.;
	ly = 32.0;
	lz = 3.5;
	
	// Data for problem setup
	slope_length = 10.2;  // horizontal extent of tank covered by sloping beach
	h_length = 10;
	height = 2.5;
	beta = .1974;
	
	
    
	SETUP_FRAMEWORK(
	    //viscosity<ARTVISC>,
		//viscosity<KINEMATICVISC>,
		viscosity<SPSVISC>,
		boundary<LJ_BOUNDARY>
		//boundary<MK_BOUNDARY>
	);
	
	m_size = make_double3(lx, ly, lz + 2.0*height);
	m_origin = make_double3(0.0, 0.0, -1.0*height);	 
	 
	addFilter(SHEPARD_FILTER, 20);
	  //MLS_FILTER

	//    Here we have npaddles paddles for a directional seastate
 
	// Add objects to the tank
    use_cyl = false;
	use_cone = false;

	// use a plane for the bottom
	use_bottom_plane = true; 

	// SPH parameters
	set_deltap(0.15);  //0.005f;
	printf("deltap = %f \n", m_deltap);
	m_simparams->dt = 0.0001;
	m_simparams->dtadaptfactor = 0.2;
	m_simparams->buildneibsfreq = 10;
	m_simparams->tend = 55.0; //seconds

//	m_simparams->vorticity = false;
	//Testpoints
//	m_simparams->testpoints = false;

	// Free surface detection
//	m_simparams->surfaceparticle = false;
//	m_simparams->savenormals = false;
/*
	//WaveGage
	add_gage(5., 16.);
	add_gage(14., 16.);
	add_gage(5.,8.); 
	add_gage(14.,8.);

*/
 

    // Physical parameters
	H = 1.0;
	m_physparams->gravity = make_float3(0.0, 0.0, -9.81);
	float g = length(m_physparams->gravity);

	float r0 = m_deltap;
	m_physparams->r0 = r0;
	
	add_fluid(1000.f);
	set_equation_of_state(0, 7.0f, 20.f);
	set_kinematic_visc(0,1.0e-6);
	
    m_physparams->artvisccoeff= 0.2f;
	m_physparams->smagfactor = 0.12*0.12*m_deltap*m_deltap;
	m_physparams->kspsfactor = (2.0/3.0)*0.0066*m_deltap*m_deltap;
	m_physparams->epsartvisc = 0.01*m_simparams->slength*m_simparams->slength;

	// BC when using LJ
	m_physparams->dcoeff = 5.0*g*H;
    //set p1coeff,p2coeff, epsxsph here if different from 12.,6., 0.5

	// BC when using MK
	m_physparams->MK_K = g*H;
	m_physparams->MK_d = 1.1*m_deltap/MK_par;
	m_physparams->MK_beta = MK_par;

	//Wave paddle definition:  location, start & stop times, stroke and frequency (2 \pi/period)
	float waveangle=(3.1415927/180.)*8.0;
	// lateral wavenumber k sin(theta)
	float k=0.58;  //fix for every case
	lambda=k*sin(waveangle);
    npaddles= 16;
	paddle_width= ly/npaddles;	 
	paddle_length = 1.5; //vertical 
	paddle_origin = make_double3(1.0,0.0, 0.);
	paddle_tstart = 0.2;
	paddle_tend = m_simparams->tend;
	// The stroke value is given at free surface level H
	float stroke = .75;
	// m_mbamplitude is the maximal angular value par paddle angle
	// Paddle angle is in [-m_mbamplitude, m_mbamplitude]
	paddle_amplitude = atan(stroke/(2.0*(H - paddle_origin.z)));
	paddle_omega = 2.0*M_PI/2.6;		// period T = 3.6145 s
	
		
		//add wave gages to measure feedback in front of each paddle  
	 //   add_gage(mbpaddledata.amplitude+5*r0,(i+0.5)*paddle_width);
  
	// Drawing and saving times
	add_writer(VTKWRITER,.1);
	
	// Name of problem used for directory creation
	m_name = "WaveTank";
}


WaveTank::~WaveTank(void)
{
	release_memory();
}


void WaveTank::release_memory(void)
{
	parts.clear();
	boundary_parts.clear();
//	test_points.clear();
}

void
WaveTank::moving_bodies_callback(const uint index, Object* object, const double t0, const double t1,
			const float3& force, const float3& torque, const KinematicData& initial_kdata,
			KinematicData& kdata, double3& dx, EulerParameters& dr)
{
    // index gives the number of the paddles
    dx= make_double3(0.0);
    kdata.lvel=make_double3(0.0f, 0.0f, 0.0f);
    lambda_y = lambda * paddle_width *(2*index - 1)/2.;
    if (t1> paddle_tstart & t1 < paddle_tend){
       kdata.avel = make_double3(0.0, paddle_amplitude*paddle_omega*sin(lambda_y+paddle_omega*(t1-paddle_tstart)),0.0);
       EulerParameters dqdt = 0.5*EulerParameters(kdata.avel)*kdata.orientation;
       dr = EulerParameters::Identity() + (t1-t0)*dqdt*kdata.orientation.Inverse();
       dr.Normalize();
	   kdata.orientation = kdata.orientation + (t1 - t0)*dqdt;
	   kdata.orientation.Normalize();
	   }
	else {
	   kdata.avel = make_double3(0.0,0.0,0.0);
	   kdata.orientation = kdata.orientation;
	   dr.Identity();
	}
	
}


int WaveTank::fill_parts()
{
    std::cout << "filling particles" <<"\n";
	const float r0 = m_physparams->r0;
	const float br = (m_simparams->boundarytype == MK_BOUNDARY ? m_deltap/MK_par : r0);
	
	boundary_parts.reserve(100);
	paddle_parts.reserve(500);
	parts.reserve(34000);
   std::cout <<"npaddles = " << npaddles <<"\n";
    experiment_box = Cube(Point(0, 0, 0), h_length + slope_length,ly, height);
//	std::cout << "paddle origins: \n" <<"\n";
    // define the paddles
      const float amplitude = -paddle_amplitude ;
    for (uint i=0; i<npaddles; i++) {	
        double y_paddle = paddle_width*i;
        lambda_y = lambda * (y_paddle+ paddle_width/2.); 
		paddle[i] = Rect(Point(paddle_origin.x, paddle_origin.y+y_paddle, paddle_origin.z), Vector(0, paddle_width-r0, 0),
				Vector(paddle_length*sin(amplitude*cos(lambda_y)), 0, paddle_length*cos(amplitude*cos(lambda_y))));
		paddle[i].SetPartMass(m_deltap, m_physparams->rho0[0]);
		paddle[i].Fill(paddle[i].GetParts(), br, true);
		add_moving_body(&paddle[i], MB_MOVING);
	    set_body_cg(&paddle[i], make_double3(paddle_origin.x, paddle_origin.y+y_paddle, paddle_origin.z));
		std::cout << "paddle[" << i <<"] defined \n";
		std::cout << "  at y = " << y_paddle << "\n";
		  
		}
	
    std::cout << "paddles defined" << "\n";

	bottom_rect = Rect(Point(h_length, 0, 0), Vector(0, ly, 0),
			Vector(slope_length/cos(beta), 0.0, slope_length*tan(beta)));
	if (!use_bottom_plane) {
	   bottom_rect.SetPartMass(m_deltap, m_physparams->rho0[0]);
	   bottom_rect.Fill(boundary_parts,br,true);
	   std::cout << "bottom rectangle defined" <<"\n";
	   }
 
	Rect fluid;
	float z = 0;
	int n = 0;
//	const float amplitude = mbpaddledata.amplitude;
//	printf("amplitude = %f\n",mbpaddledata.amplitude);
	printf("r0 = %f\n", r0);
	printf("h_length, slope_length: %f, %f\n", h_length, slope_length);
    
	while (z < H) {
	 	z = n*m_deltap + 1.5*r0;    //z = n*m_deltap + 1.5*r0;
	 	std::cout << "z = " <<z <<"\n";
		 
		for (uint i=0; i< npaddles-1; i++) {
			double arg = lambda*paddle_width*(2*i+1)/2;
		    double amplitude = paddle_amplitude*sin(arg);
		   
            float x = paddle_origin.x + (z - paddle_origin.z)*tan(amplitude) + 1.0*r0/cos(amplitude);
 		 
			float l = h_length + z/tan(beta) - 1.5*r0/sin(beta) - x;
			float y = r0 + paddle_origin.y + paddle_width*i; // paddle corner
		//	printf("i, amplitude: %d, %f \n", i, amplitude);
		//	printf("x, y, z, l: %f, %f, %f, %f \n", x, y, z, l);
			if (l <0) l = 0;
			if (l > h_length+ slope_length -x) l = h_length+slope_length -r0 -x;
			fluid = Rect(Point(x, y, z),
				Vector(0, paddle_width-r0, 0), Vector(l, 0, 0));
			fluid.SetPartMass(m_deltap, m_physparams->rho0[0]);
			fluid.Fill(parts, m_deltap, true);			
		   }
		// last paddle is r0 shorter

			double amplitude= paddle_amplitude;
			float x = paddle_origin.x + (z - paddle_origin.z)*tan(amplitude) + 1.0*r0/cos(amplitude);
			float l = h_length + z/tan(beta) - 1.5*r0/sin(beta) - x;
			float y = r0 + paddle_origin.y + paddle_width*(npaddles-1);
	//		printf("i, amplitude: %d, %f \n", npaddles, amplitude);
		//	printf("x, y, z, l: %f, %f, %f, %f \n", x, y, z, l);
			if (l <0) l = 0;
			if (l > h_length+ slope_length -x) l = h_length+slope_length -r0 -x;
		    fluid = Rect(Point(x,y,z), Vector(0, paddle_width-2*r0,0), Vector(l,0,0));
		    fluid.SetPartMass(m_deltap, m_physparams->rho0[0]);
		    fluid.Fill(parts, m_deltap,true);
		 n++;
	}
	 
	 /*
	if (m_simparams.testpoints) {
		Point pos = Point(0.5748, 0.1799, 0.2564, 0.0);
		test_points.push_back(pos);
		pos = Point(0.5748, 0.2799, 0.2564, 0.0);
		test_points.push_back(pos);
		pos = Point(1.5748, 0.2799, 0.2564, 0.0);
		test_points.push_back(pos);
	}
	*/
	if (use_cyl) {
		Point p[10];
		p[0] = Point(h_length + slope_length/(cos(beta)*10), ly/2., 0);
		p[1] = Point(h_length + slope_length/(cos(beta)*10), ly/6.,  0);
		p[2] = Point(h_length + slope_length/(cos(beta)*10), 5*ly/6, 0);
		p[3] = Point(h_length + slope_length/(cos(beta)*5), 0, 0);
		p[4] = Point(h_length + slope_length/(cos(beta)*5), ly/3, 0);
		p[5] = Point(h_length + slope_length/(cos(beta)*5), 2*ly/3, 0);
		p[6] = Point(h_length + slope_length/(cos(beta)*5), ly, 0);
		p[7] = Point(h_length + 3*slope_length/(cos(beta)*10), ly/6, 0);
		p[8] = Point(h_length + 3*slope_length/(cos(beta)*10), ly/2, 0);
		p[9] = Point(h_length+ 3*slope_length/(cos(beta)*10), 5*ly/6, 0);
		p[10] = Point(h_length+ 4*slope_length/(cos(beta)*10), ly/2, 0);

		for (int i = 0; i < 11; i++) {
			cyl[i] = Cylinder(p[i], Vector(.025, 0, 0), Vector(0, 0, height));
			cyl[i].SetPartMass(m_deltap, m_physparams->rho0[0]);
			cyl[i].FillBorder(boundary_parts, br, false, false);
			cyl[i].Unfill(parts, br);
		}
	}
	if (use_cone) {
		Point p1 = Point(h_length + slope_length/(cos(beta)*10), ly/2, 0);
		cone = Cone(p1,Vector(ly/4, 0.0, 0.0), Vector(ly/10, 0., 0.), Vector(0, 0, height));
		cone.SetPartMass(m_deltap, m_physparams->rho0[0]);
		cone.FillBorder(boundary_parts, br, false, true);
		cone.Unfill(parts, br);
    }
	
	return parts.size() + boundary_parts.size() + get_bodies_numparts() + test_points.size();
}


uint WaveTank::fill_planes()
{
    if (!use_bottom_plane) {
		return 5;
		}
	else {
		return 6;
		} //corresponds to number of planes
}


void WaveTank::copy_planes(double4 *planes)
{
	const float w = m_size.y;
	const float l = h_length + slope_length;

	//  plane is defined as a x + by +c z + d= 0
	planes[0] = make_double4(0, 0, 1.0, 0);   //bottom, where the first three numbers are the normal, and the last is d.
	planes[1] = make_double4(0, 1.0, 0, 0);   //wall
	planes[2] = make_double4(0, -1.0, 0, w); //far wall
 	planes[3] = make_double4(1.0, 0, 0, 0);  //end
 	planes[4] = make_double4(-1.0, 0, 0, l);  //one end
 	if (use_bottom_plane)  {
		planes[5] = make_double4(-sin(beta),0,cos(beta), h_length*sin(beta));  //sloping bottom starting at x=h_length
	}
}


 


void WaveTank::copy_to_array(BufferList &buffers)
{
    float4 *pos = buffers.getData<BUFFER_POS>();
	hashKey *hash = buffers.getData<BUFFER_HASH>();
	float4 *vel = buffers.getData<BUFFER_VEL>();
	particleinfo *info = buffers.getData<BUFFER_INFO>();
	
    int j = 0;
    /*
	if (test_points.size()) {
		//Testpoints
		std::cout << "\nTest points: " << test_points.size() << "\n";
		std::cout << "      " << 0  << "--" << test_points.size() << "\n";
		for (uint i = 0; i < test_points.size(); i++) {
			calc_localpos_and_hash(test_points[i], info[i], pos[i], hash[i]);
			vel[i] = make_float4(0, 0, 0, m_physparams->rho0[0]);
			info[i]= make_particleinfo(TESTPOINTSPART, 0, i);  // first is type, object, 3rd id
		}
		std::cout << "Test point mass:" << pos[j-1].w << "\n";
		j += test_points.size();
	}
*/
	std::cout << "\nBoundary parts: " << boundary_parts.size() << "\n";
		std::cout << "      "<< j  <<"--"<< j+ boundary_parts.size() << "\n";    //FIXED: 0 --> j 
	for (uint i = j; i < j + boundary_parts.size(); i++) {
		calc_localpos_and_hash(boundary_parts[i - j], info[i], pos[i], hash[i]);
		vel[i] = make_float4(0, 0, 0, m_physparams->rho0[0]);
		info[i]= make_particleinfo(PT_BOUNDARY, 0, i);  // first is type, object, 3rd id
	}
    j += boundary_parts.size();
	std::cout << "Boundary part mass:" << pos[j-1].w << "\n";

	std::cout << "\nFluid parts: " << parts.size() << "\n";
	std::cout << "      "<< j  <<"--"<< j + parts.size() << "\n";
	for (uint i = j; i < j + parts.size(); i++) {
		calc_localpos_and_hash(parts[i - j], info[i], pos[i], hash[i]);
		float rho = m_physparams->rho0[0]; // density(H - pos[i].z, 0);
		vel[i] = make_float4(0, 0, 0, rho);
	    info[i]= make_particleinfo(PT_FLUID, 0, i);
	}
	j += parts.size();
	std::cout << "Fluid part mass:" << pos[j-1].w << "\n";
	
    for (uint k = 0; k < m_bodies.size(); k++) {
			PointVect & rbparts = m_bodies[k]->object->GetParts();
			std::cout << "Rigid body " << k << ": " << rbparts.size() << " particles ";
			for (uint i = 0; i < rbparts.size(); i++) {
				uint ij = i + j;
				float ht = H - rbparts[i](2);
				if (ht < 0)
					ht = 0.0;
				float rho = density(ht, 0);
				rho = m_physparams->rho0[0];
				vel[ij] = make_float4(0, 0, 0, rho);
				uint ptype = (uint) PT_BOUNDARY;
				switch (m_bodies[k]->type) {
					case MB_ODE:
						ptype |= FG_MOVING_BOUNDARY | FG_COMPUTE_FORCE;
						break;
					case MB_FORCES_MOVING:
						ptype |= FG_COMPUTE_FORCE | FG_MOVING_BOUNDARY;
						break;
					case MB_MOVING:
						ptype |= FG_MOVING_BOUNDARY;
						break;
				}
				info[ij] = make_particleinfo(ptype, k, i );
				calc_localpos_and_hash(rbparts[i], info[ij], pos[ij], hash[ij]);
			}
			j += rbparts.size();
			std::cout << ", part mass: " << pos[j-1].w << "\n";
			std::cout << ", part type: " << type(info[j-1])<< "\n";
		}




	std::cout << "Everything uploaded" <<"\n";
}
 

#undef MK_par
