//
//  waveDispersion.h
//  jonswap-spec
//
//  Created by Munan Xu on 2015-07-09.
//  Copyright 2015 Munan Xu. All rights reserved.
//

class waveDispersion
{
public:
	waveDispersion (void);
	waveDispersion (double wave_period, double water_depth);
	
	double calc(double wave_period, double water_depth);
	
	double getK() { return k; }
	double getF() { return f; }
	
	virtual ~waveDispersion ();
	
private:
	double T, H, k, f;
	const double g = 9.81;
	void calcDisp();
};
