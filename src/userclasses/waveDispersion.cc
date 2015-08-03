#include <iostream>
#include <cmath>
#include "waveDispersion.h"

waveDispersion::waveDispersion(void) {
	
}

waveDispersion::waveDispersion(double wave_period, double water_depth) {
	T = wave_period;
	H = water_depth;
	
	calcDisp();	
}

waveDispersion::~waveDispersion() {
	return;
}

double waveDispersion::calc(double wave_period, double water_depth) {
	T = wave_period;
	H = water_depth;
	
	calcDisp();
	
	return this->k;
}

void waveDispersion::calcDisp() {
    double sig, fp;

    sig= 2*M_PI / T;
    
    double s2og=sig*sig/g;
    k = s2og / sqrt( tanh( sig*sig * H/g ) );
    for (int i=1;i<6;i++) {
        f = s2og - k * tanh( k*H );
        fp = -tanh( k*H ) - k*H*(1 - tanh( k*H ) * tanh( k*H ) );
        k = k - f/fp;
    }
}