#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double SUIpathLoss(double f, double TxH, double RxH, double d, int mode)
{
	/*
	   f = Frequency (MHz) 1900 to 11000
	   TxH =  Transmitter height (m)
	   RxH = Receiver height (m)
	   d = distance (km)
	   mode A1 = URBAN / OBSTRUCTED
	   mode B2 = SUBURBAN / PARTIALLY OBSTRUCTED
	   mode C3 = RURAL / OPEN
	   http://www.cl.cam.ac.uk/research/dtg/lce-pub/public/vsa23/VTC05_Empirical.pdf
	 */
	d = d * 1000.0;		// km to m

	// Urban (A1) is default
	double a = 4.6;
	double b = 0.0075;
	double c = 12.6;
	double s = 10.6; // Optional fading value. Max 10.6dB
	double XhCF = -10.8;

	if (mode == 2) { // Suburban
		a = 4.0;
		b = 0.0065;
		c = 17.1;
	}
	if (mode == 3) { // Rural
		a = 3.6;
		b = 0.005;
		c = 20;
		XhCF = -20;
	}
	double d0 = 100;
	double A = 20 * log10((4 * M_PI * d0) / (300.0 / f));
	double y = a - (b * TxH) + (c / TxH);
	//Correction factors
	double Xf = 6.0 * log10(f / 2000);
	double Xh = XhCF * log10(RxH / 2000);

	return A + (10 * y * log10(d / d0)) + Xf + Xh + s;
}
