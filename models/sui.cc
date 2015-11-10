#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double SUIpathLoss(float f, float TxH, float RxH, float d, int mode)
{
	/*
	   f = Frequency (MHz)
	   TxH =  Transmitter height (m)
	   RxH = Receiver height (m)
	   d = distance (km)
	   mode 1 = Hilly + trees
	   mode 2 = Flat + trees OR hilly + light foliage
	   mode 3 = Flat + light foliage
	   http://www.cl.cam.ac.uk/research/dtg/lce-pub/public/vsa23/VTC05_Empirical.pdf
	 */
	d = d * 1000;		// km to m
	if (f < 1900 || f > 11000) {
		printf("Error: SUI model frequency range 1.9-11GHz\n");
		exit(EXIT_FAILURE);
	}
	// Terrain mode A is default
	double a = 4.6;
	double b = 0.0075;
	double c = 12.6;
	double s = 10.6;	// Optional fading value
	int XhCF = -10.8;

	if (mode == 2) {
		a = 4.0;
		b = 0.0065;
		c = 17.1;
		s = 6;		// average
	}
	if (mode == 3) {
		a = 3.6;
		b = 0.005;
		c = 20;
		s = 3;		// Optimistic
		XhCF = -20;
	}
	double d0 = 100;
	double A = 20 * log10((4 * M_PI * d0) / (300 / f));
	double y = a - (b * TxH) + (c / TxH);
	double Xf = 6 * log10(f / 2000);
	double Xh = XhCF * log10(RxH / 20);

	return A + (10 * y) * (log10(d / d0)) + Xf + Xh + s;
}
