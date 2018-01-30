/******************************************************************************
 *                       Code generated with sympy 1.0                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                  This file is part of 'plumed_curvature'                   *
 ******************************************************************************/
#include "curvature_codegen_cse.h"
#include <math.h>

double curvature_radius(double r_1x, double r_1y, double r_1z, double r_2x, double r_2y, double r_2z, double r_3x, double r_3y, double r_3z) {

   double curvature_radius_result;
   curvature_radius_result = (1.0L/2.0L)*sqrt((pow(r_1x - r_3x, 2) + pow(r_1y - r_3y, 2) + pow(r_1z - r_3z, 2))/(-pow((r_1x - r_2x)*(r_2x - r_3x) + (r_1y - r_2y)*(r_2y - r_3y) + (r_1z - r_2z)*(r_2z - r_3z), 2)/((pow(r_1x - r_2x, 2) + pow(r_1y - r_2y, 2) + pow(r_1z - r_2z, 2))*(pow(r_2x - r_3x, 2) + pow(r_2y - r_3y, 2) + pow(r_2z - r_3z, 2))) + 1));
   return curvature_radius_result;

}
void curvature_radius_grad(double r_1x, double r_1y, double r_1z, 
    double r_2x, double r_2y, double r_2z,   
    double r_3x, double r_3y, double r_3z,   
    double *g1, double *g2, double *g3)    { 
	double p0 = -r_3x;
	double p1 = p0 + r_1x;
	double p2 = -r_3y;
	double p3 = p2 + r_1y;
	double p4 = -r_3z;
	double p5 = p4 + r_1z;
	double p6 = pow(p1, 2) + pow(p3, 2) + pow(p5, 2);
	double p7 = r_1x - r_2x;
	double p8 = r_1y - r_2y;
	double p9 = r_1z - r_2z;
	double p10 = 1.0/(pow(p7, 2) + pow(p8, 2) + pow(p9, 2));
	double p11 = p0 + r_2x;
	double p12 = p2 + r_2y;
	double p13 = p4 + r_2z;
	double p14 = 1.0/(pow(p11, 2) + pow(p12, 2) + pow(p13, 2));
	double p15 = p11*p7 + p12*p8 + p13*p9;
	double p16 = 1.0/(p10*p14*pow(p15, 2) - 1);
	double p17 = sqrt(-p16*p6);
	double p18 = (1.0L/2.0L)*p17/p6;
	double p19 = p10*p15;
	double p20 = p19*p7;
	double p21 = p10*p14*p15*p16*p6;
	double p22 = p19*p8;
	double p23 = p19*p9;
	double p24 = p14*p15;
	double p25 = -p11*p24;
	double p26 = (1.0L/2.0L)*p10*p14*p15*p16*p17;
	double p27 = -p12*p24;
	double p28 = -p13*p24;
g1[0] = p18*(p1 - p21*(p11 - p20));
g1[1] = p18*(-p21*(p12 - p22) + p3);
g1[2] = p18*(-p21*(p13 - p23) + p5);
g2[0] = -p26*(p20 + p25 + r_1x - 2*r_2x + r_3x);
g2[1] = -p26*(p22 + p27 + r_1y - 2*r_2y + r_3y);
g2[2] = -p26*(p23 + p28 + r_1z - 2*r_2z + r_3z);
g3[0] = -p18*(p1 - p21*(p25 + p7));
g3[1] = -p18*(-p21*(p27 + p8) + p3);
g3[2] = -p18*(-p21*(p28 + p9) + p5);
}