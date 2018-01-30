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
 double *g1, double *g2, double *g3) { 
	double tmp0 = -r_3x;
	double tmp1 = -r_3y;
	double tmp2 = -r_3z;
	double tmp3 = pow(r_1x + tmp0, 2) + pow(r_1y + tmp1, 2) + pow(r_1z + tmp2, 2);
	double tmp4 = r_1x - r_2x;
	double tmp5 = r_2x + tmp0;
	double tmp6 = r_1y - r_2y;
	double tmp7 = r_2y + tmp1;
	double tmp8 = r_1z - r_2z;
	double tmp9 = r_2z + tmp2;
	double tmp10 = tmp4*tmp5 + tmp6*tmp7 + tmp8*tmp9;
	double tmp11 = pow(tmp10, 2);
	double tmp12 = pow(tmp4, 2) + pow(tmp6, 2) + pow(tmp8, 2);
	double tmp13 = 1.0/tmp12;
	double tmp14 = pow(tmp5, 2) + pow(tmp7, 2) + pow(tmp9, 2);
	double tmp15 = 1.0/tmp14;
	double tmp16 = -tmp11*tmp13*tmp15 + 1;
	double tmp17 = 1.0/tmp16;
	double tmp18 = sqrt(tmp17*tmp3);
	double tmp19 = (1.0L/2.0L)*tmp16*tmp18/tmp3;
	double tmp20 = 2*r_1x;
	double tmp21 = 2*r_3x;
	double tmp22 = -tmp21;
	double tmp23 = (1.0L/2.0L)*tmp17;
	double tmp24 = (1.0L/2.0L)*tmp3/pow(tmp16, 2);
	double tmp25 = 2*r_2x;
	double tmp26 = tmp22 + tmp25;
	double tmp27 = tmp10*tmp13*tmp15;
	double tmp28 = -tmp20;
	double tmp29 = tmp25 + tmp28;
	double tmp30 = tmp11*tmp15/pow(tmp12, 2);
	double tmp31 = 2*r_1y;
	double tmp32 = 2*r_3y;
	double tmp33 = -tmp32;
	double tmp34 = 2*r_2y;
	double tmp35 = tmp33 + tmp34;
	double tmp36 = -tmp31;
	double tmp37 = tmp34 + tmp36;
	double tmp38 = 2*r_1z;
	double tmp39 = 2*r_3z;
	double tmp40 = -tmp39;
	double tmp41 = 2*r_2z;
	double tmp42 = tmp40 + tmp41;
	double tmp43 = -tmp38;
	double tmp44 = tmp41 + tmp43;
	double tmp45 = (1.0L/4.0L)*tmp17*tmp18;
	double tmp46 = -tmp25;
	double tmp47 = tmp11*tmp13/pow(tmp14, 2);
	double tmp48 = -tmp34;
	double tmp49 = -tmp41;
g1[0] = tmp19*(tmp23*(tmp20 + tmp22) + tmp24*(tmp26*tmp27 + tmp29*tmp30));
g1[1] = tmp19*(tmp23*(tmp31 + tmp33) + tmp24*(tmp27*tmp35 + tmp30*tmp37));
g1[2] = tmp19*(tmp23*(tmp38 + tmp40) + tmp24*(tmp27*tmp42 + tmp30*tmp44));
g2[0] = tmp45*(tmp27*(-4*r_2x + tmp20 + tmp21) + tmp30*(tmp20 + tmp46) + tmp47*(tmp21 + tmp46));
g2[1] = tmp45*(tmp27*(-4*r_2y + tmp31 + tmp32) + tmp30*(tmp31 + tmp48) + tmp47*(tmp32 + tmp48));
g2[2] = tmp45*(tmp27*(-4*r_2z + tmp38 + tmp39) + tmp30*(tmp38 + tmp49) + tmp47*(tmp39 + tmp49));
g3[0] = tmp19*(tmp23*(tmp21 + tmp28) + tmp24*(tmp26*tmp47 + tmp27*tmp29));
g3[1] = tmp19*(tmp23*(tmp32 + tmp36) + tmp24*(tmp27*tmp37 + tmp35*tmp47));
g3[2] = tmp19*(tmp23*(tmp39 + tmp43) + tmp24*(tmp27*tmp44 + tmp42*tmp47));
}