#include "curvature_codegen.h"
#include <stdio.h>
#include <math.h>
#include <time.h>

const int MAX_CALLS=1e7;

// Call radius and grad repeatedly
void timeit() {
    double dt=2*M_PI/30;
    double r=2.5;
    double s;

    double a[3]= {0.0}, b[3]= {0.0}, c[3]= {0.0};
    a[0]=r;
    b[0]=r*cos(dt);
    b[1]=r*sin(dt);
    c[0]=r*cos(2*dt);
    c[1]=r*sin(2*dt);
    double ga[3]= {0.0}, gb[3]= {0.0}, gc[3]= {0.0};

    double junk=0.0; // Just in case things are optimized out
    clock_t t=clock();

    for (int i=0; i<MAX_CALLS; i++) {
        a[0]=i;
        junk+=curvature_radius(a[0],a[1],a[2],
			       b[0],b[1],b[2],
			       c[0],c[1],c[2]);
        curvature_radius_grad(a[0],a[1],a[2],
			      b[0],b[1],b[2],
			      c[0],c[1],c[2],
			      ga, gb, gc);
        junk+=ga[0];
    }

    t=clock()-t;

    printf ("It took me %f seconds (%f μs/call) \n %lf.\n",
            ((float)t)/CLOCKS_PER_SEC,
            1.0e6*((float)t)/CLOCKS_PER_SEC/MAX_CALLS,junk);

}


// Three tests: R, dR, R(known coordinates)
void test1() {

    double dt=2*M_PI/30;
    double r=2.5;
    double s;

    double a[3]= {0.0}, b[3]= {0.0}, c[3]= {0.0};
    a[0]=r;
    b[0]=r*cos(dt);
    b[1]=r*sin(dt);
    c[0]=r*cos(2*dt);
    c[1]=r*sin(2*dt);

    s=curvature_radius(a[0],a[1],a[2],
		       b[0],b[1],b[2],
		       c[0],c[1],c[2]);

    printf("%f %f %f\n",a[0],a[1],a[2]);
    printf("%f %f %f\n",b[0],b[1],b[2]);
    printf("%f %f %f\n",c[0],c[1],c[2]);

    printf("Rr: %f (expected %f)\n",s,r);
    printf("\n");


    double ga[3]= {0.0}, gb[3]= {0.0}, gc[3]= {0.0};
    curvature_radius_grad(a[0],a[1],a[2],
			  b[0],b[1],b[2],
			  c[0],c[1],c[2],
			  ga, gb, gc);
    printf("Gradient (atan should be 0, %f and 2*%f):\n",dt,dt);
    printf("%f %f %f\n",ga[0],ga[1],ga[2]);
    printf("%f %f %f\n",gb[0],gb[1],gb[2]);
    printf("%f %f %f\n",gc[0],gc[1],gc[2]);
    printf("\n");


    // made up sort-of-hand-checked
    double rs = sqrt(395485.)/181;
    a[0]=1.0;
    a[1]=2.0;
    a[2]=3.0;
    b[0]=2.0;
    b[1]=3.0;
    b[2]=9.0;
    c[0]=5.0;
    c[1]=4.0;
    c[2]=3.0;

    s=curvature_radius(a[0],a[1],a[2],
		       b[0],b[1],b[2],
		       c[0],c[1],c[2]);

    printf("%f %f %f\n",a[0],a[1],a[2]);
    printf("%f %f %f\n",b[0],b[1],b[2]);
    printf("%f %f %f\n",c[0],c[1],c[2]);

    printf("Rr: %f (expected %f)\n",s,rs);

}




int main(int argc, char**argv) {
//    test1();
    timeit();
}
