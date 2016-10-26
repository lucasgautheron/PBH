#include <vector>
#include <stdio.h>
#include <math.h>
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
 

#define PI 3.141592654
#define CEL 299792458
#define GRAV 6.67e-11

#define PBH_MAX_MASS 100
#define PBH_MIN_MASS 1

#define HALO_DENSITY 0.24
#define HALO_RADIUS 10

// m : solar mass
// d : pc
// v : 100 km/s



// ref: https://arxiv.org/pdf/astro-ph/0201102v1.pdf
double merger_cross_section(double m1, double m2, double v_rel)
{
    // prefactor : 2*3*(85*3/(6*sqrt(2)))^(2/7) * (6.674e-11)^2 *(1.99e30)^2 / (299792458)^(10/7) / (100000)^(18/7) / (3e16)^2
    //return 2.0*pow(85*PI/(6*sqrt(2.)) * m1 * m2, 2.0/7.0) * GRAV * GRAV * pow(m1+m2, 10.0/7.0) / (pow(CEL, 10.0/7.0) * pow(v_rel, 18.0/7.0));

    return 3.3554462e-17 * pow(m1+m2, 10.0/7.0) / (pow(CEL, 10.0/7.0) * pow(v_rel, 18.0/7.0));
}

double detector_efficiency(double m1, double m2, double d)
{
    return d < 1000 ? true : false;
}

double halo_mass_profile(double r, double m, double rho, double Rs) // NFW
{
    return rho/((r/Rs) * (1+r/Rs));
}

double pbh_mass_distribution(double m)
{
   return 1.0/(PBH_MAX_MASS-PBH_MIN_MASS);
}

double avg_sigmav(double r, double m1, double m2)
{
    return merger_cross_section(m1, m2, 0.002) * 0.002;
}

Double_t rate_integrand(Double_t *x, Double_t *p)
{
    return 4.0 * PI * x[0] * x[0] * 0.5 * pbh_mass_distribution(x[1]) * pbh_mass_distribution(x[2]) * halo_mass_profile(x[0], x[1], HALO_DENSITY, HALO_RADIUS) * halo_mass_profile(x[0], x[2], HALO_DENSITY, HALO_RADIUS) * avg_sigmav(x[0], x[1], x[2]);
}

double halo_rate()
{
    // Create the function and wrap it
   TF3 f("integrand", rate_integrand, 0, HALO_RADIUS, PBH_MIN_MASS, PBH_MAX_MASS, PBH_MIN_MASS, PBH_MAX_MASS);
   ROOT::Math::WrappedMultiTF1 wf1(f);
   // Create the Integrator
   ROOT::Math::AdaptiveIntegratorMultiDim ig;
   // Set parameters of the integration
   ig.SetFunction(wf1);
   ig.SetRelTolerance(0.001);
   double xmin[] = {0, PBH_MIN_MASS, PBH_MIN_MASS};
   double xmax[] = {HALO_RADIUS, PBH_MAX_MASS, PBH_MAX_MASS};
   return ig.Integral(xmin, xmax);
}

int main()
{
    printf("halo_rate = %e\n", 86400 * 365 * halo_rate());
    return 0;
}
