#ifndef  GAUNT_FACTOR_H
#define  GAUNT_FACTOR_H
#include<gsl/gsl_spline.h>
#include"constants.h"
/*! \file gaunt_factor.h
 *  \brief Function declarations for routines for calculating
 *  free-free and bound-free Gaunt factors.
 *
 *  The free-free Gaunt factors follow Sutherland 1998, MNRAS, 300, 321-330.
 *  
 *  http://adsabs.harvard.edu/abs/1998MNRAS.300..321S
 *
 *  The bound-free Gaunt factors follow Karzas and Latter 1961, ApJS, 6, 167
 * 
 *  http://adsabs.harvard.edu/abs/1961ApJS....6..167K
 *
 *
 *   Note this code only works at E*Ry >~ 1.0e-3.
 *   Small energies cause instability in the sigma_nl_l->E's.
 */

//note this code only works at E*Ry >~ 1.0e-3
//small energies cause instability in the sigma_nl_l->E's


/*! \class Constants C_gaunt;
 *  \brief A Constants instance used for convenience.
 */
extern Constants C_gaunt;


/*! \var double ggff[41]
 *  \brief "gamma^2" table used for free-free gaunt factor.  See Sutherland 1998. */
extern double ggff[41];

/*! \var double ugff[81]
 *  \brief "u" table used for free-free gaunt factor.  See Sutherland 1998. */
extern double ugff[81];

/*! \var double gffgu[3321]
 *  \brief 2-D tabulated free-free Gaunt factor from Sutherland 1998.
 *  As a function of "u" and "gamma^2".  */
extern double gffgu[3321];


/*! \fn double g2_int[41];
 *  \brief "gamma^2" array for use with interpolated velocity-averaged free-free Gaunt factor 
 *   from Sutherland 1998*/
extern double g2_int[41];

/*! \fn double gffint[41]
 *  \brief Interpolated velocity-averaged free-free Gaunt factor from Sutherland 1998 */
extern double gffint[41];

/*! \var int gffint_flag
 *  \brief Flag to indicate whether the gff spline interpolant has been initialized. */
extern int gffint_flag;
/*! \var gsl_spline *gff_spline
 *  \brief GSL spline for velocity-averaged free-free Gaunt factor. */
extern gsl_spline *gff_spline;
/*! \var gsl_interp_accel *gff_accel
 *  \brief GSL interp accel for velocity-averaged free-free Gaunt factor. */
extern gsl_interp_accel *gff_acc;

/*! \fn void initialize_g_ff_u(void)
 *  \brief Initialization function for a spline interpolation of the
 *         frequency averaged free-free Gaunt factor.*/
void initialize_g_ff_u(void);

/*! \fn double acot(double x) 
 *  \brief Simple arccot function.
 */
double acot(double x);


/*! \fn double g_ff_u(double T, double Z)
 *  \brief Total energy free-free Gaunt factor,
 *  averaged over frequency in the Maxwell-Boltzmann distribution
 *  of electron velocities. */
double g_ff_u(double T, double Z);

/*! \fn double g_ff(double T, double Z, double nu)
 *  \brief Returns the free-free Gaunt factor as a
 *  function of temperature, ion charge, and frequency. 
 *  THIS free-free Gaunt factor is the most commonly used.*/
double g_ff(double T, double Z, double nu);

/*! \fn double g_bf(int n, int l, double E, double Z)
 *  \brief Bound-free Gaunt factor calculated following
 *  Karzas and Latter 1961. */
double g_bf(int n, int l, double E, double z);

/*! \fn double g_bf_ave(int n, double E, double Z)
 *  \brief The shell (l) - averaged bound-free Gaunt factor.*/
double g_bf_ave(int n, double E, double z); //shell-averaged

/*! \fn double b_s_G_l(int s, int l, int m, double eta, double rho)
 *  \brief Coefficients b_s for the polynomial solution to the 
 *  bound-free transition matrix element G_l.  */
double b_s_G_l(int s, int l, int m, double eta, double rho);

/*! \fn double sigma_kramers(int n, double E, double Z)
 *  \brief Kramer's semi-classical bound-free cross-section. */
double sigma_kramers(int n, double E, double Z);

/*! \fn double sigma_nlelm1(int n, int l, double E, double Z)
 *  \brief Cross-section for bound-free absorption from state
 *   n,l through a dipole transition to free energy E with angular
 *   momentum l-1 */
double sigma_nlelm1(int n, int l, double E, double Z);

/*! \fn double sigma_nlelp1(int n, int l, double E, double Z)
 *  \brief Cross-section for bound-free absorption from state
 *   n,l through a dipole transition to free energy E with angular
 *   momentum l+1*/
double sigma_nlelp1(int n, int l, double E, double Z);
/*! \fn double G_l(int l, int m, double eta, double rho)
 *  \brief The bound-free transition matrix element G_l. */
double G_l(int l, int m, double eta, double rho);

/*! \fn double g_ff_interpol(double u, double gamma2)
 *  \brief Returns the interpolated free-free Gaunt factor
 *   from the calculations of Sutherland 1998. */
double g_ff_interpol(double u, double gamma2);

#endif //GAUNT_FACTOR_H
