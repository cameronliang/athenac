#ifndef GLOBALS_H
#define GLOBALS_H
/*============================================================================*/
/*! \file globals.h
 *  \brief Contains global variables.
 *
 * PURPOSE: Contains global variables:
 *   The first occurence in this file is included in main.c and defines the
 *   variables.  The second is included everywhere else.		      */
/*============================================================================*/

#ifdef MAIN_C

Real CourNo;                 /*!< Courant, Friedrichs, & Lewy (CFL) number */
#ifdef ISOTHERMAL
Real Iso_csound;             /*!< isothermal sound speed */
Real Iso_csound2;            /*!< isothermal sound speed squared */
#elif defined ADIABATIC
Real Gamma;                  /*!< adiabatic index (ratio of specific heats) */
Real Gamma_1, Gamma_2;       /*!< (Gamma)-1 and (Gamma)-2 */
#endif
int myID_Comm_world; /*!< Rank (proc ID) in MPI_COMM_WORLD, 0 for single proc */
Real d_MIN = 1.e-2; //TINY_NUMBER;    /*!< density floor */

GravPotFun_t ExternalGravPot = NULL;
CoolingFun_t CoolingFunc = NULL;
#ifdef SELF_GRAVITY
Real four_pi_G, grav_mean_rho;    /*!< 4\pi G and mean density in domain */
#endif

#ifdef SHEARING_BOX
GravPotFun_t ShearingBoxPot = NULL;
Real Omega_0, qshear; /*!< orbital freq and shear parameter dln\Omega/dlnr */
enum SS2DCoord ShBoxCoord;
#endif

#ifdef ROTATING_FRAME
Real Omega_0, Rc; /* Rotating frequency, Radius of center of mass */
#endif

#ifdef PARTICLES
TSFun_t     get_ts    = NULL;     /*!< get the stopping time */
WeightFun_t getweight = NULL;     /*!< get weight function */
#endif

#ifdef THERMAL_CONDUCTION
KappaFun_t KappaFun_a = NULL;   /* coeff of thermal conduction */
KappaFun_t KappaFun_i = NULL;
#endif
#ifdef RESISTIVITY
Real eta_Ohm=0.0, Q_Hall=0.0, Q_AD=0.0;        /*!< diffusivities */
Real d_ind;                                    /*!< index: n_e ~ d^(d_ind) */
EtaFun_t get_myeta = NULL;       /*!< function to calculate the diffusivities */
#endif
#ifdef VISCOSITY
NuFun_t NuFun_a = NULL;               /* coeff of viscosity */
NuFun_t NuFun_i = NULL;
#endif
#ifdef STS
int N_STS;			/*!< number of super timesteps */
Real nu_STS;			/*!< parameter controlling the substeps  */
Real STS_dt;			/*!< STS time step */
int ncycle;
#endif

#ifdef CYLINDRICAL
// StaticGravAcc_t x1GravAcc = NULL;
Real *r=NULL, *ri=NULL;
#ifdef FARGO
OrbitalFun_t OrbitalProfile = NULL;
ShearFun_t ShearProfile = NULL;
#endif
#endif

/*----------------------------------------------------------------------------*/
/* definitions included everywhere except main.c  */

#else /* MAIN_C */

extern Real CourNo;
#ifdef ISOTHERMAL
extern Real Iso_csound, Iso_csound2;
#elif defined ADIABATIC
extern Real Gamma, Gamma_1, Gamma_2;
#endif
extern int myID_Comm_world;
extern Real d_MIN;

extern GravPotFun_t ExternalGravPot;
extern CoolingFun_t CoolingFunc;
#ifdef SELF_GRAVITY
extern Real four_pi_G, grav_mean_rho;
#endif

#ifdef SHEARING_BOX
extern GravPotFun_t ShearingBoxPot;
extern Real Omega_0, qshear;
extern enum SS2DCoord ShBoxCoord;
#endif

#ifdef ROTATING_FRAME
extern Real Omega_0, Rc; /* Omega_0: rotating frequency of frame; Rc: Distance between Center of Mass to Origin of Coordinate */
#endif

#ifdef PARTICLES
extern Real alamcoeff, *grrhoa;
extern TSFun_t     get_ts;
extern WeightFun_t getweight;
#endif

#ifdef THERMAL_CONDUCTION
extern KappaFun_t KappaFun_i;
extern KappaFun_t KappaFun_a;
#endif
#ifdef RESISTIVITY
extern Real eta_Ohm, Q_Hall, Q_AD;
extern Real d_ind;
extern EtaFun_t get_myeta;
#endif
#ifdef VISCOSITY
extern NuFun_t NuFun_i, NuFun_a;
#endif
#ifdef STS
extern int N_STS,ncycle;
extern Real nu_STS, STS_dt;
#endif

#ifdef CYLINDRICAL
// extern StaticGravAcc_t x1GravAcc;
extern Real *r, *ri;
#ifdef FARGO
extern OrbitalFun_t OrbitalProfile;
extern ShearFun_t ShearProfile;
#endif
#endif

#endif /* MAIN_C */
#endif /* GLOBALS_H */
