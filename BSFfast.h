/*	BSFfast-header  --  efficient BSF interpolations
 *  written by cGPT (user: Stefan Lederer)
 * -v1 created on:		2025-11-06
 *******************************************/ 

#ifndef BSFFAST_H
#define BSFFAST_H

#include <stdlib.h>

// Data arrays (globally defined in BSFfast.c)
extern double **topPlat_data;
extern double **topCut_data;
extern double **botPlat_data;
extern double **botCut_data;
extern double **noTCut_data;
extern double **noTPlat_data;

extern double *xQED_data;
extern double *xNoTQED_data;
extern double *xNoTQCD_data;

extern double *topPlat_m;
extern double *botPlat_m;
extern double *noTPlat_m;
extern double *topCut_m;
extern double *botCut_m;
extern double *noTCut_m;

extern double *topPlat_x;
extern double *botPlat_x;
extern double *noTPlat_x;
extern double *topCut_x;
extern double *botCut_x;
extern double *noTCut_x;

extern double xQED_m;
extern double xNoTQED_m;
extern double xNoTQCD_m;

extern double *xQED_x;
extern double *xNoTQED_x;
extern double *xNoTQCD_x;

extern int topPlat_x_size;
extern int topPlat_m_size;
extern int topCut_x_size;
extern int topCut_m_size;
extern int botPlat_x_size;
extern int botPlat_m_size;
extern int botCut_x_size;
extern int botCut_m_size;
extern int noTPlat_x_size;
extern int noTPlat_m_size;
extern int noTCut_x_size;
extern int noTCut_m_size;

extern int xQED_x_size;
extern int xNoTQED_x_size;
extern int xNoTQCD_x_size;

// Function declarations

// Data loading functions
//int BSFfast_loadData(const char *filename, double **m, double **x, double **f);
//void BSFfast_importData_2d(const char *filename, double **m, double **x, double ***data, int *msize, int *xsize);
//void BSFfast_importData_1d(const char *filename, double *mUnique, double **x, double **data, int *xsize);

// Interpolation functions
double BSFfast_linear_interpolate(double *x_unique, double *data_1d, double x, int x_size);
double BSFfast_bilinear_interpolate(double *m_unique, double *x_unique, double **data_2d, double m, double x, int m_size, int x_size);

// Specific data functions for evaluating on loaded data
	//collective interface function
double BSFfast_sigveff_QCD_SU_plat(double m, double x);
double BSFfast_sigveff_QCD_SU_cut(double m, double x);
double BSFfast_sigveff_QCD_SD_plat(double m, double x);
double BSFfast_sigveff_QCD_SD_cut(double m, double x);
double BSFfast_sigveff_QCD_S_plat(double m, double x);
double BSFfast_sigveff_QCD_S_cut(double m, double x);
double BSFfast_sigveff_QED_S_const(double m, double x);
double BSFfast_sigveff_dQED_S_const(double alpha, double m, double x);
double BSFfast_sigveff_dQED_S_noTr_const(double alpha, double m, double x);
double BSFfast_sigveff_dSU3_S_const(double alpha, double m, double x);

// Library initialization
__attribute__((constructor)) void init_library(void);

#endif // BSFFAST_H

