/*	BSFfast-header  --  efficient BSF interpolations
 *  written by cGPT (user: Stefan Lederer)
 * -v1 created on:		2025-11-06
 *******************************************/ 

#ifndef BSFFAST_H
#define BSFFAST_H

#include <stdlib.h>

// Data arrays (globally defined in BSFfast.c)
    //for scalars
extern double **topPlatS_data;
extern double **topCutS_data;
extern double **botPlatS_data;
extern double **botCutS_data;
extern double **noTCutS_data;
extern double **noTPlatS_data;

extern double *xQEDS_data;
extern double *xNoTQEDS_data;
extern double *xNoTQCDS_data;

extern double *topPlatS_m;
extern double *botPlatS_m;
extern double *noTPlatS_m;
extern double *topCutS_m;
extern double *botCutS_m;
extern double *noTCutS_m;

extern double *topPlatS_x;
extern double *botPlatS_x;
extern double *noTPlatS_x;
extern double *topCutS_x;
extern double *botCutS_x;
extern double *noTCutS_x;

extern double xQEDS_m;
extern double xNoTQEDS_m;
extern double xNoTQCDS_m;

extern double *xQEDS_x;
extern double *xNoTQEDS_x;
extern double *xNoTQCDS_x;

extern int topPlatS_x_size;
extern int topPlatS_m_size;
extern int topCutS_x_size;
extern int topCutS_m_size;
extern int botPlatS_x_size;
extern int botPlatS_m_size;
extern int botCutS_x_size;
extern int botCutS_m_size;
extern int noTPlatS_x_size;
extern int noTPlatS_m_size;
extern int noTCutS_x_size;
extern int noTCutS_m_size;

extern int xQEDS_x_size;
extern int xNoTQEDS_x_size;
extern int xNoTQCDS_x_size;

    // for fermions
extern double **topPlatF_data;
extern double **topCutF_data;
extern double **botPlatF_data;
extern double **botCutF_data;
extern double **noTCutF_data;
extern double **noTPlatF_data;

extern double *xQEDF_data;
extern double *xNoTQEDF_data;
extern double *xNoTQCDF_data;

extern double *topPlatF_m;
extern double *botPlatF_m;
extern double *noTPlatF_m;
extern double *topCutF_m;
extern double *botCutF_m;
extern double *noTCutF_m;

extern double *topPlatF_x;
extern double *botPlatF_x;
extern double *noTPlatF_x;
extern double *topCutF_x;
extern double *botCutF_x;
extern double *noTCutF_x;

extern double xQEDF_m;
extern double xNoTQEDF_m;
extern double xNoTQCDF_m;

extern double *xQEDF_x;
extern double *xNoTQEDF_x;
extern double *xNoTQCDF_x;

extern int topPlatF_x_size;
extern int topPlatF_m_size;
extern int topCutF_x_size;
extern int topCutF_m_size;
extern int botPlatF_x_size;
extern int botPlatF_m_size;
extern int botCutF_x_size;
extern int botCutF_m_size;
extern int noTPlatF_x_size;
extern int noTPlatF_m_size;
extern int noTCutF_x_size;
extern int noTCutF_m_size;

extern int xQEDF_x_size;
extern int xNoTQEDF_x_size;
extern int xNoTQCDF_x_size;

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
double BSFfast_sigveff_dQCD_S_const(double alpha, double m, double x);

double BSFfast_sigveff_QCD_FU_plat(double m, double x);
double BSFfast_sigveff_QCD_FU_cut(double m, double x);
double BSFfast_sigveff_QCD_FD_plat(double m, double x);
double BSFfast_sigveff_QCD_FD_cut(double m, double x);
double BSFfast_sigveff_QCD_F_plat(double m, double x);
double BSFfast_sigveff_QCD_F_cut(double m, double x);
double BSFfast_sigveff_QED_F_const(double m, double x);
double BSFfast_sigveff_dQED_F_const(double alpha, double m, double x);
double BSFfast_sigveff_dQED_F_noTr_const(double alpha, double m, double x);
double BSFfast_sigveff_dQCD_F_const(double alpha, double m, double x);

// Library initialization
__attribute__((constructor)) void init_library(void);

#endif // BSFFAST_H

