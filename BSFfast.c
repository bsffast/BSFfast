/*written by Stefan Lederer
 * 
 *			BSFfast  --  sigv_eff interpolations for BSF
 * 
 * -version: 2
 * -last modified on:	2025-12-19
 * -v2 created on:		2025-12-19
 *******************************************/ 


// BSFfast.c
#include "BSFfast.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define sigBSFeff_alpha0 0.1
#define BSFfast_Eps 1e-9  //to be used in IsRoughly

    // scalars

double **topPlatS_data = NULL;
double **topCutS_data = NULL;
double **botPlatS_data = NULL;
double **botCutS_data = NULL;
double **noTCutS_data = NULL;
double **noTPlatS_data = NULL;

double *xQEDS_data = NULL;
double *xNoTQEDS_data = NULL;
double *xNoTQCDS_data = NULL;

double *topPlatS_m = NULL;
double *botPlatS_m = NULL;
double *noTPlatS_m = NULL;
double *topCutS_m = NULL;
double *botCutS_m = NULL;
double *noTCutS_m = NULL;

double *topPlatS_x = NULL;
double *botPlatS_x = NULL;
double *noTPlatS_x = NULL;
double *topCutS_x = NULL;
double *botCutS_x = NULL;
double *noTCutS_x = NULL;

double xQEDS_m;
double xNoTQEDS_m;
double xNoTQCDS_m; 

double *xQEDS_x = NULL;
double *xNoTQEDS_x = NULL;
double *xNoTQCDS_x = NULL;

int topPlatS_x_size = 0;
int topPlatS_m_size = 0;
int topCutS_x_size = 0;
int topCutS_m_size = 0;
int botPlatS_x_size = 0;
int botPlatS_m_size = 0;
int botCutS_x_size = 0;
int botCutS_m_size = 0;
int noTPlatS_x_size = 0;
int noTPlatS_m_size = 0;
int noTCutS_x_size = 0;
int noTCutS_m_size = 0;

int xQEDS_x_size = 0;
int xNoTQEDS_x_size = 0;
int xNoTQCDS_x_size = 0;

    // fermions

double **topPlatF_data = NULL;
double **topCutF_data = NULL;
double **botPlatF_data = NULL;
double **botCutF_data = NULL;
double **noTCutF_data = NULL;
double **noTPlatF_data = NULL;

double *xQEDF_data = NULL;
double *xNoTQEDF_data = NULL;
double *xNoTQCDF_data = NULL;

double *topPlatF_m = NULL;
double *botPlatF_m = NULL;
double *noTPlatF_m = NULL;
double *topCutF_m = NULL;
double *botCutF_m = NULL;
double *noTCutF_m = NULL;

double *topPlatF_x = NULL;
double *botPlatF_x = NULL;
double *noTPlatF_x = NULL;
double *topCutF_x = NULL;
double *botCutF_x = NULL;
double *noTCutF_x = NULL;

double xQEDF_m;
double xNoTQEDF_m;
double xNoTQCDF_m; 

double *xQEDF_x = NULL;
double *xNoTQEDF_x = NULL;
double *xNoTQCDF_x = NULL;

int topPlatF_x_size = 0;
int topPlatF_m_size = 0;
int topCutF_x_size = 0;
int topCutF_m_size = 0;
int botPlatF_x_size = 0;
int botPlatF_m_size = 0;
int botCutF_x_size = 0;
int botCutF_m_size = 0;
int noTPlatF_x_size = 0;
int noTPlatF_m_size = 0;
int noTCutF_x_size = 0;
int noTCutF_m_size = 0;

int xQEDF_x_size = 0;
int xNoTQEDF_x_size = 0;
int xNoTQCDF_x_size = 0;


//
bool IsRoughly(double a, double b) {
    return fabs(a - b) <= BSFfast_Eps;
}

// write import function
static int BSFfast_loadData(const char *filename, double **m, double **x, double **f /*, int * msize, int * xsize */) {
	//*msize = 0;
	//*xsize = 0;
	
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(1);
    }
	//read out the numer of lines right now to spare infinite amounts of pain from dynamic memory management
	int line_count = 0;
	char ch;
	while (!feof(file)) {
		ch = fgetc(file);
		if (ch == '\n') {
			line_count++;
		}
	}
	rewind(file); // Reset file pointer to the start
	
    // Capacity set by line count
    int capacity = line_count; 
    double m_value, x_value, f_value;
    int data_line = 0;

    // Allocate memory for m_vals, x_vals, and f_vals
    *m = malloc(capacity * sizeof(double));
    *x = malloc(capacity * sizeof(double));
    *f = malloc(capacity * sizeof(double));
    
    /*
    *msize += capacity;
    *xsize += capacity;*/
    
    // Check for allocation failure
    // Check for allocation failure
	if (*m == NULL || *x == NULL || *f == NULL) {
		fprintf(stderr, "Memory allocation failed for m, x, or f\n");
		exit(1);
	}
    
    // Read data from the file in the form of (m, x, f(m,x))
    while ( fscanf(file, "%lf,%lf,%lf", &m_value, &x_value, &f_value) == 3 ) {
		//// Ensure we don't access memory past capacity
		//if (data_line >= capacity) {
		//	fprintf(stderr, "❌ Overflow detected: data_line=%d capacity=%d\n", data_line, capacity);
		//	break;
		//}
		
        // Store m and x values, fill f(m, x) array in 2D
        (*m)[data_line] = log(m_value);
        (*x)[data_line] = log(x_value);
        (*f)[data_line] = log(f_value);
        data_line++;
    }
    
    fclose(file);
    //returns length of data being the same as sizeof(m)/sizeof(m[0])
    //printf("file %s: read line count = %i\n",filename,data_line);
    //printf("Total lines read: %d (expected %d)\n", data_line, capacity);

    return data_line;
}
 
// write import 2D function
static void BSFfast_importData_2d(const char *filename, double **m, double **x, double ***data, int * msize, int * xsize){
	double *sv_vals = NULL;
	int n = BSFfast_loadData(filename, m, x, &sv_vals);
	//fprintf(stderr,"done importing %s with n=%i\n",filename,n);
	
	// Step 1: Reshuffle unique m and x values to the front.
    int m_size = 0, x_size = 0;
    
		// Find unique m values
    for (int i = 0; i < n; i++) {
		//fprintf(stderr,"i=%i, n=%i, m_size=%i\n",i,n,m_size); 
		
		int found = 0;
        for (int j = 0; j < m_size; j++) {
            if (IsRoughly((*m)[i] , (*m)[j] )) {
                found = 1;
                break;
            }
        }
        if (!found) {
            (*m)[m_size] = (*m)[i];
            m_size++;
        }
        //else{if(10>i || i>n-10){fprintf(stderr,"m value detected. keeping m_size=%i\n",m_size);}}
    }
		// Find unique x values
	//fprintf(stderr, "uniquifing m done. going to x");
    for (int i = 0; i < n; i++) {
		
		int found = 0;
        for (int j = 0; j < x_size; j++) {
            if (IsRoughly( (*x)[i] , (*x)[j])) {
                found = 1;
                break;
            }
        }
        if (!found) {
            (*x)[x_size] = (*x)[i];
            x_size++;
        }
        //else{if(10>i || i>n-10){fprintf(stderr,"x value detected. keeping x_size=%i\n",x_size);}}
    }
	
    (*msize)=m_size;
    (*xsize)=x_size;

    //fprintf(stderr,"msize , xsize after counting uniques = %i , %i\n",(*msize),(*xsize));
    //one could the m,x unique-ification and 2d-construction in one step but this way we know the required size of data up front which is nice.
    // Step 2: Reshape 2D data
    *data = (double **)malloc(m_size * sizeof(double *));
	for (int i = 0; i < m_size; i++) {
		(*data)[i] = (double *)malloc(x_size * sizeof(double));
	}
/*
	data = (double **)malloc(m_size * sizeof(double **));
    for (int i = 0; i < m_size; i++) {
        data[i] = (double **)calloc(x_size, sizeof(double *));
    }*/
	int m_idx = 0;
	int x_idx = 0;
	for (int pos = 0; pos < n; pos++) {
		//find m-position
		for (int i = 0; i < m_size; i++) {
            if (IsRoughly( (*m)[i] , (*m)[pos]) ) {
                m_idx = i;
                break;
            }
		}
		//find x-position
		for (int i = 0; i < x_size; i++) {
            if (IsRoughly( (*x)[i] , (*x)[pos]) ) {
                x_idx = i;
                break;
            }
		}
		
		(*data)[m_idx][x_idx] = sv_vals[pos];
	}
	free(sv_vals);
}

//write 1D import function
static void BSFfast_importData_1d(const char *filename, double *mUnique, double **x, double **data, int * xsize){
	double *m_list = NULL;
	//assumes all m are identical and x are unique. no need for 2d-grid and no need for uniqui-fication 
	*xsize = BSFfast_loadData(filename, &m_list, x, data);
	//since there's no need for interpolation, store the non-log value
	*mUnique = exp(m_list[0]);
}

// write 1D interpolation function
double BSFfast_linear_interpolate(double *x_unique, double *data_1d, double x, int x_size){
	// Find surrounding grid points
    int i_x1 = 0, i_x2 = 0;
    //int x_size = sizeof(x_unique)/sizeof(x_unique[0]);  //need to pass as argument
    
    // Find x1 and x2 (the closest x values)
    for (int j = 0; j < x_size - 2; j++) { 
        if ( x_unique[j+1] >= x ) {  //we start from lowest x, so this suffices and implements linar extrapolation at the same time.
		//without assuming ordered x-data: x_unique[j] <= x && x_unique[j + 1] >= x){ i_x1=j;i_x2=j+1;break;
            i_x1 = j;
            i_x2 = j + 1;
            break;
        }
    }
    if (i_x2==0){ 
		//if we didn't find anything, we must extrapolate to large x
		i_x1 = x_size-2; 
		i_x2 = x_size-1;
		printf( "BSFfast: using linear extrapolation at large log(x)  = %.5f\n",x);
	}  
	else { if (i_x2==1){ printf( "BSFfast: using linear extrapolation at small log(x)  = %f\n",x);} }
    // fprintf(stderr,"identified neighours: x={%i, %i} \n",i_x1,i_x2);
    
    double result = data_1d[i_x1] + (x - x_unique[i_x1]) / (x_unique[i_x2] - x_unique[i_x1]) * (data_1d[i_x2] - data_1d[i_x1]);
	return result;
}

// write 2D interpolation function
double BSFfast_bilinear_interpolate(double *m_unique, double *x_unique, double **data_2d , double m, double x, int m_size, int x_size) {
    // Find surrounding grid points
		//NOTE: I do not sort the m- & x-vectors
    int i_m1 = -1, i_m2 = -1, i_x1 = -1, i_x2 = -1;
    //int m_size = sizeof(m_unique) / sizeof(m_unique[0]); 
    //int x_size = sizeof(x_unique) / sizeof(x_unique[0]);
    
    // Find m1 and m2 (the closest m values)
    for (int i = 0; i < m_size - 1; i++) {
		//ASSUME THE DATA IS SORTED ALREADY
        if (m_unique[i] <= m && m_unique[i + 1] >= m) {
            i_m1 = i;
            i_m2 = i + 1;
            break;
        }
    }
    
    // Find x1 and x2 (the closest x values)
    for (int j = 0; j < x_size - 1; j++) {
		//ASSUME THE DATA IS SORTED ALREADY
        if (x_unique[j] <= x && x_unique[j + 1] >= x) {
            i_x1 = j;
            i_x2 = j + 1;
            break;
        }
    }
    
    if( i_m1 == -1 || i_m2 == -1 || i_x1 == -1 || i_x2 == -1){
		fprintf(stderr,"m=%e or x=%e are likely out-of- or at-bounds!",exp(m),exp(x));
		return -1.;
		}
    
    //fprintf(stderr,"identified neighours: x={%i, %i}  ;  m={%i, %i}\n",i_x1,i_x2,i_m1,i_m2);
    
    double x1 = x_unique[i_x1];
    double x2 = x_unique[i_x2];
    double m1 = m_unique[i_m1];
    double m2 = m_unique[i_m2];
	double q11 = data_2d[i_m1][i_x1];
    double q12 = data_2d[i_m2][i_x1];
    double q21 = data_2d[i_m1][i_x2];
    double q22 = data_2d[i_m2][i_x2];

    double result = (
			 q11*(x2-x)*(m2-m)
           + q21*(x-x1)*(m2-m)
           + q12*(x2-x)*(m-m1)
           + q22*(x-x1)*(m-m1) )
           / ((x2-x1)*(m2-m1));
    /*// Interpolate in the x direction for both m1 and m2
    double slope =  (data_2d[i_m1][i_x2] - data_2d[i_m1][i_x1]) / (x_unique[i_x2] - x_unique[i_x1]);
    double f_m1_x1 = data_2d[i_m1][i_x1] + (x - x_unique[i_x1]) * slope;
    double f_m2_x1 = data_2d[i_m1][i_x2] + (x - x_unique[i_x1]) * slope;

    slope = (data_2d[i_m2][i_x2] - data_2d[i_m2][i_x1]) / (x_unique[i_x2] - x_unique[i_x1]);
    double f_m2_x1 = data_2d[i_m2][i_x1] + (x - x_unique[i_x1]) * slope;
    double f_m2_x2 = data_2d[i_m2][i_x2] + (x - x_unique[i_x1]) * slope;

    // Finally, interpolate in the m direction
    slope = (m - m_unique[i_m1]) / (m_unique[i_m2] - m_unique[i_m1]);
    double f = f_m1_x1 + slope * (f_m2_x1 - f_m1_x1);
    double result =  f + slope * (f_m2_x2 - f_m1_x2);*/

    return result;
}

// define basic functions to evaluate on data
    //scalars

	//TOP Pleteau
double BSFfast_sigveff_QCD_SU_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	//not handling extrapolation for QCD models.
	if (logm>topPlatS_m[topPlatS_m_size-1] || logm<topPlatS_m[0] || logx>topPlatS_x[topPlatS_x_size-1] || logx<topPlatS_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = exp( BSFfast_bilinear_interpolate(topPlatS_m, topPlatS_x,topPlatS_data,log(m),log(x), topPlatS_m_size, topPlatS_x_size) );
    return out;
    }
	//TOP Cutoff
double BSFfast_sigveff_QCD_SU_cut(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>topCutS_m[topCutS_m_size-1] || logm<topCutS_m[0] || logx>topCutS_x[topCutS_x_size-1] || logx<topCutS_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(topCutS_m, topCutS_x,topCutS_data,log(m),log(x), topCutS_m_size, topCutS_x_size) ));
    return out;
    }

	//BOTTOM Plateau
double BSFfast_sigveff_QCD_SD_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>botPlatS_m[botPlatS_m_size-1] || logm<botPlatS_m[0] || logx>botPlatS_x[botPlatS_x_size-1] || logx<botPlatS_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
    double out = (exp( BSFfast_bilinear_interpolate(botPlatS_m, botPlatS_x,botPlatS_data,log(m),log(x), botPlatS_m_size, botPlatS_x_size) ));
    return out;
    }
	//BOTTOM Cutoff
double BSFfast_sigveff_QCD_SD_cut(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>botCutS_m[botCutS_m_size-1] || logm<botCutS_m[0] || logx>botCutS_x[botCutS_x_size-1] || logx<botCutS_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(botCutS_m, botCutS_x,botCutS_data,log(m),log(x), botCutS_m_size, botCutS_x_size) ));
    return out;
    }

	//no-Trans Plateau
double BSFfast_sigveff_QCD_S_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>noTPlatS_m[noTPlatS_m_size-1] || logm<noTPlatS_m[0] || logx>noTPlatS_x[noTPlatS_x_size-1] || logx<noTPlatS_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(noTPlatS_m, noTPlatS_x, noTPlatS_data,log(m),log(x), noTPlatS_m_size, noTPlatS_x_size) ));
    return out;}
	//no-Trans Cutoff
double BSFfast_sigveff_QCD_S_cut(double m, double x) {
	double logm = log(m), logx=log(x);
	if (logm>noTCutS_m[noTCutS_m_size-1] || logm<noTCutS_m[0] || logx>noTCutS_x[noTCutS_x_size-1] || logx<noTCutS_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(noTCutS_m, noTCutS_x,noTCutS_data,log(m),log(x), noTCutS_m_size, noTCutS_x_size) ));
    return out;
}

	//SM-QED-const
double BSFfast_sigveff_QED_S_const(double m, double x) {
	return BSFfast_sigveff_dQED_S_const(1/128.9, m, x);	
}
	//xQED-const
double BSFfast_sigveff_dQED_S_const(double alpha, double m, double x) {
	double r = alpha / sigBSFeff_alpha0; 
	double logx = log(x) + 2 * log(r);
	//DO extrapolate for constant coupling models.
    /*if ( logx > xQED_x[ xQED_x_size - 1 ] || logx < xQED_x[0] ){
		printf("OOB: %.5f > logx=%.5f > %.5f\n", xQED_x[xQED_x_size-1], logx, xQED_x[1]); 
		return -1.;
	}*/
	
    double prefactor = xQEDS_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp( BSFfast_linear_interpolate( xQEDS_x, xQEDS_data, logx , xQEDS_x_size)  );
    return prefactor;
}
	//xQED-noTrans-const
double BSFfast_sigveff_dQED_SnT_const(double alpha, double m, double x) {
	double r = alpha / sigBSFeff_alpha0; 
	double logx = log(x) + 2 * log(r);
    /*if ( logx>xNoTQED_x[ xNoTQED_x_size - 1 ] || logx<xNoTQED_x[0] ){
		printf("OOB: %.5f > logx=%.5f > %.5f\n", xNoTQED_x[xNoTQED_x_size-1], logx, xNoTQED_x[1]); 
		return -1.;
	}*/
	 
    double prefactor = xNoTQEDS_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp( BSFfast_linear_interpolate(xNoTQEDS_x, xNoTQEDS_data, logx , xNoTQEDS_x_size) 	);
    return prefactor;
}
	// xQCD-NoTrans-const
double BSFfast_sigveff_dQCD_S_const(double alpha, double m, double x) {
	double r = alpha / sigBSFeff_alpha0; 
	double logx = log(x) + 2 * log(r);
	/*if ( logx>xNoTQCD_x[ xNoTQCD_x_size - 1 ] || logx<xNoTQCD_x[0] ){
		printf("OOB: %.5f > logx=%.5f > %.5f\n", xNoTQCD_x[xNoTQCD_x_size-1], logx, xNoTQCD_x[1]); 
		return -1.;
	}*/
		
    double prefactor = xNoTQCDS_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp(BSFfast_linear_interpolate(xNoTQCDS_x, xNoTQCDS_data, logx , xNoTQCDS_x_size) );
    return prefactor;
}


    //fermions

	//TOP Pleteau
double BSFfast_sigveff_QCD_FU_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	//not handling extrapolation for QCD models.
	if (logm>topPlatF_m[topPlatF_m_size-1] || logm<topPlatF_m[0] || logx>topPlatF_x[topPlatF_x_size-1] || logx<topPlatF_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = exp( BSFfast_bilinear_interpolate(topPlatF_m, topPlatF_x,topPlatF_data,log(m),log(x), topPlatF_m_size, topPlatF_x_size) );
    return out;
    }
	//TOP Cutoff
double BSFfast_sigveff_QCD_FU_cut(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>topCutF_m[topCutF_m_size-1] || logm<topCutF_m[0] || logx>topCutF_x[topCutF_x_size-1] || logx<topCutF_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(topCutF_m, topCutF_x,topCutF_data,log(m),log(x), topCutF_m_size, topCutF_x_size) ));
    return out;
    }

	//BOTTOM Plateau
double BSFfast_sigveff_QCD_FD_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>botPlatF_m[botPlatF_m_size-1] || logm<botPlatF_m[0] || logx>botPlatF_x[botPlatF_x_size-1] || logx<botPlatF_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
    double out = (exp( BSFfast_bilinear_interpolate(botPlatF_m, botPlatF_x,botPlatF_data,log(m),log(x), botPlatF_m_size, botPlatF_x_size) ));
    return out;
    }
	//BOTTOM Cutoff
double BSFfast_sigveff_QCD_FD_cut(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>botCutF_m[botCutF_m_size-1] || logm<botCutF_m[0] || logx>botCutF_x[botCutF_x_size-1] || logx<botCutF_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(botCutF_m, botCutF_x,botCutF_data,log(m),log(x), botCutF_m_size, botCutF_x_size) ));
    return out;
    }

	//no-Trans Plateau
double BSFfast_sigveff_QCD_F_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>noTPlatF_m[noTPlatF_m_size-1] || logm<noTPlatF_m[0] || logx>noTPlatF_x[noTPlatF_x_size-1] || logx<noTPlatF_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(noTPlatF_m, noTPlatF_x, noTPlatF_data,log(m),log(x), noTPlatF_m_size, noTPlatF_x_size) ));
    return out;}
	//no-Trans Cutoff
double BSFfast_sigveff_QCD_F_cut(double m, double x) {
	double logm = log(m), logx=log(x);
	if (logm>noTCutF_m[noTCutF_m_size-1] || logm<noTCutF_m[0] || logx>noTCutF_x[noTCutF_x_size-1] || logx<noTCutF_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(noTCutF_m, noTCutF_x,noTCutF_data,log(m),log(x), noTCutF_m_size, noTCutF_x_size) ));
    return out;
}

	//SM-QED-const
double BSFfast_sigveff_QED_F_const(double m, double x) {
	return BSFfast_sigveff_dQED_F_const(1/128.9, m, x);	
}
	//xQED-const
double BSFfast_sigveff_dQED_F_const(double alpha, double m, double x) {
	double r = alpha / sigBSFeff_alpha0; 
	double logx = log(x) + 2 * log(r);
	//DO extrapolate for constant coupling models.
    /*if ( logx > xQED_x[ xQED_x_size - 1 ] || logx < xQED_x[0] ){
		printf("OOB: %.5f > logx=%.5f > %.5f\n", xQED_x[xQED_x_size-1], logx, xQED_x[1]); 
		return -1.;
	}*/
	
    double prefactor = xQEDF_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp( BSFfast_linear_interpolate( xQEDF_x, xQEDF_data, logx , xQEDF_x_size)  );
    return prefactor;
}
	//xQED-noTrans-const
double BSFfast_sigveff_dQED_FnT_const(double alpha, double m, double x) {
	double r = alpha / sigBSFeff_alpha0; 
	double logx = log(x) + 2 * log(r);
    /*if ( logx>xNoTQED_x[ xNoTQED_x_size - 1 ] || logx<xNoTQED_x[0] ){
		printf("OOB: %.5f > logx=%.5f > %.5f\n", xNoTQED_x[xNoTQED_x_size-1], logx, xNoTQED_x[1]); 
		return -1.;
	}*/
	 
    double prefactor = xNoTQEDF_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp( BSFfast_linear_interpolate(xNoTQEDF_x, xNoTQEDF_data, logx , xNoTQEDF_x_size) 	);
    return prefactor;
}
	// xQCD-NoTrans-const
double BSFfast_sigveff_dQCD_F_const(double alpha, double m, double x) {
	double r = alpha / sigBSFeff_alpha0; 
	double logx = log(x) + 2 * log(r);
	/*if ( logx>xNoTQCD_x[ xNoTQCD_x_size - 1 ] || logx<xNoTQCD_x[0] ){
		printf("OOB: %.5f > logx=%.5f > %.5f\n", xNoTQCD_x[xNoTQCD_x_size-1], logx, xNoTQCD_x[1]); 
		return -1.;
	}*/
		
    double prefactor = xNoTQCDF_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp(BSFfast_linear_interpolate(xNoTQCDF_x, xNoTQCDF_data, logx , xNoTQCDF_x_size) );
    return prefactor;
}




__attribute__((constructor))
void init_library(void) {
	// import all 9 data files into the global variables
		//function calls are still lacking the data-length information to be passed on !

    	//scalars
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-SU_plateau.csv", &topPlatS_m, &topPlatS_x, &topPlatS_data, &topPlatS_m_size, &topPlatS_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-SU_cutoff.csv", &topCutS_m, &topCutS_x, &topCutS_data, &topCutS_m_size, &topCutS_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-SD_plateau.csv", &botPlatS_m, &botPlatS_x, &botPlatS_data, &botPlatS_m_size, &botPlatS_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-SD_cutoff.csv", &botCutS_m, &botCutS_x, &botCutS_data, &botCutS_m_size, &botCutS_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-S_plateau.csv", &noTPlatS_m, &noTPlatS_x, &noTPlatS_data, &noTPlatS_m_size, &noTPlatS_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-S_cutoff.csv", &noTCutS_m, &noTCutS_x, &noTCutS_data, &noTCutS_m_size, &noTCutS_x_size);

	BSFfast_importData_1d("BSFfast_DataCSV/dQED-S.csv", &xQEDS_m, &xQEDS_x, &xQEDS_data, &xQEDS_x_size);
	BSFfast_importData_1d("BSFfast_DataCSV/dQED-SnoTr.csv", &xNoTQEDS_m, &xNoTQEDS_x, &xNoTQEDS_data, &xNoTQEDS_x_size);
	BSFfast_importData_1d("BSFfast_DataCSV/dQCD-S.csv", &xNoTQEDS_m, &xNoTQEDS_x, &xNoTQEDS_data, &xNoTQCDS_x_size);
	
	    //fermions
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-FU_plateau.csv", &topPlatF_m, &topPlatF_x, &topPlatF_data, &topPlatF_m_size, &topPlatF_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-FU_cutoff.csv", &topCutF_m, &topCutF_x, &topCutF_data, &topCutF_m_size, &topCutF_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-FD_plateau.csv", &botPlatF_m, &botPlatF_x, &botPlatF_data, &botPlatF_m_size, &botPlatF_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-FD_cutoff.csv", &botCutF_m, &botCutF_x, &botCutF_data, &botCutF_m_size, &botCutF_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-F_plateau.csv", &noTPlatF_m, &noTPlatF_x, &noTPlatF_data, &noTPlatF_m_size, &noTPlatF_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/QCD-F_cutoff.csv", &noTCutF_m, &noTCutF_x, &noTCutF_data, &noTCutF_m_size, &noTCutF_x_size);

	BSFfast_importData_1d("BSFfast_DataCSV/dQED-F.csv", &xQEDF_m, &xQEDF_x, &xQEDF_data, &xQEDF_x_size);
	BSFfast_importData_1d("BSFfast_DataCSV/dQED-FnoTr.csv", &xNoTQEDF_m, &xNoTQEDF_x, &xNoTQEDF_data, &xNoTQEDF_x_size);
	BSFfast_importData_1d("BSFfast_DataCSV/dQCD-F.csv", &xNoTQEDF_m, &xNoTQEDF_x, &xNoTQEDF_data, &xNoTQCDF_x_size);
	
	
	printf("BSFfast > imported table sizes (scalars):\n          topP= %i x %i,topC= %i x %i, botP= %i x %i, botC= %i x %i, notP= %i x %i, notC= %i x %i |\n          QED= %i , QEDnot= %i , QCD= %i\n", 
					topPlatS_m_size, topPlatS_x_size, topCutS_m_size, topCutS_x_size, 
					botPlatS_m_size, botPlatS_x_size, botCutS_m_size, botCutS_x_size, 
					noTPlatS_m_size, noTPlatS_x_size, noTCutS_m_size, noTCutS_x_size, 
					xQEDS_x_size, xNoTQEDS_x_size, xNoTQCDS_x_size );
	printf("                               (fermion):\n          topP= %i x %i,topC= %i x %i, botP= %i x %i, botC= %i x %i, notP= %i x %i, notC= %i x %i |\n          QED= %i , QEDnot= %i , QCD= %i\n\n", 
					topPlatF_m_size, topPlatF_x_size, topCutF_m_size, topCutF_x_size, 
					botPlatF_m_size, botPlatF_x_size, botCutF_m_size, botCutF_x_size, 
					noTPlatF_m_size, noTPlatF_x_size, noTCutF_m_size, noTCutF_x_size, 
					xQEDF_x_size, xNoTQEDF_x_size, xNoTQCDF_x_size );
					
}









