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


double **topPlat_data = NULL;
double **topCut_data = NULL;
double **botPlat_data = NULL;
double **botCut_data = NULL;
double **noTCut_data = NULL;
double **noTPlat_data = NULL;

double *xQED_data = NULL;
double *xNoTQED_data = NULL;
double *xNoTQCD_data = NULL;

double *topPlat_m = NULL;
double *botPlat_m = NULL;
double *noTPlat_m = NULL;
double *topCut_m = NULL;
double *botCut_m = NULL;
double *noTCut_m = NULL;

double *topPlat_x = NULL;
double *botPlat_x = NULL;
double *noTPlat_x = NULL;
double *topCut_x = NULL;
double *botCut_x = NULL;
double *noTCut_x = NULL;

double xQED_m;
double xNoTQED_m;
double xNoTQCD_m; 

double *xQED_x = NULL;
double *xNoTQED_x = NULL;
double *xNoTQCD_x = NULL;

int topPlat_x_size = 0;
int topPlat_m_size = 0;
int topCut_x_size = 0;
int topCut_m_size = 0;
int botPlat_x_size = 0;
int botPlat_m_size = 0;
int botCut_x_size = 0;
int botCut_m_size = 0;
int noTPlat_x_size = 0;
int noTPlat_m_size = 0;
int noTCut_x_size = 0;
int noTCut_m_size = 0;

int xQED_x_size = 0;
int xNoTQED_x_size = 0;
int xNoTQCD_x_size = 0;

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
	//TOP Pleteau
double BSFfast_sigveff_QCD_SU_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	//not handling extrapolation for QCD models.
	if (logm>topPlat_m[topPlat_m_size-1] || logm<topPlat_m[0] || logx>topPlat_x[topPlat_x_size-1] || logx<topPlat_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = exp( BSFfast_bilinear_interpolate(topPlat_m, topPlat_x,topPlat_data,log(m),log(x), topPlat_m_size, topPlat_x_size) );
    return out;
    }
	//TOP Cutoff
double BSFfast_sigveff_QCD_SU_cut(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>topCut_m[topCut_m_size-1] || logm<topCut_m[0] || logx>topCut_x[topCut_x_size-1] || logx<topCut_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(topCut_m, topCut_x,topCut_data,log(m),log(x), topCut_m_size, topCut_x_size) ));
    return out;
    }

	//BOTTOM Plateau
double BSFfast_sigveff_QCD_SD_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>botPlat_m[botPlat_m_size-1] || logm<botPlat_m[0] || logx>botPlat_x[botPlat_x_size-1] || logx<botPlat_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
    double out = (exp( BSFfast_bilinear_interpolate(botPlat_m, botPlat_x,botPlat_data,log(m),log(x), botPlat_m_size, botPlat_x_size) ));
    return out;
    }
	//BOTTOM Plateau
double BSFfast_sigveff_QCD_SD_cut(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>botCut_m[botCut_m_size-1] || logm<botCut_m[0] || logx>botCut_x[botCut_x_size-1] || logx<botCut_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(botCut_m, botCut_x,botCut_data,log(m),log(x), botCut_m_size, botCut_x_size) ));
    return out;
    }

	//no-Trans Plateau
double BSFfast_sigveff_QCD_S_plat(double m, double x) {
    double logm = log(m), logx=log(x);
	if (logm>noTPlat_m[noTPlat_m_size-1] || logm<noTPlat_m[0] || logx>noTPlat_x[noTPlat_x_size-1] || logx<noTPlat_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(noTPlat_m, noTPlat_x, noTPlat_data,log(m),log(x), noTPlat_m_size, noTPlat_x_size) ));
    return out;}
	//no-Trans Plateau
double BSFfast_sigveff_QCD_S_cut(double m, double x) {
	double logm = log(m), logx=log(x);
	if (logm>noTCut_m[noTCut_m_size-1] || logm<noTCut_m[0] || logx>noTCut_x[noTCut_x_size-1] || logx<noTCut_x[0] ){
		printf("OOB: logm=%.3f, logx=%.5f\n", logm, logx); 
		return -1.;
	}
	double out = (exp( BSFfast_bilinear_interpolate(noTCut_m, noTCut_x,noTCut_data,log(m),log(x), noTCut_m_size, noTCut_x_size) ));
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
	
    double prefactor = xQED_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp( BSFfast_linear_interpolate( xQED_x, xQED_data, logx , xQED_x_size)  );
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
	 
    double prefactor = xNoTQED_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp( BSFfast_linear_interpolate(xNoTQED_x, xNoTQED_data, logx , xNoTQED_x_size) 	);
    return prefactor;
}
	// xQCD-NoTrans-const
double BSFfast_sigveff_dQCD_const(double alpha, double m, double x) {
	double r = alpha / sigBSFeff_alpha0; 
	double logx = log(x) + 2 * log(r);
	/*if ( logx>xNoTQCD_x[ xNoTQCD_x_size - 1 ] || logx<xNoTQCD_x[0] ){
		printf("OOB: %.5f > logx=%.5f > %.5f\n", xNoTQCD_x[xNoTQCD_x_size-1], logx, xNoTQCD_x[1]); 
		return -1.;
	}*/
		
    double prefactor = xNoTQCD_m;
    prefactor = prefactor * prefactor * r * r / m / m; //  m0^2/m^2 * r^2
    prefactor = prefactor * exp(BSFfast_linear_interpolate(xNoTQCD_x, xNoTQCD_data, logx , xNoTQCD_x_size) );
    return prefactor;
}


__attribute__((constructor))
void init_library(void) {
	// import all 9 data files into the global variables
		//function calls are still lacking the data-length information to be passed on !
	BSFfast_importData_2d("BSFfast_DataCSV/sigBSFeff_stop_asPlateau1.RD.csv", &topPlat_m, &topPlat_x, &topPlat_data, &topPlat_m_size, &topPlat_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/sigBSFeff_stop_asCutoff1.RD.csv", &topCut_m, &topCut_x, &topCut_data, &topCut_m_size, &topCut_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/sigBSFeff_sbottom_asPlateau1.RD.csv", &botPlat_m, &botPlat_x, &botPlat_data, &botPlat_m_size, &botPlat_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/sigBSFeff_sbottom_asCutoff1.RD.csv", &botCut_m, &botCut_x, &botCut_data, &botCut_m_size, &botCut_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/sigBSFeff-NoTrans_sbottom_asPlateau1.RD.csv", &noTPlat_m, &noTPlat_x, &noTPlat_data, &noTPlat_m_size, &noTPlat_x_size);
	BSFfast_importData_2d("BSFfast_DataCSV/sigBSFeff-NoTrans_sbottom_asCutoff1.RD.csv", &noTCut_m, &noTCut_x, &noTCut_data, &noTCut_m_size, &noTCut_x_size);

	BSFfast_importData_1d("BSFfast_DataCSV/xgridSigBSFeff_QED_asConst0.1_m1GeV.csv", &xQED_m, &xQED_x, &xQED_data, &xQED_x_size);
	BSFfast_importData_1d("BSFfast_DataCSV/xgridSigBSFeff_QEDnoTrans_asConst0.1_m1GeV.csv", &xNoTQED_m, &xNoTQED_x, &xNoTQED_data, &xNoTQED_x_size);
	BSFfast_importData_1d("BSFfast_DataCSV/xgridSigBSFeff_QCDnoTrans_asConst0.1_m1GeV.csv", &xNoTQED_m, &xNoTQED_x, &xNoTQED_data, &xNoTQCD_x_size);
	
	printf("BSFfast > imported table sizes:\n          topP= %i x %i,topC= %i x %i, botP= %i x %i, botC= %i x %i, notP= %i x %i, notC= %i x %i |\n          QED= %i , QEDnot= %i , QCD= %i\n\n", 
					topPlat_m_size, topPlat_x_size, topCut_m_size, topCut_x_size, 
					botPlat_m_size, botPlat_x_size, botCut_m_size, botCut_x_size, 
					noTPlat_m_size, noTPlat_x_size, noTCut_m_size, noTCut_x_size, 
					xQED_x_size, xNoTQED_x_size, xNoTQCD_x_size );
					
}









