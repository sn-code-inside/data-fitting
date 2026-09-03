/*****************************************************************
 *
 * File....: functions_NIST.h
 * Function: proto typing for functions_NIST.c
 * Author..: Tilo Strutz
 * Date....: 23.09.2009
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 *****************************************************************/

#ifndef FUNCT_NIST_H
#define FUNCT_NIST_H

/* certified parameters*/

double b_Filip[] = {
	-1467.48961422980,
	-2772.17959193342,
	-2316.37108160893,
	-1127.97394098372,
	-354.478233703349,
	-75.1242017393757,
	-10.8753180355343,
	-1.06221498588947,
	-0.670191154593408E-01,
	-0.246781078275479E-02,
	-0.402962525080404E-04
};

double b_Longley[] = {
	-3482258.63459582,
	15.0618722713733,
	-0.358191792925910E-01,
	-2.02022980381683,
	-1.03322686717359,
	-0.511041056535807E-01,
	1829.15146461355
};

double b_Wampler1[] = {
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000
};

double b_Wampler2[] = {
		 1.00000000000000,
		 0.100000000000000,
		 0.100000000000000E-01,
		 0.100000000000000E-02,
		 0.100000000000000E-03,
		 0.100000000000000E-04
};

double b_Wampler3[] = {
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000
};

double b_Wampler4[] = {
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000
};

double b_Wampler5[] = {
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000,
		 1.00000000000000
};

double b_Bennett5[] = {
	-2.5235058043E+03,
	4.6736564644E+01,
	9.3218483193E-01
};

double b_Eckerle4[] = {
 1.5543827178E+00,
 4.0888321754E+00,
 4.5154121844E+02
};

double b_Rat43[] = {
	6.9964151270E+02,
	5.2771253025E+00,
	7.5962938329E-01,
	1.2792483859E+00
};

double b_Rat42[] = {
 7.2462237576E+01,
 2.6180768402E+00,
 6.7359200066E-02
};

double b_BoxBOD[] = {
 2.1380940889E+02,
 5.4723748542E-01 
};

double b_Thurber[] = {
1.2881396800E+03,
 1.4910792535E+03,
 5.8323836877E+02,
 7.5416644291E+01,
 9.6629502864E-01,
 3.9797285797E-01,
 4.9727297349E-02
};

double b_MGH09[] = {
	1.9280693458E-01,
	1.9128232873E-01,
	1.2305650693E-01,
	1.3606233068E-01
};

double b_MGH10[] = {
 5.6096364710E-03,
 6.1813463463E+03,
 3.4522363462E+02
};

double b_Example[] = {
 2.25
};

int getNumEqualDigits( double* a, double* b, int M );

/* linear functions not required */

/* nonlinear functions */

double fNIST_BoxBOD( unsigned int i, double *cond, double *a);
double fNIST_BoxBOD_deriv( double (*funct)(unsigned int, double*, double*), unsigned int i, unsigned int j, unsigned int M,
	double* cond, double* a );
int init_NIST_BoxBOD( unsigned int N, double *obs,
						double *cond, double *a, unsigned char *a_flag, FILE *out);
double fNIST_MGH09( unsigned int i, double *cond, double *a);
double fNIST_MGH09_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fNIST_MGH09_deriv2( double (*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, double* cond, double* a );
int init_NIST_MGH09( unsigned int N, double *obs, double *cond,
					 double *a, unsigned char *a_flag, FILE *logfile);

double fNIST_thurber( unsigned int i, double *cond, double *a);
double fNIST_thurber_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fNIST_thurber_deriv2( double(*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, double* cond, double* a );
int init_NIST_thurber( unsigned int N, double *obs,
						double *cond, double *a, unsigned char *a_flag, FILE *out);

double fNIST_Rat42( unsigned int i, double *cond, double *a);
double fNIST_Rat42_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_NIST_Rat42( unsigned int N, double *obs,
						double *cond, double *a, unsigned char *a_flag, FILE *out);

double fNIST_Rat43( unsigned int i, double *cond, double *a);
double fNIST_Rat43_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_NIST_Rat43( unsigned int N, double *obs,
						double *cond, double *a, unsigned char *a_flag, FILE *out);

double fNIST_Eckerle4( unsigned int i, double *cond, double *a);
double fNIST_Eckerle4_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
int init_NIST_Eckerle4( unsigned int N, double *obs,
						double *cond, double *a, unsigned char *a_flag, FILE *out);

double fNIST_MGH10( unsigned int i, double *cond, double *a);
int init_NIST_MGH10( unsigned int N, double *obs,
						double *cond, double *a, unsigned char *a_flag, FILE *out);
double fNIST_MGH10_deriv( double (*funct)(unsigned int,double*,double*), 
					 unsigned int i, unsigned int j, unsigned int M, double *cond, double *a);
double fNIST_MGH10_deriv2( double (*funct)(unsigned int,double*,double*), 
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, double *cond, double *a);

double fNIST_Bennett5( unsigned int i, double* cond, double* a );
int init_NIST_Bennett5( unsigned int N, double *obs,
						double *cond, double *a, unsigned char *a_flag, FILE *out);
double fNIST_Bennett5_deriv( double (*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a );
double fNIST_Bennett5_deriv2( double (*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int k, unsigned int M, double* cond, double* a );

double fNIST_Bennett5LOG( unsigned int i, double* cond, double* a );
double fNIST_Bennett5LOG_deriv( double (*funct)(unsigned int, double*, double*),
	unsigned int i, unsigned int j, unsigned int M, double* cond, double* a );
int
init_NIST_Bennett5LOG( unsigned int N, double* obs,
	double* cond, double* a, unsigned char* a_flag, FILE* out );

#endif
