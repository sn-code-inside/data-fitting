/*****************************************************************
 *
 * File....:	fitting.c
 * Function:	data fitting with least squares
 * Author..:	Tilo Strutz
 * Date....:	28.09.2009
 *
 * changes:
 * 20.08.2012 implementation of RANSAC, M-score
 * 29.04.2013 bugfix MAX_CONDITIONS vs M_MAX
 * 28.01.2014 new option cw
 * 09.12.2014 output of weighting and outlier-detection mode
 * 25.02.2021 new: f(x|a) = (a1*x)^a2 + a3 and f(x|a) = a2 * log(a1*x)
 * 05.11.2025 new: variable model function via option string (-y)
 * 05.01.2026 mu_fac based on opportunities for parameter initialistion
 * 11.05.2026 new hyperbel function
 * 19.06.2026 new option '-D'
 * 19.07.2026 scaling of conditions improved, extended; new: NIST_BENNETT5_LOG
 * 12.08.2026 new EXAMPLE
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#ifndef WIN32
#include <unistd.h>
#else
#include <Windows.h> // sysinfoapi.h> /* for GetSystemInfo() */
#endif
#include "get_option.h"
#include "errmsg.h"
#include "matrixUtils.h"
#include "functions.h"
#include "functionsNIST.h"
#include "macros.h"
#include "ls.h"
#include "prototypes.h"
#include "defines.h"
#include "tinyexpr.h"
#include "cpuInfo.h"

#define MAX_OUTPUT_LEN 1024

double x; /* single observation for tinyexpr */
double fourthRootDBLEPS, thirdRootDBLEPS;

/* #define OUTPUT_DEVIATES */
#define COMPARISON_WITH_CERTIFIED_VALUES
#undef COMPARISON_WITH_CERTIFIED_VALUES

/*---------------------------------------------------------------
 *	getNumEqualDigits()
 *--------------------------------------------------------------*/
int getNumEqualDigits( double* a, double* b, int M )
{
	/* number of decimal digits in double defined in float.h */
	int numOfEqualDigits = 0, minOfEqualDigits = OP_DBL_Digs;
	for (int j = 0; j < M; j++) /* check all parameters*/
	{
		int exponent_a, exponent_b, exponent_max;
		double diff;

		numOfEqualDigits = 0;
		if (a[j] == 0.0) exponent_a = 0;
		else exponent_a = (int)floor( log10( fabs( a[j] ) ) );
		if (b[j] == 0.0) exponent_b = 0;
		else exponent_b = (int)floor( log10( fabs( b[j] ) ) );
		exponent_max = max( exponent_a, exponent_b );

		diff = fabs( a[j] - b[j] ) / pow( 10, exponent_max );
		if (diff == 0.0) numOfEqualDigits = OP_DBL_Digs;
		else
		{
			double f = 1.;
			int digits;
			for (digits = 0; digits < OP_DBL_Digs; digits++)
			{
				if (round( diff * f ) != 0.0) 		break;
				f *= 10;
			}
			numOfEqualDigits = digits;
		}
		/* pick smallest number of equal digits */
		if (minOfEqualDigits > numOfEqualDigits) minOfEqualDigits = numOfEqualDigits;
	}
	return minOfEqualDigits;
}


/*---------------------------------------------------------------
 *	renumberVariables()
 * expressions provided by the use could have arbitrary parameter
 * varaiables a4, a2, a10; these are renumbered to a1, a2, a3, ...
 *--------------------------------------------------------------*/
static int renumberVariables( const char* input, char* output, FILE* out )
{
	int mapping[M_MAX]; /* Index 1..11 */
	int nextIndex = 1;
	int modidifiction_flag = 0;
	size_t outPos = 0;
	size_t i;

	for (i = 0; i < M_MAX; i++)
	{
		mapping[i] = 0;
	}

	i = 0;
	while (input[i] != '\0')
	{
		/* Check for a1 ... a11 */
		if (input[i] == 'a' && isdigit( (unsigned char)input[i + 1] ))
		{
			int number = 0;
			size_t j = i + 1;

			while (isdigit( (unsigned char)input[j] ))
			{
				number = number * 10 + (input[j] - '0');
				j++;
			}

			/* Only accept a1..a11 */
			if (number < M_MAX)
			{
				if (mapping[number] == 0)
				{
					if (number != nextIndex)
						modidifiction_flag = 1;
					mapping[number] = nextIndex++;
				}

				outPos += sprintf( output + outPos,
					"a%d",
					mapping[number] );

				i = j;
				continue;
			}
			else
			{
				fprintf( stderr, "\n ! Wrong parameter number: %d !", number );
				fprintf( out, "\n# ! Wrong parameter number: %d !", number );
			}
		}
		output[outPos++] = input[i++];
	}

	output[outPos] = '\0';
	return modidifiction_flag;
}

/*---------------------------------------------------------------
 *	createFinalFunction()
 *--------------------------------------------------------------*/
static void createFinalFunction( const char* input, char* output,
	double* a, FILE* out )
{
	size_t outPos = 0;
	size_t i;

	i = 0;
	while (input[i] != '\0')
	{
		/* Check for a1 ... a11 */
		if (input[i] == 'a' && isdigit( (unsigned char)input[i + 1] ))
		{
			int number = 0;
			size_t j = i + 1;

			while (isdigit( (unsigned char)input[j] ))
			{
				number = number * 10 + (input[j] - '0');
				j++;
			}

			/* Only accept a1..a11 */
			if (number > 0 && number < M_MAX)
			{
				outPos += sprintf( output + outPos,	" %.4e ", a[number - 1] );
				i = j;
				continue;
			}
			else
			{
				fprintf( stderr, "\n ! Wrong parameter number: %d !", number );
				fprintf( out, "\n# ! Wrong parameter number: %d !", number );
			}
		}

		output[outPos++] = input[i++];
	}

	output[outPos] = '\0';
}

/*---------------------------------------------------------------
 *	main()
 *--------------------------------------------------------------*/
int
main( int argc, char* argv[] )
{
	char* rtn = "main";
	char* field;	/* used for reading text files	*/
	char* inname = NULL;	/* filename of input data */
	char* outname = NULL;	/* filename of results		*/
	/* string with list of columns containing conditions */
	char* column_cond_str = NULL;
	char funcStringInput[511]; /* for model function given as string */
	char funcString[511]; /* corrected string */
	char formatString[20]; /* for final output in text file */
	char formatStringP[20]; /* for final output in text file */
	char formatStringText[20]; /* for final output in text file */
	/* used for reading text files	*/
	char line[MAXLINELENGTH + 1], * ptr;
	/* pointer to model function */
	double (*funct) (unsigned int, double*, double*) = NULL;
	/* pointer to its derivative */
	double (*funct_deriv) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, double*, double*) = NULL;
	double (*funct_deriv2) (double(*)(unsigned int, double*, double*),
		unsigned int, unsigned int, unsigned int, unsigned int, double*, double*) = NULL;
	/* pointer to initialisation function */
	int (*init) (unsigned int, double*, double*, double*,
		unsigned char*, FILE*) = NULL;
	int err = 0, res;
	unsigned int N;							  /* number of observations */
	unsigned int M = 1;					  /* number of model parameters */
	int cntInliers;		            /* counter for observations	*/
	int columnCond[MAX_CONDITIONS], col, ch;
	int columnUncertainties = 0;  /* column containing the weights */
	int columnObs = 0;       /* column containing the observations */
	unsigned int condDim = 1;		/* dimensionality of conditions */
	unsigned int obsDim = 1;		/* dimensionality of observations */
	int type = UNSPECIFIED;	/* type of model function */
	int M_flag = 0;				  /* flag for model LINEAR, POLYNOMIAL_REG	*/
	int numerical_flag = 0;	/* force numerical derivation */
	int scaling_flag = 0;		/* enable scaling of conditions */
	int forget_flag = 0;	  /* enables reset of weights after
												   * outlier removal	*/
	int numOfOutliers = 0;	/* number of outliers	*/
	int outMode = 0;		/*
											 * 0 .. no outlier removal
											 * 1 .. enable z-score + Chauvenet's
											 * 2 .. enable ClubOD
											 * 3 .. enable M-score + Chauvenet's
											 * 4 .. enable RANSAC
											 */
	int weightMode = 0; /*
											 * 0 .. use equal weights (no weighting)
											 * 1 .. enable deviates based weighting
											 * 2 .. weighting by binning
											 * 3 .. uncertainties are given along with data
											 */
	unsigned int obsPerBin = 50; /* observations per bin for
																* weightMode = 2
																*/
	int algoMode = 1;	  /*
											 * 0 .. use simple matrix inversion, M <= 5
											 * 1 .. singular value decomposition
											 * 2 .. LU decomposition
											 */
	int iterations = 0;
	int iterWmax = 120;	/* max number of iterations for weights
													 estimation */
	int iterStop_flag;	/* stop weights iteration when convergence
												 is reached */
	int iterFinal_flag;	/* last iteration after convergence is
												 reached */
	int numOutliers = 0;			/* number of detected outliers */
	int outDetect_flag;	/* indicates, whether outlier detection
												 * has been performed
												 */
												 /* for formatstring "%12.6f und %*e" */
	int digitsTotal = 12, digitsDecimal = 6, digitsIntermediate = 12;
	int argc_orig;
	double** jacob = NULL;	/* Jacobian matrix J		*/
	double** covar = NULL;	/* covariance matrix C	*/
	double* weights = NULL;	/* vector for weights		*/
	double* weightsOld = NULL;	/* vector for weights */
	double* observations = NULL;		/* observations */
	double* datac = NULL;	/* calculated values based on parameters
												 */
	double* conditions = NULL;	/* conditions X */
	double a[M_MAX];			/* parameter of model function		*/
	unsigned char a_flag[M_MAX];		/* corresponding flags	*/
	double* deviates = NULL;				/* = [observations - f(x|a) ]		*/
	double* deviatesAbsolute = NULL;	/* = abs[observations - f(x|a) ] */
	double chisq = 0., sum, variance, energy, mean, meanOfWeights;
	double chisqTargetValue = 0;
	double gfit = 0.;				/* goodness of fit	*/
	double uncertainty;	/* sigma_y, root of sample Variance, RMSE	*/
	double scalingFactor = 1.;	/* factor for scaling of conditions	*/
	LS_PARAMS lsq_params, * lsParams;
	FILE* in = NULL;
	FILE* out = stdout;

	te_expr* teExpr = NULL;	/* pointer to expression for tinyexpr */
	/* array of variables for tinyexpr */
	te_variable teVars[] = {
			{"x", &x},
			{"a1", &a[0]},
			{"a2", &a[1]},
			{"a3", &a[2]},
			{"a4", &a[3]},
			{"a5", &a[4]},
			{"a6", &a[5]},
			{"a7", &a[6]},
			{"a8", &a[7]},
			{"a9", &a[8]},
			{"a10", &a[9]},
			{"a11", &a[10]}
	};

	errno = 0; /* reset global error number, is checked in lsLinear.c	*/
	/* reset all flags, assume that initial parameters must
	 * be determined by the programm itself
	 */
	for (unsigned int j = 0; j < M_MAX; j++)
	{
		a_flag[j] = 0;
	}

	/* set pointer to object of flags and initialise LS parameters */
	lsParams = &lsq_params;
	lsParams->linear = 1;
	lsParams->svd = 1;		/* use SVD for linear systems as default	*/
	lsParams->LM = 1;		/* use Levenberg-Marquardt as default	*/
	lsParams->DM = 0;		/* use identity matrix for 
											 * Levenberg-Marquardt as default	*/
	lsParams->chisqTarget_flag = 0;
	lsParams->trueH = 0;
	lsParams->positiveParams = 0; /* for option -y limitation of 
																 * parameters domain of definition */
	lsParams->expression_flag = 0;
	lsParams->GNspecial_flag = 0;
	lsParams->modifyVanishingGradients_flag = 0;
	lsParams->uphillAcceptanceThreshold = 0.0; /* in combination 
	                                with -G and GNspecial_flag */
	lsParams->muFac = 0.001; /* factor for Levenberg-Marquardt */

	fourthRootDBLEPS = sqrt( sqrt( DBL_EPSILON ) );
	thirdRootDBLEPS = cbrt( DBL_EPSILON ); /* cubic root */

	ITERAT_MAX = 2000;	/* declared in lsNonlinear.c	*/
	funcStringInput[0] = 0; /* initialize empty string */
#ifdef TESTT
	{
		double** a, ** b;
		double	det;
		a = matrix( 5, 5 );	/*  matrix */
		b = matrix( 5, 5 );	/*  matrix */

		a[0][0] = 2.; a[0][1] = 2.; a[0][2] = 3.; a[0][3] = 4.; a[0][4] = 5.;
		a[1][0] = 2.; a[1][1] = 3.; a[1][2] = 5.; a[1][3] = 5.; a[1][4] = 5.;
		a[2][0] = 1.; a[2][1] = 4.; a[2][2] = 4.; a[2][3] = 4.; a[2][4] = 2.;
		a[3][0] = 1.; a[3][1] = 2.; a[3][2] = 1.; a[3][3] = 5.; a[3][4] = 3.;
		a[4][0] = 3.; a[4][1] = 3.; a[4][2] = 3.; a[4][3] = 2.; a[4][4] = 6.;

		det = inverse_5x5( a, b );
		fprintf( stderr, "\n det = %f\n", det );
		for (i = 0; i < 5; i++)
		{
			for (j = 0; j < 5; j++)
			{
				fprintf( stderr, " %6.2f", b[i][j] );
			}
			fprintf( stderr, "\n " );
		}

		free_matrix( &a );
		free_matrix( &b );
		exit( 1 );
	}
#endif
	argc_orig = argc; /* remember the number of arguments	*/
	/* check command-line parameters */
	while ((optstr =
		(char*)get_option( argc, (const char**)argv )) != NULL)
	{
		switch (optstr[1])
		{
			unsigned int j;
		case 'a':
			switch (optstr[2])
			{
			case '\0':
				algoMode = atoi( OptArg );
				/* 0 .. use simple matrix inversion, M <= 5
				 * 1 .. singular value decomposition
				 * 2 .. LU decomposition
				 */
				break;
			default:
				/* take number after 'a' as number of parameter
				 * since OPTIONSTRING (usage.c) does not contain 'a0'
				 * the domain of definition of j is limited to 1...9
				 */
				j = atoi( &(optstr[2]) );
				if (j < 1)
					fprintf( stderr, "\n!! parameter %s is invalid. Index must be in [1;9] !!\n", optstr );
				else
				{
					a[j - 1] = atof( OptArg );
					a_flag[j - 1] = 1;
				}
			}
			break;
		case 'b':
			obsPerBin = atoi( OptArg ); /* observations per bin */
			break;
		case 'c':
			switch (optstr[2])
			{
			case 'c':	/* string of comma-separated column numbers */
				column_cond_str = OptArg;
				break;
			case 'o':	/* column number of observations */
				columnObs = atoi( OptArg );
				break;
			case 'u':	/* column number of weights */
				columnUncertainties = atoi( OptArg );
				weightMode = 3;
				break;
			default: /* option '-c'	*/
				scaling_flag = 1;
				break;
			}
			break;
		case 'D':	/* use diagonal instead identity matrix for LM */
			lsParams->DM = 1;
			break;
		case 'e':	/* enable uphill walk and step damping for GaussNewton */
			lsParams->GNspecial_flag = 1;
			break;
		case 'f':
			forget_flag = 1;
			/* forget weights after outlier removal
			 */
			break;
		case 'G':	/* use plain Gauss-Newton */
			lsParams->LM = 0;
			break;
		case 'g':	/* allow gradient manipulation */
			lsParams->modifyVanishingGradients_flag = 1;
			break;
		case 'H':	/* use true Hessian matrix	*/
			lsParams->trueH = 1;
			break;
		case 'i':
			inname = OptArg;
			break;
		case 'o':
			outname = OptArg;
			break;
		case 'I':	/* maximum number of iterations	*/
			ITERAT_MAX = atoi( OptArg );
			break;
		case 'm':	/* model function	*/
			type = atoi( OptArg );
			break;
		case 'M':	/* number of parameters (type == LINEAR)	*/
			M = atoi( OptArg );
			M_flag = 1;
			break;
		case 'n':	/* force usage of numerical derivation	*/
			numerical_flag = 1;
			break;
		case 'p':	/* in conjunction with '-y' keep parameters positive, avoiding NAN or INF */
			lsParams->positiveParams = 1;
			break;
		case 's':
			lsParams->svd = 0;
			/* disable special SVD function for solving linear model
			 */
			break;
		case 't':
			chisqTargetValue = atof( OptArg );
			lsParams->chisqTarget_flag = 1;
			break;
		case 'w':
			weightMode = atoi( OptArg );
			/* 0 ... equal weights; 1 ... deviates based;
			 * 2 ... Bin-wise;  3 ... from input file
			 */
			break;
		case 'x':
			outMode = atoi( OptArg );
			/* 0 ... no removal;
			 * 1 ... z-score + Chauvenet's;
			 * 2 ... CluBOD;
			 * 3 ... M-score + Chauvenet's;
			 * 4 ... RANSAC;
			 */
			break;
		case 'y':	/* string for model function	*/
			strcpy_s( funcStringInput, 511, OptArg ); /* maximal string length = 511 characters */
			lsParams->expression_flag = 1;
			break;
		case 'z':	/* string for output formatting	*/
			switch (optstr[2])
			{
			case '0':	/* format string: number of digits for intermediate output */
				digitsIntermediate = atoi( OptArg );
				break;
			case '1':	/* format string: total number of digits for final output */
				digitsTotal = atoi( OptArg );
				break;
			case '2':	/* format string: number of decimal digits */
				digitsDecimal = atoi( OptArg );
				break;
			default: /* option '-z'	*/
				usage( argv[0] ); /* provides  help */
				err = errmsg( ERR_CALL, rtn, "option -z?", 0 ); ;
				goto endfunc;
			}
			break;
		case '?':
		default:
			usage( argv[0] ); /* provides  help */
			err = errmsg( ERR_CALL, rtn, "unknown", 0 ); ;
			goto endfunc;
		}
	}
	/* check, whether all mandatory options were given */
	err = check_opt( argv[0] );
	fprintf( stderr, "\n command:" );
	for (unsigned int i = 0; i < (unsigned)argc_orig; i++)
	{
		fprintf( stderr, " %s", argv[i] );
	}
	if (err)
	{
		fprintf( stderr, "\n" );
		goto endfunc;
	}

	/*
	 * evaluation of programm options
	 */
	if (outname == NULL)
	{
		fprintf( stderr, "\n Name of output file missing!" );
		usage( argv[0] );
		goto endfunc;
	}
	if (inname == NULL)
	{
		fprintf( stderr, "\n Name of input file missing!" );
		usage( argv[0] );
		goto endfunc;
	}
	if ((type == LINEAR || type == LINEAR_2 || type == POLYNOMIAL ||
		type == POLYNOMIAL_REG) && M_flag == 0)
	{
		fprintf( stderr, "\n You have chosen mode '-m %d', ", type );
		fprintf( stderr,
			"but forgotten to set the function order '-M'" );
		usage( argv[0] );
		goto endfunc;
	}
	if (lsParams->LM)
	{
		if (lsParams->GNspecial_flag)
		{
			fprintf( stderr, "\n You have not chosen mode `-G'" );
			fprintf( stderr, "\n ==> option `-e' remains unused" );
			lsParams->GNspecial_flag = 0;
		}
	}
	if (lsParams->expression_flag)
	{
		fprintf( stderr, "\n You have chosen mode `-y'" );
		fprintf( stderr, " ==> option `-m' will be ignored (even it must be provided)" );
	}

	if (digitsTotal > 12) digitsTotal = 12;
	if (digitsTotal < 6) digitsTotal = 6;
	if (digitsDecimal > 6) digitsTotal = 6;
	if (digitsDecimal > digitsTotal - 2) digitsDecimal = digitsTotal - 2;
	if (digitsDecimal < 2) digitsTotal = 2;
	sprintf( formatString, "%%%d.%df ", digitsTotal, digitsDecimal ); /* set formatting for results table */
	sprintf( formatStringP, "%%+%d.%df ", digitsTotal, digitsDecimal ); /* set output formatting */
	sprintf( formatStringText, "%%%ds ", digitsTotal - 2 ); /* set formatting for results table*/


	/* initialise default columns of conditions	 */
	for (unsigned int i = 0; i < MAX_CONDITIONS; i++)
	{
		/* conditions in increasing order */
		columnCond[i] = i + 1;
	}

	/* if comma-separated list of columns is given */
	if (column_cond_str != NULL)
	{
		unsigned int i;
		/* convert column string into numbers */
		i = col = 0;
		do
		{
			/* loop until all columns are read or
			 * maximal number of columns is reached
			 */
			ptr = &(column_cond_str[i]);
			res = sscanf( ptr, "%d", &(columnCond[col]) );
			do
			{ /* go to next number */
				i++;
				ch = column_cond_str[i];
			} while (ch != '\0' && ch != ',');
			i++;
			col++;
		} while (ch != '\0' && col < MAX_CONDITIONS);
		for (i = col; i < MAX_CONDITIONS; i++)
		{
			columnCond[i] = columnCond[i - 1] + 1;
		}
	}

	/*
	 * open the input file
	 * determine the number of data sets
	 */
	in = fopen( inname, "rt" );
	if (in == NULL)
	{
		err = errmsg( ERR_OPEN_READ, rtn, inname, 0 );
		goto endfunc;
	}
	/* open out file */
	out = fopen( outname, "wt" );
	if (out == NULL)
	{
		err = errmsg( ERR_OPEN_WRITE, rtn, outname, 0 );
		goto endfunc;
	}

	fprintf( out, "# %s ===============================", rtn );
	fprintf( out, "\n# use data file: %s", inname );

	/* determine number of observations by counting of valid lines */
	N = 0;
	while ((ptr = fgets( line, MAXLINELENGTH, in )) != NULL)
	{
		/* skip comment lines (starting with '#') and empty ones */
		if (is_data_line( line, MAXLINELENGTH ))
		{
			N++;
			if (strlen( line ) == MAXLINELENGTH - 1)
			{
				fprintf( stderr,
					"\n lines of input file are too long (>%d)",
					MAXLINELENGTH );
				fprintf( stderr, ", increase MAXLINELENGTH" );
			}
		}
	}
	fclose( in );

	fprintf( stderr, "\n datafile contains %d data points\n", N );

	/*
	 * set number of parameters and redirect pointer to functions
	 */
	fprintf( out, "\n# -m %d ", type );
	switch (type)
	{
	case CONSTANT:
		/* y = a1 */
		fprintf( out, ", constant function y = a1" );
		fprintf( stderr, "\n constant function y = a1" );
		funct_deriv = fconstant_deriv;
		M = 1;
		break;

	case LINEAR:
		/* f(x|a) = a1 + Sum_j(a_j*x_j)  */
		fprintf( out, ", a1 + Sum_j(a_j*x_j), order %d", M - 1 );
		fprintf( stderr, "\n a1 + Sum_j(a_j*x_j), order %d", M - 1 );
		funct_deriv = flin_deriv;
		/* M is set via program parameter */
		condDim = M - 1; /* first parameter a1 is just an offset
											 * w/o corresponding condition
											 */
		break;
	case LINEAR_2:
		/* f(x|a) = Sum_j(a_j*x_j)  */
		fprintf( out, ", Sum_j(a_j*x_j), order %d", M );
		fprintf( stderr, "\n Sum_j(a_j*x_j), order %d", M );
		funct_deriv = flin2_deriv;
		/* M is set via program parameter */
		condDim = M; /* first parameter a1 is not used	 */
		break;
	case POLYNOM_2NDORD:
		/*
		 * f(x|a) = a1 + a2 * x + a3 * x^2
		 */
		fprintf( out, ", a1 + a2 * x + a3 * x^2" );
		fprintf( stderr, "\n a1 + a2 * x + a3 * x^2" );
		funct_deriv = fpolynom2_deriv;
		M = 3;
		break;
	case POLYNOM_3RDORD:
		/*
		 * f(x|a) = a1 + a2 * x + a3 * x^2 + a4 * x^3
		 */
		fprintf( out, ", a1 + a2 * x + a3 * x^2 + a4 * x^3" );
		fprintf( stderr, "\n a1 + a2 * x + a3 * x^2 + a4 * x^3" );
		funct_deriv = fpolynom3_deriv;
		M = 4;
		break;
	case POLYNOMIAL:
		/*
		 * f(x|a) = a1 + a2 * x + a3 * x^2 + ...
		 */
		fprintf( out, ", polynomial of %dth order", M - 1 );
		fprintf( stderr, "\n polynomial of %dth order", M - 1 );
		funct_deriv = fpolynomial_deriv;
		/* M is set via program parameter */
		break;
	case POLYNOMIAL_REG:
		/*
		 * f(x|a) = a1 + a2 * x + a3 * x^2 + ...
		 */
		if (M == 2)
		{
			fprintf( out, ", polynomial: a1 + a2 * x" );
			fprintf( stderr, "\n polynomial: a1 + a2 * x" );
		}
		else if (M == 3)
		{
			fprintf( out, ", polynomial a1 + a2 * x + a3 * x^2" );
			fprintf( stderr, "\n polynomial a1 + a2 * x + a3 * x^2" );
		}
		else if (M == 4)
		{
			fprintf( out, ", polynomial a1 + a2 * x + a3 * x^2 + a4 * x^3" );
			fprintf( stderr, "\n polynomial a1 + a2 * x + a3 * x^2 + a4 * x^3" );
		}
		else
		{
			fprintf( out, ", polynomial of %dth order", M - 1 );
			fprintf( stderr, "\n polynomial of %dth order", M - 1 );
		}
		fprintf( out, ", regularised (nonlinear)" );
		fprintf( stderr, ", regularised (nonlinear)" );

		funct = fpolynomial;
		init = init_polynomial;
		funct_deriv = fpolynomial_deriv;
		/* M is set via program parameter */
		lsParams->linear = 0;	/* nonlinear */
		break;
	case COSINE_LIN:
		/*
		 * f(x|b) = b1 + b2 * cos( x - b3)
		 * f(x|a) = a1 + a2 * cos( x) + a3 * sin( x)
		 * a2 = b2 * cos(b3), a3 = b2 * sin(b3)
		 * a1 = b1
		 */
		fprintf( out, ", a1 + a2 * cos( x) + a3 * sin( x) " );
		fprintf( stderr, "\n# a1 + a2 * cos( x) + a3 * sin( x) " );
		funct_deriv = fcosineLin_deriv;
		M = 3;
		break;
	case COSINE_NONLIN:
		/*
		 * f(x|b) = b1 + b2 * cos( x - b3)
		 */
		fprintf( out, ", a1 + a2 * cos( x- a3) " );
		fprintf( stderr, "\n# a1 + a2 * cos( x- a3) " );
		funct = fcosineNonLin;
		init = init_cosineNonLin;
		funct_deriv = fcosineNonLin_deriv;
		M = 3;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case COSINE:
		/*
		 * f(x|a) = a1 + a2 * cos( a3*x + a4)
		 */
		fprintf( out, ", a1 + a2 * cos( a3*x + a4)" );
		fprintf( stderr, "\n a1 + a2 * cos( a3*x + a4)" );
		funct = fcosine;
		init = init_cosine;
		funct_deriv = fcosine_deriv;
		M = 4;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case COSINE2:
		/*
		 * f(x|a) = a1 + a2*cos(a3*x-a4) + a5*cos(2*a3*x-a6)
		 */
		fprintf( out, ", a1 + a2*cos(a3*x-a4) + a5*cos(2*a3*x-a6)" );
		fprintf( stderr, "\n a1 + a2*cos(a3*x-a4) + a5*cos(2*a3*x-a6)" );
		funct = ftrigonometric2;
		init = init_cosine2;
		funct_deriv = f_deriv;	/* use numerical differentiation */
		numerical_flag = 1; /* no derivation function implemented yet */
		M = 6;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case COSINE_TREND:
		/*
		 * f(x|a) = a1 + a2 * x + a3 * cos( x - a4)
		 */
		fprintf( out, "# a1 + a2 * x + a3 * cos( x - a4)" );
		printf( "\n  a1 + a2 * x + a3 * cos( x - a4)" );
		funct = fcosine_trend;
		init = init_cosine_trend;
		funct_deriv = fcosine_trend_deriv;
		M = 4;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case EXPONENTIAL:
		/* f(x|a) = a1 + a2 * exp( a3 * x) */
		fprintf( out, ", a1 + a2 * exp( a3 * x)" );
		fprintf( stderr, "\n a1 + a2 * exp( a3 * x)" );
		funct_deriv = fexponential_deriv;
		funct = fexponential;
		if (fabs( a[0] - 99999 ) < TOL && fabs( a[1] - 99999 ) < TOL
			&& fabs( a[2] - 99999 ) < TOL)
		{ /* "-a1 99999 -a2 99999 -a3 99999" triggers linear regression */
			init = init_exponentialLinearised; /* use linear regression */
		}
		else if (fabs( a[0] - 11111 ) < TOL && fabs( a[1] - 11111 ) < TOL
			&& fabs( a[2] - 11111 ) < TOL)
		{/* use dummy values based on curvature */
			init = init_exponentialdummy;
		}
		else 	init = init_exponential;/* normal initialisation */
		M = 3;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case EXPONENTIAL2:
		/* f(x|a) = a2 * exp( a3 * x) */
		fprintf( out, ", a1 * exp( a2 * x)" );
		fprintf( stderr, "\n a1 * exp( a2 * x)" );
		funct = fexpon2;
		init = init_expon2;
		funct_deriv = fexpon2_deriv;
		M = 2;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case EXPONENTIAL2_LIN:
		/* ln(f(x|a)) = ln(a2) + a3 * x */
		fprintf( out, ", ln(f(x|a)) = ln(a1) + a2 * x" );
		fprintf( stderr, "\n ln(f(x|a)) = ln(a1) + a2 * x" );
		funct_deriv = flin_deriv;
		M = 2;
		break;
	case EXP_BASE_X:
		/* f(x|a) = (a1*x)^a2 + a3 */
		fprintf( out, ", EXP_BASE_X: (a1*x)^a2 + a3" );
		fprintf( stderr, "\n EXP_BASE_X: (a1*x)^a2 + a3" );
		funct = fexpBaseX;
		init = init_fexpBaseX;
		funct_deriv = fexpBaseX_deriv;
		M = 3;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case LOGARITHM:
		/* f(x|a) = a1 + a2 * log( a3 + x) */
		fprintf( out, ", a1 + a2 * log( a3 + x)" );
		fprintf( stderr, "\n a1 + a2 * log( a3 + x)" );
		funct = flogarithmic;
		init = init_logarithmic;
		funct_deriv = flogarithmic_deriv;	/*  */
		M = 3;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case HYPERBEL:
		/* f(x|a) = a1 / (x+a2) */
		fprintf( out, ", a1 / (x+a2)" );
		fprintf( stderr, "\n a1 / (x+a2)" );
		funct = fhyperbel;
		init = init_hyperbel;
		funct_deriv = f_deriv;/*fhyperbel_deriv;*/
		numerical_flag = 1; /* no derivation function implemented yet */
		M = 2;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case HYPERBEL_LIN:
		/* f(x|a) = a1 / (x+a2) */
		fprintf( out, ", 1/y =  (x+a2)/ a1 = b1 + b2 * x" );
		fprintf( stderr, "\n 1/y =  (x+a2)/ a1 = b1 + b * x" );
		funct_deriv = flin_deriv;
		M = 2;
		lsParams->linear = 1;	/* nonlinear */
		break;
	case GAUSSIAN_1:
		/*
		 * f(x|a) = a1 * exp( a2 * (x-a3)^2)
		 */
		fprintf( out, ", a1 * exp( a2 * (x-a3)^2)" );
		fprintf( stderr, "\n a1 * exp( a2 * (x-a3)^2)" );
		funct = fgauss1;
		init = init_gauss;
		funct_deriv = fgauss_deriv;
		funct_deriv2 = fgauss_deriv2;
		M = 3;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case GAUSSIAN_2:
		/*
		 * f(x|a) = a1 * exp( a2 * (x-a3)^2) +
		 *          a4 * exp( a5 * (x-a6)^2)
		 */
		fprintf( out, ", a1 * exp( a2 * (x-a3)^2) + a4 * exp( a5 * (x - a6) ^ 2" );
		fprintf( stderr, "\n a1 * exp( a2 * (x-a3)^2) + a4 * exp( a5 * (x - a6) ^ 2" );
		funct = fgauss2;
		init = init_gauss2;
		funct_deriv = f_deriv;	/* use numerical differentiation */
		M = 6;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case GEN_LAPLACE:
		/* f(x|a) = a1 * exp( -|x|^a2 * a3) */
		fprintf( out, ", a1 * exp( -|x|^a2 * a3)" );
		fprintf( stderr, "\n a1 * exp( -|x|^a2 * a3)" );
		funct = fgen_laplace;
		init = init_gen_laplace;
		funct_deriv = fgen_laplace_deriv;
		M = 3;
		lsParams->linear = 0;	/* nonlinear */
		break;

	case QUAD_SURFACE:
		/*
		 * f(x|a) = a1 + a2*x1 + a3*x1^2 + a4*x2 +a5*x2^2
		 */
		fprintf( out, ", a1 + a2*x1 + a3*x1^2 + a4*x2 +a5*x2^2" );
		fprintf( stderr, "\n a1 + a2*x1 + a3*x1^2 + a4*x2 +a5*x2^2" );
		funct_deriv = fquadsurface_deriv;
		condDim = 2;
		M = 5;
		break;
	case COORD_TRANSF:
		/*
		 * f1(x|a) = a1 + cos(a3) * x1 - sin(a3) * x2
		 * f2(x|a) = a2 + sin(a3) * x1 + cos(a3) * x2
		 */
		fprintf( out, ", rotation f1(x|a) = a1 + cos(a3) * x1 - sin(a3) * x2" );
		fprintf( out, ", f2(x|a) = a2 + sin(a3) * x1 + cos(a3) * x2" );
		fprintf( stderr, "\n rotation" );
		funct = frotation;
		init = init_rotation;
		funct_deriv = frotation_deriv;
		M = 3;
		condDim = 2;
		obsDim = 2;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case CIRCLE:
		/*
		 * f(x|a) = 0 = (x1 - a1)^2 + (x2 - a2)^2 - a3*a3
		 */
		fprintf( out, ", circle: 0 = (x1 - a1)^2 + (x2 - a2)^2 - a3*a3" );
		fprintf( stderr, "\n circle: 0 = (x1 - a1)^2 + (x2 - a2)^2 - a3*a3" );
		funct = fcircle;
		init = init_circle;
		funct_deriv = fcircle_deriv;
		M = 3;
		condDim = 2;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case CIRCLE_LIN:
		/*
		 * f(x|a) = 0 = (x1 - a1)^2 + (x2 - a2)^2 - a3*a3
		 * f(x|b) = x1^2 + x2^2 = b1*x1 + b2*x2 - b3
		 * b1 = 2*a1, b2 = 2*a2, b3 = a1^2 + a2^2 - a3^2
		 */
		fprintf( out, ", circle, linearised: f(x|a) = 0 = (x1 - a1)^2 + (x2 - a2)^2 - a3*a3" );
		fprintf( stderr, "\n circle, linearised" );
		funct_deriv = fcirclelin_deriv;
		M = 3;
		condDim = 2;
		obsDim = 1;
		lsParams->linear = 1;	/* linear */
		break;
	case CIRCLE_TLS:
		/*
		 * f(x|a) = 0 = (sqrt[(x1 - a1)^2 + (x2 - a2)^2] - a3)^2
		 */
		fprintf( out, ", circle, TLS: 0 = (sqrt[(x1 - a1)^2 + (x2 - a2)^2] - a3)^2" );
		fprintf( stderr, "\n circle, TLS: 0 = (sqrt[(x1 - a1)^2 + (x2 - a2)^2] - a3)^2" );
		funct = fcircleTLS;
		init = init_circle;
		funct_deriv = fcircleTLS_deriv;
		M = 3;
		condDim = 2;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case CLASSIF_2D:	/* f(x|a) = a1 + a2*x1 + a3*x2 + a4*x1^2 + a5*x2^2 + a6*x1*x2 + a7*x1^3 + a8*x2^3 */
		fprintf( out, ", CLASSIF_2D: f(x|a) = a1 + a2*x1 + a3*x2 + a4*x1^2 + a5*x2^2 + a6*x1*x2 + a7*x1^3 + a8*x2^3" );
		fprintf( stderr, "\n CLASSIF_2D: f(x|a) = a1 + a2*x1 + a3*x2 + a4*x1^2 + a5*x2^2 + a6*x1*x2 + a7*x1^3 + a8*x2^3" );
		funct_deriv = fclassif2D_deriv;
		M = 8;
		condDim = 2;
		lsParams->linear = 1;	/* linear */
		break;

	case NN_3x3x1:
		/*
		 * f(x|a) = neural network 3x3x1
		 */
		fprintf( out, ", NN 3x3x1" );
		fprintf( stderr, "\n NN 3x3x1" );
		funct_deriv = f_deriv;	/* use numerical differentiation */
		funct = fNN_3_3;
		init = init_NN3x3x1;
		M = 16; /* in total 16 weights */
		condDim = 3;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NN_3x2x1:
		/*
		 * f(x|a) = neural network 3x2x1
		 */
		fprintf( out, ", NN 3x2x1" );
		fprintf( stderr, "\n NN 3x2x1" );
		funct_deriv = f_deriv;	/* use numerical differentiation */
		funct = fNN_3_2;
		init = init_NN;
		M = 11;
		condDim = 3;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NN_1x2x1:
		/*
		 * f(x|a) = neural network 1x2x1
		 */
		fprintf( out, ", NN 1x2x1" );
		fprintf( stderr, "\n NN 1x2x1" );
		funct_deriv = f_deriv;	/* use numerical differentiation */
		funct = fNN_1_2;
		init = init_NN;
		M = 7;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NN_2x2x1:
		/*
		 * f(x|a) = neural network 2x2x1
		 */
		fprintf( out, ", NN 2x2x1" );
		fprintf( stderr, "\n NN 2x2x1" );
		funct_deriv = f_deriv;	/* use numerical differentiation */
		funct = fNN_2_2;
		init = init_NN;
		M = 9;
		condDim = 2;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NN_1x3x1:
		/*
		 * f(x|a) = neural network 1x2x1
		 */
		fprintf( out, ", NN 1x3x1" );
		fprintf( stderr, "\n NN 1x3x1" );
		funct_deriv = f_deriv;	/* use numerical differentiation */
		funct = fNN_1_3;
		init = init_NN1x3x1;
		M = 10;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;

	case NIST_THURBER:
		/*
		 * f(x|a) =(a1 + a2*x + a3*x**2 + a4*x**3) /
		 *          (1 + a5*x + a6*x**2 + a7*x**3)
		 */
		fprintf( out, ", NIST_THURBER: (a1 + a2*x + a3*x**2 + a4*x**3) / (1 + a5*x + a6*x**2 + a7*x**3)" );
		fprintf( stderr, "\n NIST_THURBER: (a1 + a2*x + a3*x**2 + a4*x**3) / (1 + a5*x + a6*x**2 + a7*x**3)" );
		funct = fNIST_thurber;
		funct_deriv = fNIST_thurber_deriv;
		funct_deriv2 = fNIST_thurber_deriv2;
		init = init_NIST_thurber;
		M = 7;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NIST_MGH09:
		/*
		 * f(x|a) =a1 * (x**2 + a2*x) / (x*x + a3*x + a4)
		 */
		fprintf( out, ", NIST_MGH09: a1 * (x**2 + a2*x) / (x*x + a3*x + a4) " );
		fprintf( stderr, "\n NIST_MGH09: a1 * (x**2 + a2*x) / (x*x + a3*x + a4) " );
		funct = fNIST_MGH09;
		funct_deriv = fNIST_MGH09_deriv;
		funct_deriv2 = fNIST_MGH09_deriv2;
		init = init_NIST_MGH09;
		M = 4;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NIST_RAT42:
		/*
		 * f(x|a) = a1 / (1 + exp(a2 - a3*x))
		 */
		fprintf( out, ", NIST_Rat42: a1 / (1 + exp(a2 - a3*x)) " );
		fprintf( stderr, "\n NIST_Rat42: a1 / (1 + exp(a2 - a3*x)) " );
		funct_deriv = fNIST_Rat42_deriv;
		funct = fNIST_Rat42;
		init = init_NIST_Rat42;
		M = 3;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NIST_RAT43:
		/*
		 * f(x|a) = a1 / [1 + exp(a2 - a3*x)]^(1/a4)
		 * *  generalized logistic curve (Richards curve)
		 */
		fprintf( out, ", NIST_Rat43: a1 / [1 + exp(a2 - a3*x)]^(1/a4) " );
		fprintf( stderr, "\n NIST_Rat43: a1 / [1 + exp(a2 - a3*x)]^(1/a4) " );
		funct_deriv = fNIST_Rat43_deriv;
		funct = fNIST_Rat43;
		init = init_NIST_Rat43;
		M = 4;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NIST_ECKERLE4:
		/*
		 * f(x|a) = a1 / a2 * exp(-0.5*((x -a3)/ a2)^2)
		 */
		fprintf( out, ", NIST_ECKERLE4: a1 / a2 * exp(-0.5*((x -a3)/ a2)^2) " );
		fprintf( stderr, "\n NIST_ECKERLE4: a1 / a2 * exp(-0.5*((x -a3)/ a2)^2) " );
		funct_deriv = fNIST_Eckerle4_deriv;
		funct = fNIST_Eckerle4;
		init = init_NIST_Eckerle4;
		M = 3;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NIST_MGH10:
		/*
		 * f(x|a) = a1 * exp( a2 / (x+a3))
		 */
		fprintf( out, ", NIST_MGH10: a1 * exp( a2 / (x+a3)) " );
		fprintf( stderr, "\n NIST_MGH10: a1 * exp( a2 / (x+a3)) " );
		funct_deriv = fNIST_MGH10_deriv;
		funct_deriv2 = fNIST_MGH10_deriv2;
		funct = fNIST_MGH10;
		init = init_NIST_MGH10;
		M = 3;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NIST_BENNETT5:
		/*
		 * f(x|a) = a1 * (x+a2)^(-1/a3)
		 */
		fprintf( out, ", NIST_BENNETT5: a1 * (x+a2)^(-1/a3) " );
		fprintf( stderr, "\n NIST_BENNETT5: a1 * (x+a2)^(-1/a3) " );
		funct = fNIST_Bennett5;
		funct_deriv = fNIST_Bennett5_deriv;
		funct_deriv2 = fNIST_Bennett5_deriv2;
		init = init_NIST_Bennett5;
		M = 3;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NIST_BENNETT5_LOG:
		/*
		 * log(-f(x|a)) = log(-a1) + (-1/a3) * log(x+a2)
		 */
		fprintf( out, ", NIST_BENNETT5_LOG: log(-f(x|a)) = log(-a1) + (-1/a3) * log(x+a2) " );
		fprintf( stderr, "\n NIST_BENNETT5_LOG: log(-f(x|a)) = log(-a1) + (-1/a3) * log(x+a2) " );
		funct_deriv = fNIST_Bennett5LOG_deriv;	/* */
		funct = fNIST_Bennett5LOG;
		init = init_NIST_Bennett5LOG;
		M = 3;
		condDim = 1;
		obsDim = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case NIST_BOXBOD:
		/* f(x|a) = a1 *(1 - exp( -a2 * x) */
		fprintf( out, ", NIST_BOXBOD: a1 *(1 - exp( -a2 * x)" );
		fprintf( stderr, "\n NIST_BOXBOD: a1 *(1 - exp( -a2 * x)" );
		funct_deriv = fNIST_BoxBOD_deriv;
		funct = fNIST_BoxBOD;
		init = init_NIST_BoxBOD;
		M = 2;
		lsParams->linear = 0;	/* nonlinear */
		break;
	case EXAMPLE:
		/* f(x|a) = a1 *(1 - exp( -a2 * x) */
		fprintf( out, ", EXAMPLE: (a1*a1-1.5*a1)*(x-2)*(x-3)*(x-4)*(x-5)/24 - a1*(x-1)*(x-3)*(x-4)*(x-5)/ 6" );
		fprintf( stderr, "\n EXAMPLE: (a1*a1-1.5*a1)*(x-2)*(x-3)*(x-4)*(x-5)/24 - a1*(x-1)*(x-3)*(x-4)*(x-5)/ 6" );
		funct_deriv2 = fExample_deriv2;
		funct_deriv = fExample_deriv;
		funct = fExample;
		init = init_Example;
		M = 1;
		lsParams->linear = 0;	/* nonlinear */
		break;
	default:
		err = errmsg( ERR_NOT_DEFINED, rtn, "-m ", type );
		usage( argv[0] );
		goto endfunc;
	}

	if (lsParams->expression_flag)
	{
		fprintf( out, ", TinyExpr: variable function, %s", funcStringInput );
		if (renumberVariables( funcStringInput, funcString, out ))
		{
			fprintf( out, "\n# corrected function string: %s", funcString );
		}
		if (lsParams->positiveParams)
		{
			fprintf( out, "\n# !Attention: if option '-p' is used, function must be defined" );
			fprintf( out, "\n#             such that all parameters aj are non-negative!" );
		}
		fprintf( stderr, "\n# TinyExpr: variable function, %s", funcString );
		{
			const double* vars[M_MAX];
			int count = 0;
			int errPos;
			teExpr = te_compile( funcString, teVars, sizeof( teVars ) / sizeof( te_variable ), &errPos );
			if (!teExpr)
			{
				fprintf( stderr, "\n te_compile: error at position %d", errPos );
				fprintf( stderr, "\n %s", funcString );
				fprintf( stderr, "\n %*s\nError occurred here.\n", errPos, "^" );
				fprintf( out, "\n# te_compile: error at position %d", errPos );
				fprintf( out, "\n# %s", funcString );
				fprintf( out, "\n# %*s\nError occurred here.\n", errPos, "^" );
				fclose( out );
				return 1;
			}

			collect_vars( teExpr, vars, &x, &count );

			//teExpr->M = count;
			M = count;
		}
		lsParams->linear = isLinearExpr( teExpr ).isLinear;
		if (numerical_flag == 0)
		{
			fprintf( stderr, "\n model function is defined using option `-y'" );
			fprintf( stderr, "\n ==> derivatives must be calculated numerically" );
			fprintf( out, "\n# model function is defined using option `-y'" );
			fprintf( out, "\n# ==> derivatives must be calculated numerically" );
			numerical_flag = 1; /* overwrite other settings */
		}
		//if (lsParams->trueH)
		//{
		//	fprintf( stderr, "\n numerical calculation of second derivative is not recommended" );
		//	fprintf( stderr, "\n option '-H' is deactivated" );
		//	fprintf( out, "\n# numerical calculation of second derivative is not recommended" );
		//	fprintf( out, "\n# option '-H' is deactivated" );
		//	lsParams->trueH = 0; /* second derivative is not supported */
		//}
		// f_derivTinyExpr() /* use numerical differentiation */
		init = init_tinyExp; /* initialization of parameters */
	}

	if (numerical_flag && lsParams->linear)
	{
		if (lsParams->expression_flag == 0)
		{
			/* the model functions of linear problems do not exist */
			fprintf( stderr, "\n Numerical derivation is not implemented" );
			fprintf( stderr, "   for linear model functions! -------\n" );
			fprintf( out, "\n# Numerical derivation is not implemented" );
			fprintf( out, "   for linear model functions! -------\n#" );
			numerical_flag = 0;
		}
		else
		{/* initialisation of parameters required !! */
			for (unsigned j = 0; j < M; j++)
			{
				a[j] = 1.0; /* may not be zero */
			}
		}
	}
	if (numerical_flag)
	{
		funct_deriv = f_deriv;/* use numerical differentiation */
		funct_deriv2 = f_deriv2;/* use numerical differentiation */
	}

	if (funct_deriv == NULL)
	{
		fprintf( out, "\n# funct_deriv() not defined \n" );
		fprintf( stderr, "\n funct_deriv() not defined \n" );
		err = ERR_ANY;
		goto endfunc;
	}

	if (lsParams->linear)
	{
		if (!lsParams->LM)
		{
			fprintf( stderr, "\n The optimization problem is linear" );
			fprintf( stderr, "\n ==> option `-G' will be ignored" );
			if (lsParams->GNspecial_flag)
			{
				fprintf( stderr, "\n ==> option `-e' will be ignored" );
				lsParams->GNspecial_flag = 0;
			}
			if (lsParams->DM)
			{
				fprintf( stderr, "\n ==> option `-D' will be ignored" );
				lsParams->DM = 0;
			}
		}
	}
	else /* is nonlinear */
	{
		if (!lsParams->svd)
		{
			fprintf( stderr, "\n The optimization problem is nonlinear" );
			fprintf( stderr, "\n ==> option `-s' will be ignored" );
		}
		if (lsParams->LM)
		{
			if (lsParams->GNspecial_flag)
			{
				fprintf( stderr, "\n Levenberg-Marquardt is active" );
				fprintf( stderr, "\n ==> option `-e' will be ignored" );
				lsParams->GNspecial_flag = 0;
			}
		}
		else
		{
			if (lsParams->DM)
			{
				fprintf( stderr, "\n Levenberg-Marquardt is NOT active ('-G')" );
				fprintf( stderr, "\n ==> option `-D' will be ignored" );
				lsParams->DM = 0;
			}
		}
	}

	/* if column for observation is not given explicitely by a
	 * command-line parameter, then assume the column following
	 * the conditions
	 */
	if (columnObs == 0) columnObs = condDim + 1;
	fprintf( stderr, "\n" );
	fflush( stdout );

	if (M > 5 && algoMode == 0)
	{
		err = errmsg( ERR_TOO_MANY_PARAMETERS, rtn, "for standard matrix inversion", M );
		usage( argv[0] );
		goto endfunc;
	}

	if (M > MAX_CONDITIONS)
	{
		err = errmsg( ERR_TOO_MANY_PARAMETERS, rtn, "maximum is exceeded", M );
		goto endfunc;
	}

	if (N <= M)
	{
		err = errmsg( ERR_TOO_MANY_PARAMETERS, rtn, "compared to number of observations", M );
		goto endfunc;
	}
	fprintf( out, "\n# Number of observations: %d", N );
	fprintf( out, "\n# Number of parameters  : %d", M );

	if (lsParams->linear && lsParams->svd && algoMode != 1)
	{
		fprintf( stderr, "\n# option '-a %d' is ignored, ", algoMode );
		fprintf( stderr, "since special SVD approach is used!\n" );
		algoMode = 1;
	}
	if (lsParams->trueH && funct_deriv2 == NULL)
	{
		fprintf( stderr,
			"\n### function for 2nd derivativ was not initialised!" );
		err = errmsg( ERR_ANY, rtn, "funct_deriv2()", 0 );
		goto endfunc;
	}

	/*
	 * allocate memory
	 */
	jacob = matrix( N * obsDim, M );	/* Jacobian */
	covar = matrix( M, M );	/* covariance matrix */

	observations = vector( N * obsDim );	/* observations */
	datac = vector( N * obsDim );	/* calculated data using
																		 f(x|a) */
	conditions = vector( N * condDim );	/* conditions x */
	weights = vector( N * obsDim );	/* weights */
	weightsOld = vector( N * obsDim );	/* weights one step back
																				 */
	deviates = vector( N * obsDim );	/* remaining differences */
	/* remaining absolute differences */
	deviatesAbsolute = vector( N * obsDim );
	//deltasq = vector( N * obsDim );	/* remaining squared
																			 //differences */

																			 /* open input file again */
	in = fopen( inname, "rt" );
	if (in == NULL)
	{
		err = errmsg( ERR_OPEN_READ, rtn, inname, 0 );
		perror( "\nReason" );
		fprintf( out, "\n#  Cannot open %s \n", inname );
		goto endfunc;
	}

	fprintf( out, "\n# -cc , condition columns: " );
	for (unsigned int d = 0; d < condDim; d++)
	{
		fprintf( out, "%d ", columnCond[d] );
	}

	fprintf( out, "\n# -co , observations column: " );
	for (unsigned int d = 0; d < obsDim; d++)
	{
		fprintf( out, "%d ", columnObs + d );
	}

	if (weightMode == 3)
	{
		if (columnUncertainties == 0)
		{
			fprintf( out, "\n# -cu missing for weighting mode -w 3 !!!" );
			fprintf( out, "\n#     fallback to weighting mode -w 0" );
			weightMode = 0;
			fprintf( stderr, "\n#### column of uncertainties must be > 0 for weighting mode -w 3 ! #######\n" );
			fprintf( stderr, "\n#### Option -cu must be set correctly !                            #######\n" );
		}
		else
		{
			fprintf( out, "\n# -cu , column of uncertainties: %d", columnUncertainties );
		}
	}

	fprintf( out, "\n#" );
	/* write LS flags in output	*/
	if (lsParams->linear)
	{
		fprintf( out, "\n# fitting a linear system" );
		if (lsParams->svd)
			fprintf( out, ", use special SVD based algorithm" );
		else
			fprintf( out, ", -s do not use special SVD based algorithm" );
	}
	else
	{
		fprintf( out, "\n# fitting a nonlinear system" );
		if (lsParams->LM)
		{
			fprintf( out, "\n# use Levenberg-Marquardt method" );
			if (lsParams->DM)
				fprintf( out, "\t -D , with diag(J^T*W*J)" );
			else
				fprintf( out, "\t with identity matrix" );
		}
		else
		{
			fprintf( out, "\n# -G , use Gauss-Newton method" );
			if (lsParams->GNspecial_flag) fprintf( out, " (-e , enhanced)" );
		}
		if (lsParams->trueH)
		{
			fprintf( out, "\n# -H , use true Hessian matrix" );
		}
		if (lsParams->modifyVanishingGradients_flag)
		{
			fprintf( out, "\n# -g , allow gradient manipulation at saddle points" );
		}
		if (lsParams->chisqTarget_flag)
			fprintf( out, "\n# -t %f , chisq must be lower than given value",
				chisqTargetValue );
		if (numerical_flag)
		{
			fprintf( out, "\n# -n , numerical derivation activated" );
		}
	}

	fprintf( out, "\n# -a %d , algorithm for inversion: ", algoMode );
	if (lsParams->linear && lsParams->svd)
		fprintf( out, "not needed" );
	else if (algoMode == 0)
		fprintf( out, "Cofactor" );
	else if (algoMode == 1)
		fprintf( out, "SVD" );
	else if (algoMode == 2)
		fprintf( out, "LU decomposition" );
	//fprintf( out, "\n# -u %.2f , probability of uphill step: ", lsParams->uphillProbability );

	/* put statement about weighting and outlier detection scheme	*/
	fprintf( out, "\n# -w %d , mode of weighting:", weightMode );
	switch (weightMode)
	{
	case 0:		fprintf( out, "  equal weights" ); break;
	case 1:		fprintf( out,
		"  estimate weights based on deviates" ); break;
	case 2:		fprintf( out,
		"  estimate weights based on binning" ); break;
	case 3:		fprintf( out, "  uncertainties have been provided" );
		break;
	default: break;
	}
	if (forget_flag)
		fprintf( out, "\n# reset weights after outlier removal" );
	fprintf( out, "\n# -x %d , mode of outlier detection:", outMode );
	switch (outMode)
	{
	case 0:	fprintf( out, "  no outlier detection" ); break;
	case 1:	fprintf( out,
		"  based on z-score and Chauvenet's criterion" );
		break;
	case 2:	fprintf( out, "  based on ClubOD" ); break;
	case 3:	fprintf( out,
		"  based on M-score and Chauvenet's criterion" );
		break;
	case 4:	fprintf( out, "  based on RANSAC" ); break;
	default: break;
	}


	/* read the conditions and observations */
	for (unsigned int i = 0; i < N; i++)
	{
		/* jump over comments and empty lines */
		do
		{
			ptr = fgets( line, MAXLINELENGTH, in );
		} while (!is_data_line( line, MAXLINELENGTH ));

		/* loop over all conditions */
		for (unsigned int d = 0; d < condDim; d++)
		{
			/* if 0 was given, then assume that conditions are just
			 * serial numbers 1,2,3,...
			 */
			if (columnCond[d] == 0)
			{
				if (condDim > 1)
				{
					fprintf( stderr,
						"\n There is more than one condition (%d)!", condDim );
					fprintf( stderr,
						"\n Columns of conditions must be given via '-cc'!\n" );
					goto endfunc;
				}
				conditions[condDim * i + d] = i + 1;
			}
			else
			{
				/* get string starting from desired column */
				field = get_nth_field( line, columnCond[d] );
				if (field != NULL)
				{
					/* multidimensional conditions are stored one after each
						 other */
					res = sscanf( field, "%lf", &(conditions[condDim * i + d]) );
				}
				else
				{
					fprintf( stderr, "\n\n === %d th column does not exist",
						columnCond[d] );
					err = errmsg( ERR_COLUMN_DOES_NOT_EXIST, rtn, "option -cc ", columnCond[d] );
					goto endfunc;
				}
			}
		}
		/* loop over all observations */
		for (unsigned int d = 0; d < obsDim; d++)
		{
			/* get string starting from desired column */
			field = get_nth_field( line, columnObs + d );
			if (field != NULL)
			{
				res = sscanf( field, "%lf", &(observations[obsDim * i + d]) );
			}
			else
			{
				err = errmsg( ERR_COLUMN_DOES_NOT_EXIST, rtn, "option -co ", columnObs );
				fprintf( stderr, "### see line %d\n", i );
				goto endfunc;
			}
		}

		if (columnUncertainties > 0)
		{
			/* get string starting from desired column */
			field = get_nth_field( line, columnUncertainties );
			if (field != NULL)
			{
				double colUncertainty;
				res = sscanf( field, "%lf", &(colUncertainty) );

				/* all observations per row i share the same uncertainty */
				for (unsigned int d = 0; d < obsDim; d++)
				{
					weights[obsDim * i + d] = 1 / (colUncertainty * colUncertainty);
				}
			}
			else
			{
				fprintf( stderr,
					"\n\n === %d th column does not exist in line %d",
					columnUncertainties, i );
				err = errmsg( ERR_COLUMN_DOES_NOT_EXIST, rtn, "option -cu ", columnUncertainties );
				fprintf( stderr, "### see line %d\n", i );
				goto endfunc;
			}
		}
	}
	fclose( in );

	/* sort the conditions and observations */
	if (condDim < 2)
	{ /* sorting is not possible for more than one dimension of conditions */
		long* idx = lvector( N );
		double* ptrObs = NULL;
		double* newObs = vector( N );
		heap_sort_d_( N, conditions, idx );
		for (unsigned int i = 0; i < N; i++)
		{
			newObs[i] = observations[idx[i]];
		}
		ptrObs = observations;
		observations = newObs; /* use sorted vector of observations */
		free_lvector( &idx );
		free_vector( &ptrObs ); /* free original vector */
	}

	/*---------------------------------------------
	 * scaling of conditions, if enabled
	 */
	if (scaling_flag)
	{
		scalingFactor = conditionScaling( N, condDim, conditions, type, out );
	}

	/* initialize weights */
	if (columnUncertainties == 0) /* no uncertainties given	*/
	{
		for (unsigned int i = 0; i < N * obsDim; i++)
		{
			/* sum of all weights must be equal to N minus number of
			 * outliers
			 */
			weights[i] = 1.0;
			weightsOld[i] = 1.0;
		}
		meanOfWeights = 1; // initialization ist required
	}
	else
	{
		meanOfWeights = 0;
		for (unsigned int i = 0; i < N; i++)
		{
			meanOfWeights += weights[i]; // initialization ist required
		}
		meanOfWeights /= N;
	}

	/* data preparation for some model functions */
	specialPreparations( N, conditions,
		observations, weights, type );

	/* logging information */
	fprintf( out, "\n# DBL_DECIMAL_DIG: %d", DBL_DECIMAL_DIG );
	fprintf( out, "\n# DBL_EPSILON....: %.*e", DBL_DECIMAL_DIG, DBL_EPSILON );

	/*
	 * set initial parameters for nonlinear functions
	 * parameter values given on command line wont be changed
	 */
	if (!lsParams->linear)
	{
		init( N, observations, conditions, a, a_flag, out );
		if (lsParams->expression_flag)
		{ /* try to correct parameters if necessary */
			double tmp;
			int stop = 1;
			unsigned int j = 0;
			do
			{
				tmp = te_eval( teExpr );
				if (!IsFiniteNumber( tmp ))
				{
					a[j] = -a[j];
					stop = 0;
				}
				j++;
			} while (j < M && !stop);
		}

		fprintf( out, "\n# initial Parameters\n# " );
		/* write initial parameters to output */
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "a%d=%.14G  ", j + 1, a[j] );
		}
		for (unsigned int j = M; j < M_MAX; j++)
		{
			/* zero out unnecessary parameters
			 * required for POLYNOMIAL_REG
			 */
			a[j] = 0.;
		}

		/* If scaling is enabled, all **initial** parameters, which are set
		 * independently on the condition values, must be corrected
		 */
		 //if (scaling_flag)
		 //{
		 //	condition_unscalingN( M, a, scalingFactor, type, out );
		 //}
	} /* non-linear */


	/*
	 *  estimation of weights if required
	 */
	iterStop_flag = 0;
	if (columnUncertainties > 0) /* weights given	*/
	{
		if (weightMode != 3)
			fprintf( out,
				"\n! weights have been read from file. set weight_mode = 3 !\n" );
		weightMode = 3; /* overwrite weights mode	*/
	}

	if (weightMode == 0 || weightMode == 3)
	{
		iterStop_flag = 1; /* only one run */
		iterFinal_flag = 1; /* only one run for lsLinear and outlier removal */
	}
	else if (weightMode == 2)
	{
		/* weights can be estimated beforehand via binning */
		est_weights2( N * obsDim, conditions, observations, weights,
			obsPerBin, out );
		iterStop_flag = 1; /* only one run */
		iterFinal_flag = 1;	/* only one run of least squares */
	}
	else
		iterFinal_flag = 0;

	/* outlier detection has not be performed yet */
	outDetect_flag = 0;

	//if (errno)
	//{
	//	perror( "\n### " );
	//}
	/* ----------------------------------------------------- */
	/* loop for weights estimation */
	for (int iterW = 0;
		(iterW <= iterWmax /* cont. as long max. number of
													iterations has not reached			*/
			&& !iterStop_flag)	/* the stop flag is not set		*/
		|| (iterFinal_flag);/* or it is the last round		*/
		iterW++)
	{
		if (iterFinal_flag && iterStop_flag)
			iterFinal_flag = 0; /* if final round reached, then reset flag */

		/* feedback on console */
		fprintf( stderr, "\r weights iterations: %3d\n", iterW );

		/* estimate weights in all but the last iteration */
		fprintf( out, "\n#\n#====================================" );
		if (!iterStop_flag)
			fprintf( out, "\n# %s: weights iteration #: %d", rtn, iterW );
		else if (weightMode)
		{
			fprintf( out,
				"\n# %s: approximation with final weights #: %d", rtn, iterW );
		}
		if (errno)
		{
			char str[511];
			snprintf( str, 510, "\n### %s: ", rtn );
			perror( str );
		}

		/* do the least squares approximation		*/
		if (lsParams->linear)
		{
			err =
				lsLinear( funct, funct_deriv,
					teExpr, N * obsDim, M,
					observations, conditions, jacob, weights,
					a, algoMode, lsParams,
					covar, out );
		}
		else
		{
			err =
				lsNonlinear( funct, funct_deriv, funct_deriv2,
					teExpr, N * obsDim, M,
					observations, conditions, jacob, weights,
					a, algoMode, lsParams,
					chisqTargetValue, covar, &iterations, out,
					digitsIntermediate );
		}
		if (err)
		{
			errmsg( err, rtn, "ups", 0 );
		}
		/* compute weighted and squared differences, chi-squared */
		chisq = energy = mean = 0.0;
		cntInliers = 0;
		if (lsParams->linear)
		{
			/* separate for linear and nonlinear, because funct() is
				 not defined for linear models */
			for (unsigned int i = 0; i < N * obsDim; i++)
			{
				double deltasq;
				if (lsParams->expression_flag)
				{ /* is more exact than using the Jacobian, since it was computed numerically */
					x = conditions[i];
					datac[i] = te_eval( teExpr );
				}
				else
				{
					/* get calculated data points dependent on current
						 parameters */
					datac[i] = 0.0; /* model functions are not implemented */
					for (unsigned int j = 0; j < M; j++)
					{
						//datac[i] += a[j] * jacob[i][j];
						datac[i] = fma( a[j], jacob[i][j], datac[i] );
					}
				}
				deviates[i] = observations[i] - datac[i];
				deviatesAbsolute[i] = fabs( deviates[i] );
				/* weighted and squared differences */
				deltasq = deviates[i] * deviates[i];
				if (weights[i] > 0.)
				{
					/* exclude outliers, i.e. weigths == 0 */
					//chisq += weights[i] * deltasq[i];
					chisq = fma( weights[i], deltasq, chisq );
					energy += deltasq;
					mean += deviatesAbsolute[i];
					cntInliers++;
				}
			}
		}
		else /* if not linear */
		{
			for (unsigned int i = 0; i < N * obsDim; i++)
			{
				double deltasq;
				/* get calculated data points dependent on current
					 parameters */
				if (lsParams->expression_flag)
				{
					x = conditions[i];
					datac[i] = te_eval( teExpr );
				}
				else
				{
					datac[i] = funct( i, conditions, a );
				}
				if (errno)
				{
					char str[511];
					snprintf( str, 510, "\n### %s: ", rtn );
					perror( str );
				}

				deviates[i] = observations[i] - datac[i];
				deviatesAbsolute[i] = fabs( deviates[i] );
				/* weighted and squared differences */
				deltasq = deviates[i] * deviates[i];
				if (weights[i] > 0.)
				{
					/* exclude outliers in final iteration*/
					//chisq += weights[i] * deltasq[i];
					chisq = fma(weights[i], deltasq, chisq);
					energy += deltasq;
					mean += deviatesAbsolute[i];
					cntInliers++;
				}
			}
		}
		numOutliers = N * obsDim - cntInliers;

		/* estimate of k,  w_i= k/sigma^2_i */
		gfit = chisq / (double)(cntInliers - M);	/* goodness of fit	*/
		mean = mean / (double)(cntInliers);
		variance = energy / (double)cntInliers - mean * mean;

		if (weightMode && !iterStop_flag)
		{ /* intermediate output */
			fprintf( out, "\n#\n# %s\n# Parameters:  ", rtn );
			for (unsigned int j = 0; j < M; j++)
			{
				fprintf( out, "a%d=%.6f, ", j + 1, a[j] );
			}

			fprintf( out, "\n#\n# | chi-square........: %f", chisq );
			fprintf( out, "\n# | mean of |deviates|: %f", mean );
			fprintf( out, "\n# | variance of |deviates|: %f", variance );
			fprintf( out, "\n# | goodness of fit...: %f", gfit );
			fprintf( out, "\n# | number of outliers: %d", numOutliers );
		}

		/* estimate weights, but not in last iteration */
		if (!iterStop_flag && weightMode)
		{
			fprintf( out, "\n#\n# enter weight estimation" );
			/* estimation of weights based on absolute deviates */
			if (weightMode == 1)
			{
				est_weights1( N * obsDim, deviatesAbsolute, weights, out );
			}
			/* no iterative weighting in weightMode == 2 | 3 */

			fprintf( out, "\n#\n# %s\n#  i       observ         ", rtn );
			fprintf( out, "calc      deviates    weights" );

			/* get mean of weights and output current values */
			meanOfWeights = 0;
			cntInliers = 0;
			for (unsigned int i = 0; i < N * obsDim; i++)
			{
				meanOfWeights += weights[i];	/* compute sum of all weights */
				if (i < MAX_LINES) /* limits the number of output lines */
				{
					fprintf( out, "\n# %2d  %12.5f %12.5f %12.5f %14f", i,
						observations[i], datac[i], deviatesAbsolute[i], weights[i] );
				}
				if (weights[i] > 0.0)
				{
					cntInliers++;	/* count used observations	*/
				}
			}
			/* mean of all weights > 0.; it is used later on */
			meanOfWeights /= (double)cntInliers;

			/* compare new and old weights */
			sum = 0;
			for (unsigned int i = 0; i < N * obsDim; i++)
			{
				/* watch changes of weights */
				sum += fabs( weights[i] - weightsOld[i] );
				weightsOld[i] = weights[i];	/* remember for next
																				 iteration */
			}
			sum /= (double)N * obsDim;	/* mean difference in
																			 weights */
																			 /* criterion of convergence */
			if (sum < 0.0001)
			{
				if (!iterStop_flag)
				{
					iterFinal_flag = 1;	/* go to last iteration */
				}
				iterStop_flag = 1; /* last iteration has been performed */

				fprintf( out, "\n#\n# convergence of weights" );
			} /* if (sum < 0.0001)*/

			if (iterW == iterWmax && !iterStop_flag)
			{
				fprintf( out,
					"\n#\n# maximum number of iterations reached" );
				fprintf( out, "\n# no convergence of weights" );
				iterFinal_flag = 1;	/* go to last iteration */
			} /* if (iterW == iterWmax && !iterStop_flag)*/
		}	 /* if (!iterStop_flag && weightMode) */

#ifdef OUTPUT_DEVIATES
		/* write deviates in separate file */
		if (iter_stop == 1)
		{
			char dev_name[500];
			int i, len;
			FILE* out_dev;

			len = strlen( outname );
			/* copy filename w/o extension */
			for (i = 0; i < len; i++)
			{
				if (outname[i] == '.') break;
				dev_name[i] = outname[i];
			}
			dev_name[i] = '_';
			dev_name[i + 1] = 'd';
			dev_name[i + 2] = 'e';
			dev_name[i + 3] = 'v';
			dev_name[i + 4] = '.';
			dev_name[i + 5] = 'x';
			dev_name[i + 6] = 'y';
			dev_name[i + 7] = '\0';
			/* open out file */
			out_dev = fopen( dev_name, "wt" );
			if (out_dev == NULL)
			{
				err = errmsg( ERR_OPEN_WRITE, rtn, dev_name, 0 );
				goto endfunc;
			}
			fprintf( out_dev, "#deviates for file " );
			fprintf( out_dev, "%s ", outname );
			fprintf( out_dev, "\n# before outlier detection" );
			fprintf( out_dev, "\n# i          deviate" );
			if (type == COORD_TRANSF)
			{
				fprintf( out_dev, "X     deviateY" );
				for (i = 0; i < N * obs_dim; i += 2)
				{
					fprintf( out_dev, "\n%4d \t %.4e\t %.4e", i,
						deviate[i], deviate[i + 1] );
				}
			}
			else
			{
				for (i = 0; i < N * obs_dim; i++)
				{
					fprintf( out_dev, "\n%4d \t %.4e", i, deviate[i] );
				}
			}
			fclose( out_dev );
		} /* end of deviates output */
#endif

		/* enter outlier detection only, if it was not done before */
		if (iterStop_flag && outMode && (!outDetect_flag))
		{
			numOfOutliers = 0;
			outDetect_flag = 1;
			/* do the outlier detection */
			if (outMode == 1)
			{
				/* do the outlier detection via z-score */
				/* we can exploit the value of gfit,
				 * because sigma^2 = gfit / mean(weigts)
				 */
				numOfOutliers = outlier_detection1( N * obsDim,
					sqrt( gfit / meanOfWeights ), deviatesAbsolute,
					weights, 0.15,
					out );
			}
			else if (outMode == 2)
			{
				/* do the cluster-based outlier detection */
				numOfOutliers = outlier_detection2( N * obsDim,
					deviatesAbsolute, weights, out );
			}
			else if (outMode == 3)
			{
				/* do the MAD outlier detection */
				numOfOutliers = outlier_detection3( N * obsDim,
					deviatesAbsolute, weights, 0.15, out );
			}
			else if (outMode == 4)
			{
				/* do the RANSAC outlier detection based on
				 * least-squares approximation	on subsets
				 */
				if (obsDim == 1)
				{
					numOfOutliers =
						ransac( funct, funct_deriv, funct_deriv2,
							teExpr, N, M,
							observations, conditions, jacob, weights, a, algoMode,
							lsParams, chisqTargetValue, covar, out, deviatesAbsolute,
							condDim, obsDim, digitsIntermediate );
					/* do a last round of lsLinear in order to get a correct
					 * covariance matrix
					 */
					iterFinal_flag = 1;
				}
				else
				{
					fprintf( out,
						"# RANSAC is not implemented for multi-variate data\n" );
				}
			}
			/* if outliers have been found, do an additional final
			 * least square approx.
			 */
			if (numOfOutliers > 0)
			{
				iterFinal_flag = 1;
				if (forget_flag == 1)
				{
					fprintf( out, "\n#forget weights" );
					/* set  all weights to 1. */
					for (unsigned int i = 0; i < N * obsDim; i++)
					{
						if (weights[i] > 0.0)
							weights[i] = (float)1.0;
					}
					meanOfWeights = 1.0;
				}
			}
		} /* if outMode */
	}		 /* for iterW */
	/* ----------------------------------------------------- */

	/*
	 * evaluation of results
	 */

	 /*
		* since the weights are already normalised to their
		* mean value (w = k/sigma^2), gfit is also equal to
		* the variance of observations
		*/

		/* uncertainty in observations, root of sample variance, RMSE */
	uncertainty = sqrt( gfit / meanOfWeights );

	fprintf( out, "\n#\n# evaluation of results" );
	fprintf( out, "\n#  number of outliers: %d", numOutliers );
	fprintf( out, "\n#  parameters:" );
	for (unsigned int j = 0; j < M; j++)
	{
		fprintf( out, " a%d=%f", j + 1, a[j] );
	}
	fprintf( out, "\n#  chi-square.................: %.*e", OP_DBL_Digs, chisq );
	fprintf( out, "\n#  goodness of fit............: %.12G", gfit );
	fprintf( out, "\n#  uncertainty in observations: %.12G  (weighted) RMSE",
		uncertainty );
	fprintf( out, "\n#\n#  (co)variance of parameters:" );
	{
		int flag = 0;
		for (unsigned int i = 0; i < M; i++)
		{
			fprintf( out, "\n#" );
			for (unsigned int j = 0; j < M; j++)
			{
				covar[i][j] *= gfit;
				fprintf( out, " %15.9G", covar[i][j] );
				if ((i == j) && (covar[i][j] < 0.)) flag = 1;
			}
		}
		if (flag)
		{
			fprintf( stderr, "\n#  negative parameter variance" );
			fprintf( stderr, "\n#  probably ill-conditioned problem" );
			if (!lsParams->svd)
			{
				fprintf( stderr, "\n#  solution is probably wrong" );
				fprintf( stderr, "\n#  disable option '-s'" );
				fprintf( out, "\n####  solution is probably wrong       #" );
				fprintf( out, "\n####  disable option '-s'" );
			}
			fprintf( out, "\n####  negative parameter variance      #" );
			fprintf( out, "\n####  probably ill-conditioned problem #" );
		}
		/* check correlation between parameters */
		for (unsigned int i = 0; i < M; i++)
		{
			for (unsigned int j = i + 1; j < M; j++)
			{
				double rho = covar[i][j] / sqrt( covar[i][i] * covar[j][j] );
				if (rho > 0.9)
				{
					fprintf( stderr, "\n#  there is suspicious correlation between parameters a%d and a%d", i + 1, j + 1 );
					fprintf( stderr, "\n#  check model function!" );
					fprintf( out, "\n####  there is suspicious correlation between parameters a%d and a%d", i + 1, j + 1 );
					fprintf( out, "\n####  check model function!" );
				}
			}
		}

	}

	/*-------------------------------------------------------
	 * post-processing
	 */
	specialPostprocessing( M, a, type, out );

	/* check the uncertainty in parameters */
	checkParameterUncertainties( M, a, covar, out, outname );

	if (numOutliers)
	{
		fprintf( stdout, "\n# detection of %d outliers!", numOutliers );
	}

	/* unscale the parameters	*/
	if (scaling_flag)
	{
		condition_unscaling( N, conditions, condDim, M, a, scalingFactor, type, out );
	}

	fprintf( out, "\n#\n# Final_Parameters " );
	for (unsigned int j = 0; j < M; j++)
	{
		//fprintf( out, "a%d= %.12G ", j+1, a[j]);
		/* output with maximum precision */
		fprintf( out, "\n# a%d= %.*e", j + 1, OP_DBL_Digs, a[j] );
	}
	if (lsParams->expression_flag)
	{
		char finalFunction[511];
		createFinalFunction( funcString, finalFunction, a, out );
		fprintf( out, "\n# fy(x) = %s", finalFunction );
	}
	else
	{
		if (type == POLYNOMIAL)
		{
			fprintf( out, "\n# f%d(x) = %.14G", M, a[0] );
			if (M > 1)
			{
				fprintf( out, " + %.14G*x", a[1] );
			}
			for (unsigned int j = 2; j < M; j++)
			{
				fprintf( out, " + %.14G * x**%d", a[j], j );
			}
		}
	}

	/* map parameters to original model function */
	if (type == COSINE_LIN)
	{
		double r0, phi0;

		/* f(x|b) = b1 + r0 * cos( x - phi0)
		 * f(x|a) = a1 + a2 * cos( x) + a3 * sin( x)
		 */
		r0 = sqrt( a[1] * a[1] + a[2] * a[2] );	/* radius a2 */
		phi0 = 0;
		if (r0 > 0.)
		{
			double pc1, pc2, ps1, ps2;
			double d11, d12, d21, d22;

			/* solve ambiguity of angles */
			pc1 = acos( a[1] / r0 );	/* 1st solution */
			pc2 = -pc1;	/* 2nd solution */
			ps1 = asin( a[2] / r0 );	/* 3rd solution */
			if (ps1 < 0)	/* 4th solution */
				ps2 = -M_PI - ps1;
			else
				ps2 = M_PI - ps1;

			d11 = fabs( pc1 - ps1 );	/* two of four must be equal */
			d12 = fabs( pc1 - ps2 );	/* take differences */
			d21 = fabs( pc2 - ps1 );
			d22 = fabs( pc2 - ps2 );
			/* look for smallest difference */
			if (d11 < d12 && d11 < d21 && d11 < d22)
			{
				phi0 = 0.5 * (pc1 + ps1);
			}
			if (d12 < d11 && d12 < d21 && d12 < d22)
			{
				phi0 = 0.5 * (pc1 + ps2);
			}
			if (d21 < d11 && d21 < d12 && d21 < d22)
			{
				phi0 = 0.5 * (pc2 + ps1);
			}
			if (d22 < d11 && d22 < d12 && d22 < d21)
			{
				phi0 = 0.5 * (pc2 + ps2);
			}
			a[1] = r0;
			double tmp = phi0 * 180. / M_PI;
			// does not work: a[2] = atan2( a[2], a[1] ) * 180. / M_PI;
			a[2] = tmp;
			//if (fabs( 1000 * tmp - 1000 * a[2] ) > 0.0)
			//{
			//	fprintf( out, "\n#\n# Aaaarrrgh! phi0=%f", phi0 );
			//	fprintf( stderr, "\n#\n# Aaaarrrgh! phi0=%f", phi0 );
			//}
		}
		else
		{
			fprintf( out, "\n#\n# function is constant! ######## " );
			fprintf( stderr, "\n#\n# function is constant! ######## " );
			a[1] = 0.;
			a[2] = 0.;
		}
		fprintf( out, "\n#----------------------------------------" );
		fprintf( out, "\n# corrected Parameters " );
		fprintf( out, "according to f(x)=b1+b2*cos(x-b3)\n# " );
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "b%d=%.9f, ", j + 1, a[j] );
		}
		fprintf( out, "\n#" );
	}
	else if (type == EXPONENTIAL2_LIN)
	{
		/* convert observations back */
		for (unsigned int i = 0; i < N; i++)
		{
			observations[i] = exp( observations[i] );
		}
		a[0] = exp( a[0] );

		fprintf( out, "\n#\n# corrected Parameters\n# " );
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "a%d=%.9f, ", j + 1, a[j] );
		}
		fprintf( out, "\n#\n# uncertainties of corrected Parameters " );
		fprintf( out, "are not available yet\n# " );

		for (unsigned int i = 0; i < N; i++)
		{
			/* get calculated data points dependent on corrected
				 parameters */
			datac[i] = fexpon2( i, conditions, a );
		}
	}
	else if (type == HYPERBEL_LIN)
	{
		/* 1/y = (a2 + x) / a1 */
		/* convert observations back */
		for (unsigned int i = 0; i < N; i++)
		{
			observations[i] = 1. / observations[i];
		}

		double b0 = a[0];
		double b1 = a[1];
		a[0] = 1. / b1; /* a1 = 1/b2 */
		a[1] = a[0] * b0;  /* a2 = a1*b1 */
		fprintf( out, "\n#\n# corrected Parameters\n# " );
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "a%d=%.9f, ", j + 1, a[j] );
		}
		for (unsigned int i = 0; i < N; i++)
		{
			/* get calculated data points dependent on corrected
				 parameters */
			datac[i] = fhyperbel( i, conditions, a );
		}
	}
	else if (type == COORD_TRANSF)
	{
		/* print angle in degrees */

		fprintf( out, "  (%.4f degrees)", a[2] * 180. / M_PI );
	}
	else if (type == CIRCLE || type == CIRCLE_TLS)
	{
		/* evaluate result in terms of mean squared distance
		 * of points to circle
		 */
		double eval = 0, d, delta;

		cntInliers = 0;
		for (unsigned int i = 0; i < N; i++)
		{
			if (weights[i] > 0.)
			{
				//d = euclid_dist( (conditions[2 * i] - a[0]), (conditions[2 * i + 1] - a[1]) );
				d = hypot( (conditions[2 * i] - a[0]), (conditions[2 * i + 1] - a[1]) );
				/* orthogonal distance to curve of circle */
				delta = d - a[2];
				eval += delta * delta; /* sum up squared distances */
				cntInliers++;
			}
		}
		fprintf( out, "\n# mean squared distance to circle:  %.6f\n#",
			eval / cntInliers );
	}
	else if (type == CIRCLE_LIN)
	{
		double a1, a2, a3;
		double eval = 0, d, delta;
		/* convert parameters back	*/
		a1 = 0.5 * a[0];
		a2 = 0.5 * a[1];
		a3 = sqrt( a1 * a1 + a2 * a2 - a[2] );

		fprintf( out, "\n#\n# corrected Parameters\n# " );
		fprintf( out, "a1= %.*e", OP_DBL_Digs, a1 );
		fprintf( out, "a2= %.*e", OP_DBL_Digs, a2 );
		fprintf( out, "a3= %.*e", OP_DBL_Digs, a3 );
		fprintf( out, "\n#" );

		/* evaluate result in terms of mean squared distance
		 * of points to circle
		 */
		cntInliers = 0;
		for (unsigned int i = 0; i < N; i++)
		{
			if (weights[i] > 0.)
			{
				//d = euclid_dist( (conditions[2 * i] - a1), (conditions[2 * i + 1] - a2) );
				d = hypot( (conditions[2 * i] - a1), (conditions[2 * i + 1] - a2) );
				/* orthogonal distance to curve of circle */
				delta = d - a3;
				eval += delta * delta; /* sum up squared distances */
				cntInliers++;
			}
		}
		fprintf( out, "\n# mean squared distance to circle:  %.6f\n#",
			eval / cntInliers );
	}
	else if (type == NIST_BENNETT5_LOG)
	{
		/* convert observations back */
		for (unsigned int i = 0; i < N; i++)
		{
			observations[i] = -exp( observations[i] );
			datac[i] = -exp( datac[i] );
		}
		/* b1 = log(-a1)  => exp(b1)  = -a1 */
		/* b3 = 1/a3  => a3  = 1/b3 */
		a[0] = -exp( a[0] );
		a[2] = 1. / a[2];

		fprintf( out, "\n#\n# corrected Parameters\n# " );
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "a%d=%.9f, ", j + 1, a[j] );
		}
		fprintf( out, "\n#\n# uncertainties of corrected Parameters " );
		fprintf( out, "are not available yet\n# " );
	}

	/* ----------------------------------------- */
	fprintf( out, "\n#" );
	for (unsigned int d = 0; d < condDim; d++)
	{
		fprintf( out, formatStringText, "cond" );
		fprintf( out, "%d ", d + 1 );
	}
	fprintf( out, formatStringText, "observed" );
	fprintf( out, formatStringText, "fitted" );
	if (weightMode)
	{ /* skip output if no weights have been used */
		fprintf( out, formatStringText, "weight" );
		fprintf( out, formatStringText, "uncertainty" );
	}
	if (condDim > 1) 	fprintf( out, formatStringText, "diff." );
	else fprintf( out, formatStringText, "difference" );

	/*
	 * estimated weight should be w = 1/sigma^2
	 * thus estimated uncertainty is sigma = 1/ sqrt(weight)
	 * make nice output
	 */
	if (type == COORD_TRANSF)
	{
		fprintf( out, formatStringText, "observed2" );
		fprintf( out, formatStringText, "fitted2" );
		if (weightMode)
		{
			fprintf( out, formatStringText, "weight2" );
			fprintf( out, formatStringText, "uncert.2" );
		}
		fprintf( out, formatStringText, "diff.2" );
	}

	/*
	 * output of final results
	 */
	for (unsigned int i = 0; i < N * obsDim; i += obsDim)
	{
		double uncert;

		fprintf( out, "\n" );
		for (unsigned int d = 0; d < condDim; d++)
		{
			fprintf( out, formatString,
				conditions[condDim * (i / obsDim) + d] );
		}
		/* put both fitting results for fixed conditions in same row
		 */
		for (int d = 0; d < (int)obsDim; d++)
		{
			int k;
			k = i + d;
			if (weights[k] > 0)
				uncert = 1. / sqrt( weights[k] );
			else
				uncert = 9999.;
			fprintf( out, formatString, observations[k] );
			fprintf( out, formatString, datac[k] );
			if (weightMode)
			{
				fprintf( out, formatString, weights[k] );
				fprintf( out, formatString, uncert );
			}
			fprintf( out, formatStringP, (observations[k] - datac[k]) );
		}
		/* make a empty line to enforce plotting a mesh */
		if (type == LINEAR && M == 3 &&
			(conditions[condDim * i] != conditions[condDim * (i + 1)]))
		{
			/* plane approximation */
			fprintf( out, "\n" );
		}
		if ((type == QUAD_SURFACE || type == CLASSIF_2D) &&
			(conditions[condDim * i] != conditions[condDim * (i + 1)]))
		{
			/* surface approximation */
			fprintf( out, "\n" );
		}
	}
	fprintf( out, "\n#" );

#ifdef COMPARISON_WITH_CERTIFIED_VALUES
	int numEqualDigits;
	double* b = NULL;
	switch (type)
	{
	case LINEAR:
		if (M == 7) b = b_Longley;
		break;
	case POLYNOMIAL:
		if (M == 11) b = b_Filip;
		else if (M == 6)
		{
			if (strstr( inname, "1" ) != NULL)
				b = b_Wampler1;
			if (strstr( inname, "2" ) != NULL)
				b = b_Wampler2;
			if (strstr( inname, "3" ) != NULL)
				b = b_Wampler3;
			if (strstr( inname, "4" ) != NULL)
				b = b_Wampler4;
			if (strstr( inname, "5" ) != NULL)
				b = b_Wampler5;
		}
		break;
	case POLYNOMIAL_REG:
		if (M == 6)
		{
			if (strstr( inname, "1" ) != NULL)
				b = b_Wampler1;
		}
		break;
	case NIST_MGH09:
		b = b_MGH09;	break;
	case NIST_MGH10:
		b = b_MGH10;	break;
	case NIST_THURBER:
		b = b_Thurber;	break;
	case NIST_BOXBOD:
		b = b_BoxBOD;	break;
	case NIST_RAT42:
		b = b_Rat42;	break;
	case NIST_RAT43:
		b = b_Rat43;	break;
	case NIST_ECKERLE4:
		b = b_Eckerle4;	break;
	case NIST_BENNETT5:
	case NIST_BENNETT5_LOG:
		b = b_Bennett5;	break;
	case EXAMPLE:
		b = b_Example;	break;
	default:
		break;
	}
	if (b != NULL)
	{ /* reference values (certified parameters) are available*/
		numEqualDigits = getNumEqualDigits( a, b, M );
		for (unsigned int j = 0; j < M; j++)
		{
			fprintf( out, "\n# a%d= %+.*e", j + 1, OP_DBL_Digs, a[j] );
			fprintf( out, "\n# b%d= %+.*e", j + 1, OP_DBL_Digs, b[j] );
			int exponent_a = (int)floor( log10( fabs( a[j] ) ) );
			int exponent_b = (int)floor( log10( fabs( b[j] ) ) );
			int exponent_max = max( exponent_a, exponent_b );
			double aN = a[j] / pow( 10, exponent_max );
			double bN = b[j] / pow( 10, exponent_max );
			double diff = fabs( aN - bN );
			fprintf( out, "\t diff: %+.*lf", numEqualDigits + 2, diff );
			fprintf( out, "\t (diff: %+.*lf)\n", numEqualDigits, diff );
		}
		fprintf( out, "\n# number of equal parameter digits: %d\n", numEqualDigits );
		if (err)
		{
			fprintf( stdout, "&failed/%d/", err );
		}
		else
		{
			fprintf( stdout, "&\t%d", numEqualDigits );
		}
		if (!lsParams->linear)
		{
			fprintf( out, "\n#\n#   %d iterations", iterations );
			fprintf( stdout, " (%d)\t", iterations );
		}
	}
	//if (lsParams->linear && lsParams->expression_flag)
	//	fprintf( out, "\n#\n#   %d iterations", iterations );

#endif
	{
		time_t sekunden = time( NULL );
		struct tm* uhr = localtime( &sekunden );
		fprintf( out, "\n# %02d.%02d.%d", uhr->tm_mday, uhr->tm_mon + 1, uhr->tm_year + 1900 );
		fprintf( out, " %02d:%02d:%02d", uhr->tm_hour, uhr->tm_min + 1, uhr->tm_sec );
	}

	/* output of computer parameters*/
	CPUINFO ci;

	if (GetCPUInfo( &ci ))
	{
		fprintf( out, "\n# " );
		fprintf( out, "Vendor: %s, ", ci.vendor );
		fprintf( out, "Model: %s, ", ci.brand );
		fprintf( out, "Family: %u, ", ci.family );
		fprintf( out, "Model ID: %u, ", ci.model );
		fprintf( out, "Stepping: %u, ", ci.stepping );
		fprintf( out, "\n# " );
		fprintf( out, "Physical cores: %u, ", ci.physicalCores );
		fprintf( out, "Logical cores: %u, ", ci.logicalProcessors );
		fprintf( out, "CPU Speed: %.0f MHz, ", ci.cpuSpeedMHz );
		fprintf( out, "Available Memory: %.2f GB\n",
			ci.availableMemoryBytes / (1024.0 * 1024.0 * 1024.0) );
	}

endfunc:
	fprintf( stderr, "\n" );
	fflush( stdout );

	free_vector( &weights );
	free_vector( &weightsOld );
	free_vector( &observations );
	free_vector( &datac );
	free_vector( &conditions );
	free_vector( &deviates );
	free_vector( &deviatesAbsolute );
	//free_vector( &deltasq );
	free_matrix( &jacob );
	free_matrix( &covar );

	/* frree tinyexpr handle */
	if (teExpr != NULL)  te_free( teExpr );

	if (argc > 2)
	{
		if (out != NULL)
			fclose( out );
	}
	if (in != NULL)
		fclose( in );

	if (err)
	{
		fprintf( stderr, "\n failed.\n" );
		return err;
	}
	else
	{
		fprintf( stderr, "\n ready.\n" );
		return 0;
	}
}

/*---------------------------------------------------------------
 *	is_data_line()
 *--------------------------------------------------------------*/
int
is_data_line( char* line, unsigned int N )
{
	unsigned int idx;
	/* scan leading white spaces	*/
	for (idx = 0; idx < N; idx++)
	{
		if (line[idx] != ' '			/* space	*/
			&& line[idx] != '\t'	/* tab		*/
			)
			break;
	}
	/* check next character	*/
	if (line[idx] == '#')		/* comment */
	{
		return 0;
	}
	if (line[idx] == '\n')	/* newline */
	{
		return 0;
	}
	if (line[idx] == '\r')	/* carrige return (Windows)  */
	{
		return 0;
	}

	return 1;	/* is data line */
}

/*---------------------------------------------------------------
 *	get_nth_field()
 * scans a string up to the desired column (field)
 *--------------------------------------------------------------*/
char*
get_nth_field( char* line, int n )
{
	char* ptr = NULL;
	int ch, i = 0, loop_flag, cnt, field_flag;

	if (line == NULL)
		return ptr;

	loop_flag = 1;
	field_flag = 0;
	cnt = 0;
	do
	{
		ch = line[i];
		if (ch == '\0')
		{
			loop_flag = 0;
			break;
		}
		if (!field_flag)
		{
			if (ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n')
			{
				field_flag = 1;
				cnt++;
				if (cnt == n)
				{
					/* desired field is found */
					ptr = &(line[i]);
					loop_flag = 0;
				}
			}
		}
		else
		{
			/* search next white space */
			if (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n')
			{
				field_flag = 0;
			}
		}
		i++;
	} while (loop_flag);
	return ptr;
}
