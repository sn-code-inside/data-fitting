/*****************************************************************
 *
 * File........:	errmsg.h
 * Function....:	error messages
 * Author......:	Tilo Strutz
 * last changes:	25.01.2010, 1.4.2011, 28.07.2026
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 *****************************************************************/
#ifndef ERRMSG_H
#define ERRMSG_H

#define ERR_IS_SINGULAR	1
#define ERR_IS_SINGULAR_MSG	\
	"\n### %s: Matrix %s is singular\n"
#define ERR_HIGH_CONDITIONNUMBER	2
#define ERR_HIGH_CONDITIONNUMBER_MSG	\
	"\n### %s: matrix ill-conditioned (%s)\n"
#define ERR_IS_ZERO	3
#define ERR_IS_ZERO_MSG	\
	"\n\n### %s: Value for %s is zero  ###\n"
#define ERR_CALL	4
#define ERR_CALL_MSG \
	"\n### %s: Wrong command-line parameters. \n  %s\n"
#define ERR_OPEN_READ	5
#define ERR_OPEN_READ_MSG	\
	"\n### %s: Cannot open %s for reading\n"
#define ERR_OPEN_WRITE	6
#define ERR_OPEN_WRITE_MSG	\
	"\n### %s: Cannot open %s for writing\n"
#define ERR_ALLOCATE		7
#define ERR_ALLOCATE_MSG	\
	"\n### %s: Cannot allocate %s\n"
#define ERR_NOT_DEFINED	8
#define ERR_NOT_DEFINED_MSG	\
	"\n### %s: Value %d for %s is not defined\n"
#define ERR_SVD_NO_CONVERGENCE 9
#define ERR_SVD_NO_CONVERGENCE_MSG \
	"\n### %s: SVD did not converge\n"
#define ERR_SVD_STRANGE_INPUT 10
#define ERR_SVD_STRANGE_INPUT_MSG \
	"\n### %s: SVD got strange input\n"
#define ERR_IS_INFINITE	11
#define ERR_IS_INFINITE_MSG	\
	"\n### %s: Variable %s is infinite\n"
#define ERR_IS_NAN	12
#define ERR_IS_NAN_MSG	\
	"\n### %s: Variable %s is not a number (NAN)\n"
#define ERR_TOO_MANY_PARAMETERS	13
#define ERR_TOO_MANY_PARAMETERS_MSG	\
	"\n### %s: too many parameters (%d), %s\n"
#define ERR_COLUMN_DOES_NOT_EXIST	14
#define ERR_COLUMN_DOES_NOT_EXIST_MSG	\
	"\n### %s: column (%d) does not exist (%s)\n"
#define ERR_OUT_OF_RANGE	15
#define ERR_OUT_OF_RANGE_MSG	\
	"\n### %s: value is out of range\n"
#define ERR_DATA_PROBLEMS	16
#define ERR_DATA_PROBLEMS_MSG	\
	"\n### %s: input data are suspicious\n"
#define ERR_ANY	99
#define ERR_ANY_MSG	\
	"\n### %s: Variable %s is not defined \n"

int errmsg( int err, char *rtn, char *text, int value);

#endif
