/****************************************************************
 *
 * File........:	get_option.c
 * Function....:	reading and analysing of
 *                     command-line parameters/options
 * Author......:	Tilo Strutz
 * last changes:	15.08.2006, 28.07.2026
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as 
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *
 ****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "get_option.h"
#include "errmsg.h"

/* contains argument of option, if OptArg != NULL	*/
char *OptArg = NULL;
char CheckStr[256];	
/* CheckStr[] will be initialised with NEEDEDOPTIONS
 * all used optionen are deleted; if at least one option remains, 
 * an error message is output */
char *optstr;
int opt_num = 1;

/*---------------------------------------------------------------
 *	check_opt()
 *---------------------------------------------------------------*/
int
check_opt( const char *name)
{
	char *ptr;
	int err = 0;
	unsigned int i, len;

	len = strlen( CheckStr);
	for (i = 0; i < len; i++)
	{
		if (( CheckStr[i] != ':') && (CheckStr[i] != ' '))
		{
			ptr = &CheckStr[i];
			ptr = (char*)strpbrk( ptr, ";:");
			ptr[0] = '\0';
			err = ERR_CALL;
			break;
		}
	}
	if (err)
	{
		fprintf( stderr, "\n Missing Option for (-%s)!", &CheckStr[i]);
		usage( name);
	}
	return err;
}

/*---------------------------------------------------------------
 *	get_option()
 *
 *	opt_num is number of option to be read
 *	result: option string 
 * 	required: global string containing all options
 * 	at first call opt_num must be equal to 1 !
 *
 *---------------------------------------------------------------*/
char *
get_option( int argc, const char *argv[])
{
	char optstring[256], *ptr, c, d, string[256];
	char *gerrstr="#";
	int len, i, num;

	if (opt_num == 1)
	{
		strcpy( CheckStr, NEEDEDOPTIONS);
	}
	if (opt_num > (argc - 1))
	{
		return (NULL);
	}
	else if (argv[opt_num][0] == '+')
	{
		/* + signals end of parameter list */
		opt_num++;
		return (NULL);
	}
	else if (argv[opt_num][0] != '-')
	{
		fprintf( stderr, "\n Option-Error !!  *********** ");
		fprintf( stderr,
			"\n every Option must start with '-' (%s)!", argv[opt_num]);
		/* usage( argv[0]); */
		return gerrstr;
	}

	/***** copy without '-'	*****/
	num = opt_num;
	strcpy( optstring, argv[num]);
	strcpy( string, &optstring[1]);

	len = (int)strlen( string);
	if (len == 0)
	{
		/*  single '-'	 */
		fprintf( stderr, "\n Option-Error !!  *********** ");
		fprintf( stderr, "\n lonely dash !");
		/* usage( argv[0]); */
		return gerrstr;
	}
	ptr = OPTIONSTRING;
	do
	{
		/* search option string in OPTIONSTRING */
		ptr = (char*)strstr( ptr, string);
		if (ptr == NULL)
		{
			fprintf( stderr, "\n Option-Error !!  *********** ");
			fprintf( stderr, "\n Unknown Option (%s)!", optstring);
			/* usage( argv[0]); */
			return gerrstr;
		}
		c = ptr[len];	/* remember subsequent character */
		d = ptr[-1];	/* remember predecessor */

		/* skip this entry by searching for next ':' or ';' */
		ptr = (char*)strpbrk( ptr, ";:.");
	} while (( (c != ':') && (c != ';') && (c != '.')) || (( d != ':')
		&& (d != ';')&& (d != '.')));

	if (c == ';')	/* info, whether argument follows */
	{
		OptArg = NULL;
		opt_num++;
	}
	else
	{
		opt_num++;
		if (opt_num > (argc - 1))
		{
			fprintf( stderr, "\n Option-Error !!  *********** ");
			fprintf( stderr, "\n Missing Argument for (%s)!", optstring);
			/* usage( argv[0]); */
			return gerrstr;
		}
		else if (argv[opt_num][0] == '-'  && c == ':')
		{
			fprintf( stderr, "\n Option-Error !!  *********** ");
			fprintf( stderr, "\n Missing Argument for (%s)!", optstring);
			/* usage( argv[0]); */
			return gerrstr;
		}
		else 
		{
				/* if c == '.' then negativ parameter are allowed */
		}
		OptArg = (char*)argv[opt_num];
		opt_num++;
	}
	strcpy( string, ":");
	strcat( string, &optstring[1]);
	strcat( string, ":");
	ptr = (char*)strstr( CheckStr, string);
	if (ptr != NULL)
		for (i = 0; i < len; i++)
			ptr[i + 1] = ' ';

	return (( char*)argv[num]);
}
