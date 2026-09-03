/*****************************************************************
 *
 * File....:	cpuID.h
 * Function:	gets information about cpu and memory
 * Author..:	Tilo Strutz
 * Date....:	27.07.2026
 *
 * changes:
 *
 * LICENCE DETAILS: see software manual
 *	free academic use
 *  cite source as
 * "Strutz, T.: Data Fitting and Uncertainty. 3rd edition,
 *   Springer Nature, 2026"
 *****************************************************************/
#ifndef CPUINFO_H
#define CPUINFO_H

#include <stdint.h>

typedef struct
{
  char vendor[32];
  char brand[128];

  unsigned int family;
  unsigned int model;
  unsigned int stepping;

  unsigned int logicalProcessors;
  unsigned int physicalCores;

  double cpuSpeedMHz;

  uint64_t availableMemoryBytes;

} CPUINFO;

int GetCPUInfo( CPUINFO* info );

#endif
