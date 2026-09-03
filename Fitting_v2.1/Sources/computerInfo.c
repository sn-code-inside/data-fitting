/*****************************************************************
 *
 * File....:	computerInfo.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cpuInfo.h"

#ifdef _WIN32
  #include <windows.h>
  #include <winreg.h>
  #include <intrin.h>

  static void cpuidex(
    int cpuinfo[4],
    int function_id,
    int subfunction_id )
  {
    __cpuidex( cpuinfo, function_id, subfunction_id );
  }
#else // Linux
  #include <cpuid.h>
  #include <unistd.h>
  #include <sys/sysinfo.h>
  #include <dirent.h>


  static void cpuidex(
    int cpuinfo[4],
    int function_id,
    int subfunction_id )
  {
    __cpuid_count(
      function_id,
      subfunction_id,
      cpuinfo[0],
      cpuinfo[1],
      cpuinfo[2],
      cpuinfo[3] );
  }
#endif

#include <string.h>
#include <ctype.h>

/*---------------------------------------------------
 * TrimTrailingSpaces()
 *---------------------------------------------------*/
static void TrimTrailingSpaces( char* str )
{
  size_t len;

  if (str == NULL)  return;

  len = strlen( str );
  while (len > 0)
  {
    unsigned char c = (unsigned char)str[len - 1];

    if (!isspace( c )) break;

    str[len - 1] = '\0';
    --len;
  }
}

 /*---------------------------------------------------
  * GetVendor()
  *---------------------------------------------------*/
  static void GetVendor( char* vendor, size_t vendorSize )
  {
    int info[4];
    if (vendor == NULL || vendorSize < 13)
      return;

    cpuidex( info, 0, 0 );

    memcpy( vendor + 0, &info[1], 4 );
    memcpy( vendor + 4, &info[3], 4 );
    memcpy( vendor + 8, &info[2], 4 );

    vendor[12] = 0;
  }

  /*---------------------------------------------------
   * GetBrandString()
   *---------------------------------------------------*/
  static void GetBrandString( char* brand )
  {
    int info[4];

    memset( brand, 0, 128 );
    cpuidex( info, 0x80000000, 0 );
    if ((unsigned int)info[0] < 0x80000004)
      return;

    cpuidex( (int*)(brand + 0), 0x80000002, 0 );
    cpuidex( (int*)(brand + 16), 0x80000003, 0 );
    cpuidex( (int*)(brand + 32), 0x80000004, 0 );
  }

  /*---------------------------------------------------
  * GetFamilyModelStepping()
  *---------------------------------------------------*/
  static void GetFamilyModelStepping(
    unsigned int* family,
    unsigned int* model,
    unsigned int* stepping )
  {
    int info[4];

    cpuidex( info, 1, 0 );

    unsigned int eax = (unsigned int)info[0];
    unsigned int baseFamily =   (eax >> 8) & 0xF;
    unsigned int baseModel =    (eax >> 4) & 0xF;
    unsigned int extFamily =    (eax >> 20) & 0xFF;
    unsigned int extModel =     (eax >> 16) & 0xF;
    *stepping = eax & 0xF;

    if (baseFamily == 0xF)
      *family = baseFamily + extFamily;
    else
      *family = baseFamily;

    if (baseFamily == 0x6 || baseFamily == 0xF)
      *model = baseModel + (extModel << 4);
    else
      *model = baseModel;
  }

#ifdef _WIN32
  /*---------------------------------------------------
  * GetCpuSpeedMHz()
  *---------------------------------------------------*/
  static double GetCpuSpeedMHz( void )
  {
    HKEY hKey;
    DWORD data;
    DWORD dataSize = sizeof( data );

    if (RegOpenKeyExA(
      HKEY_LOCAL_MACHINE,
      "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0",
      0,  KEY_READ,
      &hKey ) == ERROR_SUCCESS)
    {
      if (RegQueryValueExA(
        hKey, "~MHz",
        NULL,  NULL,
        (LPBYTE)&data,
        &dataSize ) == ERROR_SUCCESS)
      {
        RegCloseKey( hKey );
        return (double)data;
      }
      RegCloseKey( hKey );
    }
    return 0.0;
  }
#else
  static double GetCpuSpeedMHz( void )
  {
    FILE* f;

    char line[256];
    double mhz = 0;

    f = fopen( "/proc/cpuinfo", "r" );

    if (!f)   return 0;

    while (fgets( line, sizeof( line ), f ))
    {
      if (sscanf( line, "cpu MHz : %lf", &mhz ) == 1)
        break;
    }
    fclose( f );
    return mhz;
  }
#endif

/*---------------------------------------------------
 * GetAvailableMemory()
 *---------------------------------------------------*/
static uint64_t GetAvailableMemory( void )
{
#ifdef _WIN32

  MEMORYSTATUSEX ms;

  memset( &ms, 0, sizeof( ms ) );
  ms.dwLength = sizeof( ms );

  if (GlobalMemoryStatusEx( &ms ))
  {
    return (uint64_t)ms.ullAvailPhys;
  }
  return 0;
#else
    FILE* f;
    char line[256];
    uint64_t memAvailableKB = 0;

    f = fopen( "/proc/meminfo", "r" );
    if (f == NULL)
    {
      return 0;
    }

    while (fgets( line, sizeof( line ), f ) != NULL)
    {
      if (sscanf( line, "MemAvailable: %llu kB",
        &memAvailableKB ) == 1)
      {
        fclose( f );
        return memAvailableKB * 1024ULL;
      }
    }
    fclose( f );
    return 0;
#endif
}

/*---------------------------------------------------
 * GetPhysicalCoreCount()
 *---------------------------------------------------*/
static unsigned int GetPhysicalCoreCount( void )
#ifdef _WIN32
{
  DWORD length = 0;
  PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX buffer;
  PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX ptr;
  unsigned int coreCount = 0;

  GetLogicalProcessorInformationEx(
    RelationProcessorCore,  NULL,   &length );

  buffer =
    (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)
    malloc( length );

  if (buffer == NULL)
  {
    return 0;
  }

  if (!GetLogicalProcessorInformationEx(
    RelationProcessorCore,   buffer,   &length ))
  {
    free( buffer );
    return 0;
  }

  ptr = buffer;
  while ((char*)ptr < ((char*)buffer + length))
  {
    if (ptr->Relationship == RelationProcessorCore)
    {
      ++coreCount;
    }
    ptr =
      (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)
      ((char*)ptr + ptr->Size);
  }
  free( buffer );

  return coreCount;
}
#else
{
  unsigned int coreCount = 0;
  FILE* f;
  char path[256];
  char line[64];

  int seenPhysical[1024];
  int seenCore[1024];
  int used = 0;

  memset( seenPhysical, -1, sizeof( seenPhysical ) );
  memset( seenCore, -1, sizeof( seenCore ) );

  for (int cpu = 0; cpu < 1024; ++cpu)
  {
    int physicalId = -1;
    int coreId = -1;

    snprintf(  path,     sizeof( path ),
      "/sys/devices/system/cpu/cpu%d/topology/physical_package_id",
      cpu );

    f = fopen( path, "r" );
    if (f == NULL)
    {
      continue;
    }

    if (fgets( line, sizeof( line ), f ))
    {
      physicalId = atoi( line );
    }
        fclose( f );

    snprintf(
      path,
      sizeof( path ),
      "/sys/devices/system/cpu/cpu%d/topology/core_id",
      cpu );

    f = fopen( path, "r" );
    if (f == NULL)
    {
      continue;
    }

    if (fgets( line, sizeof( line ), f ))
    {
      coreId = atoi( line );
    }
    fclose( f );

    int found = 0;
    for (int i = 0; i < used; ++i)
    {
      if (seenPhysical[i] == physicalId &&
        seenCore[i] == coreId)
      {
        found = 1;        break;
      }
    }

    if (!found)
    {
      seenPhysical[used] = physicalId;
      seenCore[used] = coreId;
      ++used;
      ++coreCount;
    }
  }
  return coreCount;
}
#endif


/*---------------------------------------------------
 * GetCPUInfo()
 *---------------------------------------------------*/
int GetCPUInfo( CPUINFO* info )
{
  if (info == NULL)  return 0;

  memset( info, 0, sizeof( CPUINFO ) );

  /*
   * CPU vendor
   * Example:
   *   GenuineIntel,  AuthenticAMD
   */
  GetVendor( info->vendor, sizeof( info->vendor ) );

  /*
   * Full processor name
   * Example:   Intel(R) Core(TM) i7-13700K
   */
  GetBrandString( info->brand );
  TrimTrailingSpaces( info->brand );
  /*
   * Family / Model / Stepping from CPUID
   */
  GetFamilyModelStepping(
    &info->family,  &info->model,  &info->stepping );

  /*
   * Logical processor count
   */
#ifdef _WIN32
  {
    SYSTEM_INFO systemInfo;

    memset( &systemInfo, 0, sizeof( systemInfo ) );
    GetNativeSystemInfo( &systemInfo );
    info->logicalProcessors =
      (unsigned int)systemInfo.dwNumberOfProcessors;
  }

#else
  {
    long cpuCount;

    cpuCount = sysconf( _SC_NPROCESSORS_ONLN );
    if (cpuCount > 0)
    {
      info->logicalProcessors =
        (unsigned int)cpuCount;
    }
  }
#endif
  info->physicalCores = GetPhysicalCoreCount();

  /*
   * Nominal CPU frequency
   */
  info->cpuSpeedMHz = GetCpuSpeedMHz();

  /*
   * Available physical memory
   */
  info->availableMemoryBytes =
    GetAvailableMemory();

  return 1;
}