
/*
 * Including libraries, for avoiding repetitive include.
 */

#ifndef __INCLUDE_PP
#define __INCLUDE_PP


#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <map>
#include <random>


#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>


#include <omp.h>
#include <pthread.h>


#ifdef __linux__
    #define      __stderrp stderr
    #define      __stdoutp stdout
    #define      __stdinp  stdin
#endif


/*
 * Default file access permissions for new files.
 */
#define	FILE_MODE	(S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

#endif /* __INCLUDE_PP */