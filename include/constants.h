#ifndef CONSTANTS_H
#define CONSTANTS_H

#define PI 3.141592654
#define FOUR_PI 12.566370614

#ifdef _WIN32
    #define PATH_SEPARATOR '\\'
    #define OTHER_SEPARATOR '/'
#else
    #define PATH_SEPARATOR '/'
    #define OTHER_SEPARATOR '\\'
#endif

#endif