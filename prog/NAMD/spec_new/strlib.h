/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   strlib contains a number of useful routines for doing file I/O
   and some basic string manipulation.  These routines are used for 
   reading in the parameter and .psf files
*/

#ifndef STRLIB_H

#define STRLIB_H

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <vector>

#ifdef WIN32
#define strcasecmp(s,t) stricmp(s,t)
#define strncasecmp(s,t,n) strnicmp(s,t,n)
#else
#include <strings.h>
#endif

#include "common.h"

void	NAMD_truncate(char *);		//  Remove trailing spaces from
					//  a string
int	NAMD_read_line(FILE *, char *, int bufsize=512); //  Read in a line from a file
int	NAMD_find_word(const char *, const char *); //  Check for given word in a
					//  string
int	NAMD_blank_string(char *);	//  Check to see if a string
					//  is blank
void	NAMD_find_first_word(char *, char *);
					//  Find the first word in a string
int	NAMD_read_int(FILE *, const char *);  //  Read an integer from a file
void	NAMD_pad(char *, size_t);	//  Pad a string with leading spaces
void	NAMD_remove_comment(char *);	//  Remove comments at the end of
					//  a line demarked by !


// Remove leading and trailing spaces from a string
inline std::string NAMD_strip(const std::string str) {
	// skip leading spaces
	const char *s = str.c_str();
	for ( ; *s != STRINGNULL && isspace(*s); s++ )
	  ;
	// remove trailing spaces
	int i;
	for ( i = strlen(s); i >= 0 && isspace(s[i-1]); i-- )
	  ;
	return std::string(s, i);
}

// Split a string into an array of strings
inline std::vector<std::string> NAMD_splitstr(const std::string s, char delim) {
	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> arr;
	while ( std::getline(ss, item, delim) )
	  arr.push_back(item);
	return arr;
}

#endif

