//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "BpResource.h"

#define RESOURCES_FILENAME ".BpResource"

int BpGetInt(const char *resource_name, int val)
{
    FILE *fp;
    char line[82];
    char *line_ptr;
    int length;

    // open resource file
    fp = fopen(RESOURCES_FILENAME, "r");
    if (fp == NULL)
        return val;

    while (1)
    {
        line_ptr = fgets(line, 82, fp);
        if (line_ptr == NULL)
	{
	    fclose(fp);
            return val;
	}

	// length up to colon
        length = strcspn(line_ptr, ":");

	// names must match in length
	if (length != strlen(resource_name))
	    continue;

        line_ptr[length] = '\0';
        if (strcmp(line_ptr, resource_name) == 0)
	{
	    fclose(fp);
	    line_ptr = strtok(line_ptr+length+1, COMMENT_CHARS);
            return atoi(line_ptr);
	}
    }
}

double BpGetDouble(const char *resource_name, double val)
{
    FILE *fp;
    char line[82];
    char *line_ptr;
    int length;

    // open resource file
    fp = fopen(RESOURCES_FILENAME, "r");
    if (fp == NULL)
        return val;

    while (1)
    {
        line_ptr = fgets(line, 82, fp);
        if (line_ptr == NULL)
	{
	    fclose(fp);
            return val;
	}

	// length up to colon
        length = strcspn(line_ptr, ":");

	// names must match in length
	if (length != strlen(resource_name))
	    continue;

        line_ptr[length] = '\0';
        if (strcmp(line_ptr, resource_name) == 0)
	{
	    fclose(fp);
	    line_ptr = strtok(line_ptr+length+1, COMMENT_CHARS);
	    return atof(line_ptr);
	}
    }
}

// Note: calling function must free the string that is returned.
//
char *BpGetStr(const char *resource_name, char *val)
{
    FILE *fp;
    char line[82];
    char *line_ptr;
    int length;

    // open resource file
    fp = fopen(RESOURCES_FILENAME, "r");
    if (fp == NULL)
        return val;

    while (1)
    {
        line_ptr = fgets(line, 82, fp);
        if (line_ptr == NULL)
	{
	    fclose(fp);
            return val;
	}

	// length up to colon
        length = strcspn(line_ptr, ":");

	// names must match in length
	if (length != strlen(resource_name))
	    continue;

        line_ptr[length] = '\0';
        if (strcmp(line_ptr, resource_name) == 0)
	{
	    line_ptr = strtok(line_ptr+length+1, COMMENT_CHARS);

	    fclose(fp);
	    length = strlen(line_ptr); // length of resource value
	    return strcpy(new char[length+1], line_ptr); // need to free
	}
    }
}
