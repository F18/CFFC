//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BPRESOURCE_H_
#define _BPRESOURCE_H_

#define COMMENT_CHARS "%"

int BpGetInt(const char *resource_name, int val);
double BpGetDouble(const char *resource_name, double val);
char *BpGetStr(const char *resource_name, char *val);

#endif // _BPRESOURCE_H_
