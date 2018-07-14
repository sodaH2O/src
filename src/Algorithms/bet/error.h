/* error.h
   error handling definition

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   07-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: error.h 91 2007-05-06 17:47:45Z bakker $
*/


#ifndef _ERROR_DEF
#define _ERROR_DEF

#include <stdio.h>

void initErrorMsg(int level = 0, const char *fbase = NULL);
void initErrorFilename(const char *fbase);

void writeError();

void setReportLevel(int level);

#endif
