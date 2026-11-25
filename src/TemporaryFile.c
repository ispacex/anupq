/****************************************************************************
**
*A  TemporaryFile.c             ANUPQ source                   Eamonn O'Brien
**
*Y  Copyright 1995-2001,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
*Y  Copyright 1995-2001,  School of Mathematical Sciences, ANU,     Australia
**
*/

#include "pq_defs.h"
#include "pq_functions.h"
#include "constants.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* set up a temporary file and return an appropriate FILE * indicator;
   if in Unix environment, open temporary file in directory specified
   by value of environment variable TMPDIR, else on /var/tmp */

FILE *TemporaryFile(void)
{
   FILE *file = NULL;
   const char *dir = getenv("TMPDIR");

   /* glibc's tmpfile() ignores TMPDIR and puts its files in /tmp, often a
      small tmpfs (MATH-525): honour TMPDIR, the file being unlinked at once
      as tmpfile()'s is */
   if (dir != NULL && *dir != '\0' && strlen(dir) < 4000) {
      char name[4096];
      int fd;
      sprintf(name, "%s/pqtmpXXXXXX", dir);
      fd = mkstemp(name);
      if (fd != -1) {
         unlink(name);
         file = fdopen(fd, "w+");
         if (file == NULL)
            close(fd);
      }
   }
   if (file == NULL)
      file = tmpfile();

   if (file == NULL) {
      perror("Cannot open temporary file");
      exit(FAILURE);
   }

   return file;
}
