/****************************************************************************
**
*A  GAP_library_file.c          ANUPQ source
**
**  The descendants pq writes for GAP, optionally gzip-compressed.
**
**  With the environment variable PQ_GZ set to 1, the group descriptions and
**  automorphisms that setup_reps writes for GAP go, instead of the file
**  GAP_library (opened in mode "a+" for every item), into a pipe to one
**  "gzip -1" process appending to GAP_library.gz in the same directory.  The
**  pipe stays open while the descendants are constructed and is closed (gzip
**  ends its member and exits; pclose waits for it) when the construction menu
**  option returns, i.e. before pq prompts GAP again, and at exit.  So the file
**  is a complete gzip file at the point GAP reads it; concatenated members
**  (several constructions in one pq process) are again a gzip file, which GAP
**  reads with InputTextFile on the .gz name.  Decompressed, it is byte for
**  byte what the descendants part of GAP_library would have been.  The header
**  line (ANUPQmagic) and the standard presentation output stay in GAP_library.
**
**  Without PQ_GZ nothing changes: OpenGAPLibrary is OpenFile(name, "a+") and
**  CloseGAPLibrary is CloseFile.
**
**  On the heavy vertices of the dimension-subrings 2^12 search GAP_library is
**  tens of GB of text (MATH-525).
*/

#include "pq_defs.h"
#include "constants.h"
#include "pq_functions.h"

#include <fcntl.h>
#include <unistd.h>
#include <sys/wait.h>

static int gz_mode = -1;
static FILE *gz_pipe = NULL;
static char *gz_name = NULL;

int GAPLibraryCompressed(void)
{
   if (gz_mode < 0) {
      const char *e = getenv("PQ_GZ");
      gz_mode = (e != NULL && strcmp(e, "1") == 0);
   }
   return gz_mode;
}

static char *gz_file_name(const char *name)
{
   char *s = malloc(strlen(name) + 4);
   if (s == NULL) {
      perror(NULL);
      exit(FAILURE);
   }
   sprintf(s, "%s.gz", name);
   return s;
}

/* remove a stale compressed library (GAP empties GAP_library at the start
   of each pq process in the same way) */

void RemoveGAPLibraryCompressed(const char *name)
{
   char *s;
   if (!GAPLibraryCompressed())
      return;
   s = gz_file_name(name);
   unlink(s);
   free(s);
}

static int finish(void)
{
   int r = 0;
   if (gz_pipe != NULL) {
      int status = pclose(gz_pipe);
      gz_pipe = NULL;
      if (status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
         fprintf(stderr, "pq: gzip -1 >> %s failed (status %d)\n", gz_name,
                 status);
         r = -1;
      }
   }
   return r;
}

static void finish_at_exit(void)
{
   (void)finish();
}

/* open the GAP library file name for appending one item */

FILE *OpenGAPLibrary(const char *name)
{
   char *s, *cmd;
   const char *q;

   if (!GAPLibraryCompressed())
      return OpenFile(name, "a+");

   s = gz_file_name(name);
   if (gz_pipe != NULL && strcmp(s, gz_name) == 0) {
      free(s);
      return gz_pipe;
   }
   FinishGAPLibrary();
   free(gz_name);
   gz_name = s;

   for (q = gz_name; *q; ++q)
      if (*q == '\'') {
         fprintf(stderr, "pq: PQ_GZ: quote in file name %s\n", gz_name);
         exit(FAILURE);
      }
   cmd = malloc(strlen(gz_name) + 32);
   if (cmd == NULL) {
      perror(NULL);
      exit(FAILURE);
   }
   sprintf(cmd, "exec gzip -1 >> '%s'", gz_name);
   fflush(NULL);
   gz_pipe = popen(cmd, "w");
   free(cmd);
   if (gz_pipe == NULL) {
      perror("pq: PQ_GZ: popen gzip");
      exit(FAILURE);
   }
   /* not inherited by the processes pq starts later (GAP for stabilisers) */
   fcntl(fileno(gz_pipe), F_SETFD, FD_CLOEXEC);
   {
      static int registered = 0;
      if (!registered) {
         atexit(finish_at_exit);
         registered = 1;
      }
   }
   return gz_pipe;
}

/* close the item opened by OpenGAPLibrary */

void CloseGAPLibrary(FILE *fp)
{
   if (!GAPLibraryCompressed()) {
      CloseFile(fp);
      return;
   }
   if (ferror(fp)) {
      fprintf(stderr, "pq: PQ_GZ: write to gzip -1 >> %s failed\n", gz_name);
      exit(FAILURE);
   }
}

/* end the gzip member written so far (a no-op without PQ_GZ): called when
   a construction of descendants is complete, before GAP reads the file */

void FinishGAPLibrary(void)
{
   if (finish() != 0)
      exit(FAILURE);
}
