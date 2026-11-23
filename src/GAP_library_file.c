/****************************************************************************
**
*A  GAP_library_file.c          ANUPQ source
**
**  The descendants pq writes for GAP, optionally gzip-compressed.
**
**  With the environment variable PQ_GZ set to 1, the group descriptions and
**  automorphisms that setup_reps writes for GAP go, instead of the file
**  GAP_library (opened in mode "a+" for every descendant), to GAP_library.gz
**  in the same directory: each write is collected in memory and appended to
**  one zlib stream (level 1) that stays open while the descendants are
**  constructed and is closed (gzip trailer written) when the construction
**  menu option returns, and at exit.  The file is a valid gzip file at the
**  point GAP reads it, and concatenated members (several constructions in one
**  pq process) are again a valid gzip file; GAP reads it with InputTextFile on
**  the .gz name.  Its content, decompressed, is byte for byte what the
**  descendants part of GAP_library would have been.  The header line
**  (ANUPQmagic) and the standard presentation output stay in GAP_library.
**
**  Without PQ_GZ nothing changes: OpenGAPLibrary is OpenFile(name, "a+").
**
**  On a heavy vertex (MATH-525) GAP_library is tens of GB of text that
**  compresses about tenfold.
*/

#include "pq_defs.h"
#include "constants.h"
#include "pq_functions.h"

#if defined(HAVE_ZLIB)
#include <zlib.h>
#endif

static int gz_mode = -1;

int GAPLibraryCompressed(void)
{
   if (gz_mode < 0) {
      const char *e = getenv("PQ_GZ");
      gz_mode = (e != NULL && strcmp(e, "1") == 0);
#if !defined(HAVE_ZLIB)
      if (gz_mode) {
         fprintf(stderr, "PQ_GZ=1 but pq was built without zlib\n");
         exit(FAILURE);
      }
#endif
   }
   return gz_mode;
}

#if defined(HAVE_ZLIB)

static gzFile gz_out = NULL;
static char *gz_name = NULL;
static char *gz_buf = NULL;
static size_t gz_len = 0;

static void finish_at_exit(void);

static void gz_fail(const char *what)
{
   int err = 0;
   const char *msg = gz_out != NULL ? gzerror(gz_out, &err) : "";
   fprintf(stderr, "pq: %s %s failed: %s\n", what, gz_name ? gz_name : "", msg);
   perror(NULL);
   exit(FAILURE);
}

/* remove a stale compressed library (GAP empties GAP_library at the start
   of each pq process in the same way) */

void RemoveGAPLibraryCompressed(const char *name)
{
   char *s;
   if (!GAPLibraryCompressed())
      return;
   s = malloc(strlen(name) + 4);
   if (s == NULL) {
      perror(NULL);
      exit(FAILURE);
   }
   sprintf(s, "%s.gz", name);
   unlink(s);
   free(s);
}

#else

void RemoveGAPLibraryCompressed(const char *name)
{
   (void)GAPLibraryCompressed();
}

#endif

/* open the GAP library file name for appending one item */

FILE *OpenGAPLibrary(const char *name)
{
   if (!GAPLibraryCompressed())
      return OpenFile(name, "a+");

#if defined(HAVE_ZLIB)
   if (gz_name == NULL || strlen(gz_name) != strlen(name) + 3 ||
       strncmp(gz_name, name, strlen(name)) != 0) {
      FinishGAPLibrary();
      free(gz_name);
      gz_name = malloc(strlen(name) + 4);
      if (gz_name == NULL) {
         perror(NULL);
         exit(FAILURE);
      }
      sprintf(gz_name, "%s.gz", name);
   }
   {
      FILE *fp = open_memstream(&gz_buf, &gz_len);
      if (fp == NULL) {
         perror("open_memstream");
         exit(FAILURE);
      }
      return fp;
   }
#else
   return NULL;
#endif
}

/* close the item opened by OpenGAPLibrary: append it to the library */

void CloseGAPLibrary(FILE *fp)
{
   if (!GAPLibraryCompressed()) {
      CloseFile(fp);
      return;
   }

#if defined(HAVE_ZLIB)
   {
      size_t done = 0;
      CloseFile(fp); /* fixes gz_buf and gz_len */
      if (gz_len > 0) {
         if (gz_out == NULL) {
            static int registered = 0;
            gz_out = gzopen(gz_name, "ab1");
            if (gz_out == NULL)
               gz_fail("gzopen");
            if (!registered) {
               atexit(finish_at_exit);
               registered = 1;
            }
         }
         while (done < gz_len) {
            size_t n = gz_len - done;
            int w;
            if (n > (1u << 30))
               n = 1u << 30;
            w = gzwrite(gz_out, gz_buf + done, (unsigned)n);
            if (w <= 0)
               gz_fail("gzwrite");
            done += (size_t)w;
         }
      }
      free(gz_buf);
      gz_buf = NULL;
      gz_len = 0;
   }
#endif
}

/* end the gzip member written so far (a no-op without PQ_GZ): called when
   a construction of descendants is complete, before GAP reads the file */

void FinishGAPLibrary(void)
{
#if defined(HAVE_ZLIB)
   if (gz_out != NULL) {
      int r = gzclose(gz_out);
      gz_out = NULL;
      if (r != Z_OK) {
         fprintf(stderr, "pq: gzclose %s failed (%d)\n", gz_name, r);
         exit(FAILURE);
      }
   }
#endif
}

#if defined(HAVE_ZLIB)
/* at exit: the same, but without calling exit again */
static void finish_at_exit(void)
{
   if (gz_out != NULL) {
      int r = gzclose(gz_out);
      gz_out = NULL;
      if (r != Z_OK)
         fprintf(stderr, "pq: gzclose %s failed (%d)\n", gz_name, r);
   }
}
#endif
