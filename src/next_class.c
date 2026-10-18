/****************************************************************************
**
*A  next_class.c                ANUPQ source                   Eamonn O'Brien
**
*Y  Copyright 1995-2001,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
*Y  Copyright 1995-2001,  School of Mathematical Sciences, ANU,     Australia
**
*/

#include "pq_defs.h"
#include "constants.h"
#include "pcp_vars.h"
#include "pq_functions.h"
#include "exp_vars.h"
#define BOTH_TAILS 0

#if defined(GROUP)

/* calculate the next class of the group layer by layer */

void next_class(Logical report, int **head, int **list, struct pcp_vars *pcp)
{
   register int *y = y_address;

   register int class;
   struct exp_vars exp_flag;
   int prev = 1, new;

   /* if class 1 computation, setup has already been done --
      before relations are read */
   if (y[pcp->clend + 1] != 0)
      setup(pcp);

   if (pcp->overflow || (pcp->complete != 0 && !pcp->multiplicator))
      return;

   if (pcp->extra_relations != 0 || pcp->m != 0)
      initialise_exponent(&exp_flag, pcp);

   for (class = pcp->cc; class > 1; --class) {
      tails(BOTH_TAILS, class, pcp->cc, 1, pcp);
      if (pcp->overflow)
         return;
      /* the queue holds redundant generators of the current class, found
         by consistency and later by enforce_exponent; now that tails has
         defined them, make room for all of them */
      if (pcp->m != 0) {
         new = pcp->lastg - pcp->ccbeg + 1;
         exp_flag.queue =
             reallocate_vector(exp_flag.queue, prev, new, 1, FALSE);
         prev = new;
      }
      if (class != 2) {
         consistency(0, exp_flag.queue, &exp_flag.queue_length, class, pcp);
         if (pcp->overflow || (pcp->complete != 0 && !pcp->multiplicator))
            return;
      }
   }

   if (!pcp->multiplicator) {
      if (pcp->cc > 1) {
         update_generators(pcp);
         if (pcp->overflow)
            return;

         /* this call is not necessary if there are no automorphisms
            present -- however, it may be useful in reducing space
            requirements and hence improve efficiency -- this needs
            to be investigated further -- EO'B October 1991 */

         if (pcp->m != 0) {
            eliminate(0, pcp);
            /* an elimination has been performed -- must reset queue length */
            exp_flag.queue_length = 0;
            if (pcp->overflow)
               return;
         }
      }

      collect_relations(pcp);
      if (pcp->overflow || pcp->complete != 0 || !pcp->valid)
         return;
      if (pcp->extra_relations != 0) {
         if (pcp->m == 0) {
            exp_flag.list = ALL_WORDS;
            extra_relations(&exp_flag, pcp);
         } else {
            exp_flag.list = REDUCED_LIST;
            enforce_exponent(report, &exp_flag, head, list, pcp);
            if (pcp->overflow)
               return;
         }
      }

      if (pcp->overflow || pcp->complete != 0 || !pcp->valid)
         return;
   }


   /* if the multiplicator flag is set and there are redundant
      generators, then we must be careful about elimination --
      update_generators has not yet been performed; see note
      in code of that procedure */

   if (pcp->multiplicator && pcp->ndgen > y[pcp->clend + 1])
      eliminate(TRUE, pcp);
   else
      eliminate(FALSE, pcp);
}

#endif


/* largest p^(d^2) and |GL(d, p)| for which induce_full_gl decides */
#define GL_MAX_CODES (1L << 27)
#define GL_MAX_ORDER (1L << 24)

/* return TRUE if the automorphisms, whose action is described by head
   and list, induce the full general linear group GL(d, p) on the Frattini
   quotient of rank d; return FALSE if they do not, or if this is too
   expensive to decide; each d x d matrix over GF(p) is coded as an
   integer in base p, and the products of the automorphism matrices are
   enumerated until there are |GL(d, p)| of them */

static Logical induce_full_gl(int *head, int *list, struct pcp_vars *pcp)
{
   register int *y = y_address;

   int d = y[pcp->clend + 1];
   int p = pcp->p;
   int n = d * d;
   int m = pcp->m;
   long nmr_codes = 1, gl_order = 1, pd = 1, pi = 1, count, next, code;
   int *mats, *a, *b, *orbit = NULL;
   unsigned char *seen = NULL;
   int alpha, i, j, k, r, ptr, length, gen;
   long inv, t;
   Logical result = FALSE;

#include "access.h"

   if (d == 0 || m == 0)
      return FALSE;

   for (i = 0; i < n; ++i)
      if ((nmr_codes *= p) > GL_MAX_CODES)
         return FALSE;
   for (i = 0; i < d; ++i)
      pd *= p;
   for (i = 0; i < d; ++i, pi *= p)
      if ((gl_order *= pd - pi) > GL_MAX_ORDER)
         return FALSE;

   /* the matrices of the automorphisms, then two work matrices */
   mats = (int *)calloc((size_t)(m + 2) * n, sizeof(int));
   if (mats == NULL)
      return FALSE;
   a = mats + m * n;
   b = a + n;

   for (alpha = 0; alpha < m; ++alpha) {
      for (i = 0; i < d; ++i) {
         ptr = head[alpha * head[0] + i + 1];
         length = list[ptr + 1];
         for (k = ptr + 2; k <= ptr + length + 1; ++k) {
            gen = FIELD2(list[k]);
            if (gen <= d)
               mats[alpha * n + i * d + gen - 1] = FIELD1(list[k]) % p;
         }
      }

      /* each matrix must be invertible */
      for (i = 0; i < n; ++i)
         a[i] = mats[alpha * n + i];
      for (j = 0; j < d; ++j) {
         for (r = j; r < d && a[r * d + j] == 0; ++r)
            ;
         if (r == d)
            goto done;
         for (k = 0; k < d; ++k) {
            t = a[j * d + k];
            a[j * d + k] = a[r * d + k];
            a[r * d + k] = t;
         }
         for (inv = 1; inv * a[j * d + j] % p != 1; ++inv)
            ;
         for (r = j + 1; r < d; ++r) {
            t = a[r * d + j] * inv % p;
            for (k = j; k < d; ++k)
               a[r * d + k] = (a[r * d + k] + (p - t) * a[j * d + k]) % p;
         }
      }
   }

   seen = (unsigned char *)calloc((size_t)(nmr_codes / 8 + 1), 1);
   orbit = (int *)malloc((size_t)gl_order * sizeof(int));
   if (seen == NULL || orbit == NULL)
      goto done;

   /* start with the identity matrix */
   for (code = 0, i = d - 1; i >= 0; --i)
      for (j = d - 1; j >= 0; --j)
         code = code * p + (i == j);
   orbit[0] = code;
   seen[code >> 3] |= 1 << (code & 7);
   count = 1;

   for (next = 0; next < count && count < gl_order; ++next) {
      for (code = orbit[next], i = 0; i < n; ++i, code /= p)
         a[i] = code % p;
      for (alpha = 0; alpha < m && count < gl_order; ++alpha) {
         for (i = 0; i < d; ++i)
            for (j = 0; j < d; ++j) {
               for (t = 0, k = 0; k < d; ++k)
                  t += a[i * d + k] * mats[alpha * n + k * d + j];
               b[i * d + j] = t % p;
            }
         for (code = 0, i = n - 1; i >= 0; --i)
            code = code * p + b[i];
         if (!(seen[code >> 3] & (1 << (code & 7)))) {
            seen[code >> 3] |= 1 << (code & 7);
            orbit[count++] = code;
         }
      }
   }

   /* the matrices are invertible, so they generate a subgroup of
      GL(d, p) with at least |GL(d, p)| elements */
   result = (count == gl_order);

done:
   free(orbit);
   free(seen);
   free(mats);
   return result;
}


/* when automorphisms are supplied, enforce exponent law by setting
   up a queue of redundant generators and then close this queue */

void enforce_exponent(Logical report,
                      struct exp_vars *exp_flag,
                      int **head,
                      int **list,
                      struct pcp_vars *pcp)
{
   register int *y = y_address;

#if defined(TIME)
   int t;
#endif
   int factor;
   int limit = 0;

   int *queue, queue_length;
   int list_length;
   char *s;

   if (pcp->m != 0) {
      read_value(TRUE, "Input queue factor: ", &factor, 0);
      limit = factor * (pcp->lastg - pcp->ccbeg + 1) / 100;
   }

#if defined(TIME)
   t = runTime();
#endif

   Extend_Auts(head, list, y[pcp->clend + 1] + 1, pcp);
   if (pcp->overflow)
      return;

#if defined(TIME)
   t = runTime() - t;
   printf("Time to extend automorphisms is %.2f seconds\n", t * CLK_SCALE);
#endif

#if defined(TIME)
   t = runTime();
#endif

   list_length = pcp->lastg - pcp->ccbeg + 1;

   /* the reduced list of test words, and the bound 2c/3 on their weight,
      are valid only if the group is relatively free and the automorphisms
      induce the full general linear group on its Frattini quotient;
      otherwise test all words */
   if (exp_flag->list == REDUCED_LIST &&
       (pcp->ndrel != 0 || pcp->metabelian ||
        !induce_full_gl(*head, *list, pcp)))
      exp_flag->list = ALL_WORDS;

   if (exp_flag->list == REDUCED_LIST) {
      pcp->start_wt = 1;
      pcp->end_wt = (2 * pcp->cc) / 3;
   }

   extra_relations(exp_flag, pcp);
   if (pcp->overflow)
      return;
   queue = exp_flag->queue;
   queue_length = exp_flag->queue_length;

   pcp->start_wt = 0;
   pcp->end_wt = 0;
#if defined(TIME)
   t = runTime() - t;
#endif

   s = (queue_length == 1) ? "y" : "ies";
   if (report || pcp->fullop || pcp->diagn)
      printf("Exponent checks gave %d redundanc%s\n", queue_length, s);

#if defined(TIME)
   printf("Time to check exponents is %.2f seconds\n", t * CLK_SCALE);
#endif

   close_queue(
       report, list_length, limit, *head, *list, queue, queue_length, pcp);
}

/* close the queue of redundant generators under the action
   of the automorphisms */

void close_queue(Logical report,
                 int list_length,
                 int limit,
                 int *head,
                 int *list,
                 int *queue,
                 int queue_length,
                 struct pcp_vars *pcp)
{
   int *long_queue, long_queue_length;
#if defined(TIME)
   int t;
#endif

   long_queue = allocate_vector(list_length, 1, 0);
   long_queue_length = 0;

#if defined(TIME)
   t = runTime();
#endif

   if (!pcp->complete) {
      close_relations(report,
                      limit,
                      1,
                      head,
                      list,
                      queue,
                      queue_length,
                      long_queue,
                      &long_queue_length,
                      pcp);
      if (report || pcp->fullop || pcp->diagn)
         printf("Length of long queue after closing short queue is %d\n",
                long_queue_length);
   }

   if (!pcp->complete)
      close_relations(report,
                      limit,
                      2,
                      head,
                      list,
                      long_queue,
                      long_queue_length,
                      long_queue,
                      &long_queue_length,
                      pcp);

   if (report || pcp->fullop || pcp->diagn)
      printf("Final long queue length is %d\n", long_queue_length);

#if defined(TIME)
   t = runTime() - t;
   printf("Time to close under action of automorphisms is %.2f seconds\n",
          t * CLK_SCALE);
#endif

   free_vector(queue, 1);
   free_vector(long_queue, 1);
}

int rearrange_queues(int limit,
                     int *queue,
                     int *queue_length,
                     int *long_queue,
                     int long_queue_length,
                     struct pcp_vars *pcp)
{
   register int *y = y_address;

   int gen;
   int i, p1;

   for (i = 1; i <= long_queue_length; ++i) {
      gen = long_queue[i];
      p1 = -y[pcp->structure + gen];
      if (y[p1 + 1] < limit) {
         queue[++*queue_length] = gen;
         long_queue[i] = 0;
      }
   }
   return 0;
}

/* sort queue of redundant generators according to increasing length
   of those relations which imply that each generator is redundant */

void bubble_sort(int *x, int len, struct pcp_vars *pcp)
{
   register int *y = y_address;

   register int i, j, pointer, temp;
   Logical swap = TRUE;
   register int structure = pcp->structure;
   int *length = allocate_vector(len, 1, 1);

   /* set up the length of the relations as an array */
   for (i = 1; i <= len; ++i) {
      if (x[i] != 0) {
         pointer = -y[structure + x[i]];
         length[i] = y[pointer + 1];
      }
   }

   /* now sort the queue of generators */
   for (i = 1; i <= len && swap; ++i) {
      swap = FALSE;
      for (j = len; j > i; --j)
         if (length[j] < length[j - 1]) {
            temp = x[j];
            x[j] = x[j - 1];
            x[j - 1] = temp;
            temp = length[j];
            length[j] = length[j - 1];
            length[j - 1] = temp;
            swap = TRUE;
         }
   }

   free_vector(length, 1);
}
