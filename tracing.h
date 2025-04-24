// vim:ft=c
#ifndef TRACING_H
#define TRACING_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#ifndef _USECFMT
#ifdef __APPLE__
#define _USECFMT "d"
#else
#define _USECFMT "ld"
#endif
#endif /* ifndef _USECFMT */

/*
 * Add `y` to `x`
 */
static inline void timeval_add(struct timeval *x, struct timeval *y) {
  x->tv_sec += y->tv_sec;
  x->tv_usec += y->tv_usec;
  if (x->tv_usec >= 1000000) {
    x->tv_sec++;
    x->tv_usec -= 1000000;
  } else if (x->tv_usec <= -1000000) {
    x->tv_sec--;
    x->tv_usec += 1000000;
  }
}

static inline struct timeval timeval_subtract(struct timeval *x,
                                              struct timeval *y) {
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
  return (struct timeval){.tv_sec = x->tv_sec - y->tv_sec,
                          .tv_usec = x->tv_usec - y->tv_usec};
}

typedef struct {
  clock_t c;
  struct timeval t;
} trace_t;

/*
 * Create a new trace containing current time
 */
static inline trace_t trace_start() {
  clock_t c = clock();
  struct timeval t;
  gettimeofday(&t, NULL);
  return (trace_t){
      .c = c,
      .t = t,
  };
}
/*
 * Increment `dest` trace with the values computed from `now - start`
 */
static inline void trace_incnow(trace_t *dest, trace_t *start) {
  trace_t now = trace_start();
  struct timeval tdiff = timeval_subtract(&now.t, &start->t);
  timeval_add(&dest->t, &tdiff);
  dest->c += now.c - start->c;
}
/*
 * Print the tracing of `t - now` with `name` as recognizer
 */
static inline void printtrace_diffnow(trace_t *t, const char *name) {
  struct timeval t1, d;
  clock_t c1 = clock();
  gettimeofday(&t1, NULL);
  d = timeval_subtract(&t1, &t->t);
  fprintf(stderr, "[%s] time: %8ld:%04" _USECFMT " s:us\t CPU: %8ld s\n", name,
          d.tv_sec, d.tv_usec, (c1 - t->c) / CLOCKS_PER_SEC);
}
/*
 * Print the tracing `t`, meant to be used with incremental tracing
 * `trace_inc*()` functions
 */
static inline void printtrace(trace_t *t, const char *name) {
  fprintf(stderr, "[%s] time: %8ld:%04" _USECFMT " s:us\t CPU: %8ld s\n", name,
          t->t.tv_sec, t->t.tv_usec, t->c / CLOCKS_PER_SEC);
}

/*static inline void _pTRACE(struct timeval d, char *name) {*/
/*  fprintf(stderr, "[%s] time: %8ld:%04ld s:us\t CPU: %8ld s\n", "" #f "",*/
/*          d.tv_sec, d.tv_usec, (c1 - c0) / CLOCKS_PER_SEC);*/
/*}*/

/*
 * Trace computation time of function `f`, if the function has
 * a return value, it should be passed in the variadic argument.
 *
 * Examples:
 * void fun_noret(int x, int y);
 * int fun_ret(double z);
 *
 * TRACE(fun_noret(1, 2));
 * int res;
 * TRACE(fun_ret(0.1), res);
 * // res contains the return value of `fun_ret`
 *
 *
 */
#define TRACE(f, ...)                                                          \
  do {                                                                         \
    struct timeval t0, t1, d;                                                  \
    clock_t c0 = clock(), c1;                                                  \
    gettimeofday(&t0, NULL);                                                   \
    __VA_OPT__(__VA_ARGS__ =)(f);                                              \
    gettimeofday(&t1, NULL);                                                   \
    d = timeval_subtract(&t1, &t0);                                            \
    c1 = clock();                                                              \
    fprintf(stderr, "[%s] time: %8ld:%04" _USECFMT " s:us\t CPU: %8ld s\n",    \
            "" #f "", d.tv_sec, d.tv_usec, (c1 - c0) / CLOCKS_PER_SEC);        \
  } while (0)

#endif // !TRACING_H
