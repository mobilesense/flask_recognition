#include <stdlib.h>
#include <stdio.h>
#include "fast_alloc.h"

typedef struct block block;

typedef unsigned char bool;

#define TRUE 1
#define FALSE 0
#define NEW(type) ((type*)malloc(sizeof(type)))

struct block {
  block *next;
  bool dirty;
  char data[0];
};


struct fast_alloc_t {
  int us, upb;
  block *bl;
  char *fst;
};

fast_alloc_t *fa_alloc (int us, int upb)
{
  fast_alloc_t *fa;
  if (us < sizeof (char *))
    return NULL;
  fa = NEW (fast_alloc_t);
  if (!fa)
    return NULL;
  fa->us = us;
  fa->upb = upb;
  fa->bl = NULL;
  fa->fst = NULL;
  return fa;
}

static void fillin_block (fast_alloc_t * fa, block * bl)
{
  int u;
  char *p, *ul = fa->fst;
  bl->dirty = FALSE;
  for (u = 0, p = bl->data; u < fa->upb; u++, p += fa->us) {
    *(char **) p = ul;
    ul = p;
  }
  fa->fst = ul;
}

void *fa_allocunit (fast_alloc_t * fa)
{
  char *res;
  if (!fa->fst) {
    block **wp;
    for (wp = &fa->bl; *wp; wp = &(**wp).next)
      if ((**wp).dirty)
        break;
    if (!*wp) {
      *wp = malloc (sizeof (block) + fa->us * fa->upb);
      if (!*wp) {
        fprintf (stderr, "fa_allocunit: out of memory \n");
        exit (1);
      }
      (**wp).next = NULL;
    }
    fillin_block (fa, *wp);
  }
  res = fa->fst;
  fa->fst = *(char **) res;
  return res;
}

void fa_freeunit (fast_alloc_t * fa, void *unit)
{
  *(char **) unit = fa->fst;
  fa->fst = (char *) unit;
}

void fa_clear (fast_alloc_t * fa)
{
  block *bl;
  for (bl = fa->bl; bl; bl = bl->next) {
    if (bl->dirty)
      break;
    bl->dirty = TRUE;
  }
  fa->fst = NULL;
}

void fa_free (fast_alloc_t * fa)
{
  block *bl;
  for (bl = fa->bl; bl;) {
    block *prov = bl;
    bl = bl->next;
    free (prov);
  }
  free (fa);
}


void fa_describe (fast_alloc_t * fa)
{
  int nb, ndb, rm;
  block *bl;
  char *i;
  if (!fa) {
    printf ("fast_alloc_t is NULL\n");
    return;
  }
  printf ("fast_alloc_t, unit size=%d, unit/block=%d\n", fa->us, fa->upb);
  nb = 0;
  ndb = 0;
  for (bl = fa->bl; bl; bl = bl->next) {
    nb++;
    if (bl->dirty)
      ndb++;
  }
  printf (" %d used + %d dirty = %d total blocks\n", nb - ndb, ndb, nb);
  rm = 0;
  for (i = fa->fst; i; i = *(char **) i)
    rm++;
  printf (" %d units ready to serve\n", rm);
}

void fa_stats (fast_alloc_t * fa, int *n_block)
{
  if (n_block) {
    int nbl;
    block *bl;
    nbl = 0;
    for (bl = fa->bl; bl; bl = bl->next) {
      nbl++;
    }
    *n_block = nbl;
  }
}
