/*
 *  hash_table.c
 *  Sokoban
 *
 *  Created by Matthijs on Sat Feb 15 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#include "generic_hash_table.h"



#define NEW(type) (type*)malloc(sizeof(type))
#define NEWA(n,type) (type*)malloc(sizeof(type)*(n))



typedef struct {
  void *key;
  long hash;
  void *val;
  int link;
} entry_t;
/* 
states: 
- empty: key==NULL
- full: key!=NULL
linked: link!=-1
*/

#define INIT_LOG2_N 1

struct hash_table_t {
  int log2n;
  int p;
  entry_t *entries;

  long (*hash_key) (void *key);
  int (*cmp_key) (void *key0, void *key1);

};


hash_table_t *hash_table_init (long (*hash_key) (void *key),
                               int (*cmp_key) (void *key0, void *key1))
{
  hash_table_t *ret = NEW (hash_table_t);
  int i, n;
  ret->hash_key = hash_key;
  ret->cmp_key = cmp_key;

  ret->log2n = INIT_LOG2_N;
  n = 1 << INIT_LOG2_N;
  ret->p = 0;
  ret->entries = NEWA (n, entry_t);
  for (i = 0; i < n; i++)
    ret->entries[i].key = NULL;
  return ret;
}


static void copy_entry (entry_t * src, entry_t * dest)
{
  *dest = *src;
}

static void swap_entry (entry_t * e0, entry_t * e1)
{
  entry_t tmp = *e0;
  *e0 = *e1;
  *e1 = tmp;
}

#define BIGPRIME 47

static void double_size (hash_table_t * ht)
{
  int i;
  int n = 1 << ht->log2n;
  int old_mask = n - 1;
  int new_mask = 2 * n - 1;
  ht->entries = realloc (ht->entries, sizeof (entry_t) * n * 2);
  for (i = n; i < 2 * n; i++) {
    ht->entries[i].key = NULL;
  }
  /* first pass: handle normal keys, mark collisions with link=-2 */
  for (i = 0; i < n; i++) {
    entry_t *e = &ht->entries[i];
    if (e->key == NULL)
      continue;
    if ((e->hash & old_mask) == i) {    /* normal */
      if ((e->hash & new_mask) == i) {  /* don't move */
        e->link = -1;
      } else {                  /* move by n */
        entry_t *e1 = &ht->entries[i + n];
        copy_entry (e, e1);
        e1->link = -1;
        e->key = NULL;
      }
    } else {                    /* collided */
      e->link = -2;
    }
  }

  /* second pass: handle collided keys */
  for (i = 0; i < n; i++) {
    entry_t *e = &ht->entries[i], *e1;
    int j;
    if (e->key == NULL || e->link != -2)
      continue;
  redo:
    j = e->hash & new_mask;
    e1 = &ht->entries[j];
    if (e1->key == NULL) {
      copy_entry (e, e1);
      e1->link = -1;
      e->key = NULL;
    } else if (e1->link == -2) {
      swap_entry (e, e1);
      e1->link = -1;
      if ((e->hash & new_mask) == j) {
        /* very special case: j had to be swapped */
        e->link = -1;
      } else {
        e->link = -2;
        goto redo;
      }
    } else {                    /* follow list */
      while (e1->link != -1) {
        assert (e1->key != NULL && e1->link != -2);
        j = e1->link;
        e1 = &ht->entries[j];
      }
      /* insert myself into list */
      e1->link = i;
      e->link = -1;
    }
  }
  ht->log2n++;


}

int cmp_ptr (const void *ap, const void *bp)
{
  unsigned long a = *(unsigned long *) ap;
  unsigned long b = *(unsigned long *) bp;
  return a < b;
}

void *hash_table_insert (hash_table_t * ht, void *key, void *val)
{
  long hash = ht->hash_key (key);
  long mask = (1 << ht->log2n) - 1;
  void *old_val;
  int j;
  entry_t *e;

  j = hash & mask;
  e = &ht->entries[j];
  if (e->key == NULL) {
    e->link = -1;
    old_val = NULL;
  } else {
    for (;;) {
      assert (e->key != NULL);
      if (e->hash == hash && !ht->cmp_key (e->key, key)) {
        old_val = e->val;
        break;
      }
      if (e->link == -1) {
        entry_t *e1;
        do {
          j = (j + BIGPRIME) & mask;
          e1 = &ht->entries[j];
        } while (e1->key != NULL);
        e->link = j;
        e = e1;
        old_val = NULL;
        e->link = -1;
        goto e_ready;
      }
      /* else try next one */
      j = e->link;
      e = &ht->entries[j];
    }
  }
e_ready:
  e->key = key;
  e->hash = hash;
  e->val = val;

  if (old_val == NULL)
    ht->p++;

  if (ht->p * 3 > (1 << ht->log2n) * 2) {
    double_size (ht);
  }

  return old_val;
}

void *hash_table_insert_weak (hash_table_t * ht, void *key, void *val)
{
  long hash = ht->hash_key (key);
  long mask = (1 << ht->log2n) - 1;
  void *old_val;
  int j;
  entry_t *e;

  j = hash & mask;
  e = &ht->entries[j];
  if (e->key == NULL) {
    e->link = -1;
    old_val = NULL;
  } else {
    for (;;) {
      assert (e->key != NULL);
      if (e->hash == hash && !ht->cmp_key (e->key, key)) {
        old_val = e->val;
        break;
      }
      if (e->link == -1) {
        entry_t *e1;
        do {
          j = (j + BIGPRIME) & mask;
          e1 = &ht->entries[j];
        } while (e1->key != NULL);
        e->link = j;
        e = e1;
        old_val = NULL;
        e->link = -1;
        goto e_ready;
      }
      /* else try next one */
      j = e->link;
      e = &ht->entries[j];
    }
  }
e_ready:

  if (old_val == NULL) {
    e->key = key;
    e->hash = hash;
    e->val = val;
    ht->p++;
  }
  if (ht->p * 3 > (1 << ht->log2n) * 2) {
    /*    hash_table_check(ht); */
    double_size (ht);
    /* hash_table_check(ht); */
  }

  return old_val;
}


void hash_table_check (hash_table_t * ht)
{
  /* check that all keys are distinct */
  void *kt[ht->p];
  int i, p = 0;

  for (i = 0; i < (1 << ht->log2n); i++) {
    entry_t *e = ht->entries + i;
    if (e->key) {
      assert (hash_table_lookup (ht, e->key) == e->val);
      kt[p++] = e;
    }
  }
  assert (p == ht->p);
  qsort (kt, p, sizeof (void *), cmp_ptr);
  for (i = 1; i < p; i++)
    assert (kt[i] != kt[i - 1]);

}


void *hash_table_lookup (hash_table_t * ht, void *key)
{
  long hash = ht->hash_key (key);
  long mask = (1 << ht->log2n) - 1;
  int j;
  entry_t *e;
  j = hash & mask;
  do {
    e = &ht->entries[j];
    if (e->key == NULL)
      return NULL;
    if (e->hash == hash && !ht->cmp_key (e->key, key))
      return e->val;
    j = e->link;
  } while (j != -1);
  return NULL;
}


void dump_hash_table (hash_table_t * ht)
{
  int i, n = 1 << ht->log2n;
  int ncoll, ll_max, ll_sum;
  int mask = (1 << ht->log2n) - 1;
  printf ("hash_table {\n");
  printf ("  log2n=%d (n=%d)\n" "  p=%d\n", ht->log2n, 1 << ht->log2n, ht->p);
  ncoll = ll_max = ll_sum = 0;
  for (i = 0; i < n; i++) {
    entry_t *e = &ht->entries[i];
    printf ("  entry[%d]={\n", i);
    if (e->key == NULL)
      printf ("    key=NULL\n");
    else {
      printf ("    key=%p\n", e->key);
      printf ("    hash=%08lx\n", e->hash);
      printf ("    val=%p\n", e->val);
      printf ("    link=%d\n", e->link);
      if ((e->hash & mask) != i) {
        int j, ll;
        ncoll++;
        ll = 0;
        for (;;) {
          j = e->link;
          if (j == -1)
            break;
          e = &ht->entries[j];
          ll++;
        }
        if (ll > ll_max)
          ll_max = ll;
        ll_sum += ll;
      }
    }
    printf ("  }\n");
  }
  printf ("  [collisions: %d]\n", ncoll);
  printf ("  [link len: max=%d avg=%g]\n", ll_max, ll_sum / (float) ncoll);
  printf ("}\n");

}


void hash_table_stats (hash_table_t * ht, int *size, int *n_val)
{
  if (size)
    *size = 1 << ht->log2n;
  if (n_val)
    *n_val = ht->p;
}


void hash_table_dump (hash_table_t * ht, char *fname)
{
  int i, n = 1 << ht->log2n;
  FILE *f = fopen (fname, "w");
  assert (f);

  for (i = 0; i < n; i++) {
    if (ht->entries[i].key)
      fprintf (f, "%p\n", ht->entries[i].key);
  }
  fclose (f);

}


typedef struct {
  hash_table_it_t cur;          /* current valid entry */
  int i;                        /* index of cur */
} it_t;

hash_table_it_t *hash_table_iterate (hash_table_t * ht, hash_table_it_t * hit)
{
  it_t *it = (it_t *) hit;
  int n = 1 << ht->log2n;
  if (!it) {
    it = NEW (it_t);
    it->i = -1;
  }
  for (;;) {
    it->i++;
    if (it->i == n) {
      free (it);
      return NULL;
    }
    entry_t *e = ht->entries + it->i;
    if (e->key != NULL) {
      it->cur.key = e->key;
      it->cur.val = e->val;
      return &it->cur;
    }
  }
}


void hash_table_delete (hash_table_t * ht)
{
  free (ht->entries);
  free (ht);
}



#ifdef TEST


long hash_key (char *key)
{
  int l = strlen (key);
  int i;
  unsigned long hash = 0x1528364;
  for (i = 0; i < l; i++) {
    hash = ((hash << 13) | (hash >> 19)) ^ key[i];
  }
  return hash;
}

int cmp_key (char *key0, char *key1)
{
  return strcmp (key0, key1);
}



int main ()
{
  int i;
  hash_table_t *ht = hash_table_init ();
  for (i = 0; i < 2000; i++) {
    char *key = malloc (10);
    void *val = (void *) (i + 50);
    sprintf (key, "%d %x", i, i);
    printf ("%d -> %08x\n", i, hash_key (key));
    hash_table_insert (ht, key, val);
  }
  dump_hash_table (ht);
  for (i = 0; i < 2000; i++) {
    char *key = malloc (10);
    void *val;
    sprintf (key, "%d %x", i, i);
    val = hash_table_lookup (ht, key);
    assert ((long) val == i + 50);
    printf ("%d -> \"%s\" -> %p\n", i, key, val);

    free (key);
  }
  return 0;
}



#endif
