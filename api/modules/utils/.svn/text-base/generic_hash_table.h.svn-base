/*---------------------------------------------------------------------------*/
/*
 *  generic_hash_table.h
 *
 *  Created by Matthijs on Sat Feb 15 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */
/*---------------------------------------------------------------------------*/
/*! @defgroup utils Utils
 *
 *  @brief Define a set of functions wich are needed in the other part of the
 *  project.
 *
 *  @{
 */
/*---------------------------------------------------------------------------*/

#ifndef HASH_TABLE_INCLUDED
#define HASH_TABLE_INCLUDED

/*---------------------------------------------------------------------------
Generic hash table, with lots of void*'s and casts the memory for keys
and values is owned by the caller.  You can make a set of it by using
val==key. The hash_key should return a hash value for this key. The
cmp_key should return 0 for the same (like strcmp) 
---------------------------------------------------------------------------*/

typedef struct hash_table_t hash_table_t;

/*---------------------------------------------------------------------------*/

/*! @brief Initialize hash table */
hash_table_t *hash_table_init(long (*hash_key)(void *key),
                              int (*cmp_key)(void *key0,void *key1));

/*! @brief Insert element in hash table
 *
 * - if key already in ht, return associated val \n
 * - else, insert and return NULL
 */
void *hash_table_insert(hash_table_t *,void *key,void *val);

/*! @brief Insert element in hash table
 *
 * - if key already in ht, return associated val \n
 * - else, insert and return NULL
 */
void *hash_table_insert_weak(hash_table_t *,void *key,void *val);

/*! @brief Look up in the hash table
 *
 * return value or NULL if not found
 */
void *hash_table_lookup(hash_table_t *,void *key);

/*! @brief Get stats
 *
 * size=# of allocated cells \n
 * n_val=# of actual values 
 */
void hash_table_stats(hash_table_t *,int *size,int *n_val);

/*! @brief Dump hash table into file */
void hash_table_dump(hash_table_t *,char *fname);

/*---------------------------------------------------------------------------*/
/*
 * Use like
 * 
 * hash_table_it_t *it=hash_table_iterate(ht,NULL);
 * 
 * while(it) {
 *   ... use it->key and it->val ...
 *   it=hash_table_iterate(ht,it);
 * }
 * Iterator is dealloc'ed at end of iteration
 */
/*---------------------------------------------------------------------------*/

typedef struct {
  void *key;
  void *val;
} hash_table_it_t;

/*---------------------------------------------------------------------------*/

/*! @brief Iterator */
hash_table_it_t *hash_table_iterate(hash_table_t *,hash_table_it_t *);

/*! @brief Delete hash table
 *
 * If needed, keys & values may be deleted with hash_table_iterate.
 */
void hash_table_delete(hash_table_t *);

/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
