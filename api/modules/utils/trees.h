#ifndef TREES_H_INCLUDED
#define TREES_H_INCLUDED



/*---------------------------------------------------------------------------
  AVL tree funtions. The trees store an long->float mapping (key to
  value). Find, insert, remove, max are all in O(log(n)) 
 ---------------------------------------------------------------------------*/

typedef struct avl_tree_t avl_tree_t;

avl_tree_t* avl_tree_new(void);

void avl_tree_delete(avl_tree_t*);

void avl_tree_add(avl_tree_t*,long k,float v);

void avl_tree_remove(avl_tree_t*,long k,float v);

void avl_tree_print(const avl_tree_t*);

void avl_tree_check(const avl_tree_t*);

void avl_tree_max(const avl_tree_t*,long *k_out,float *v_out);


#endif
