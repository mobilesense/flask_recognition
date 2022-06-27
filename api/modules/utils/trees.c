
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "trees.h"
#include "fast_alloc.h"


#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(sizeof(type),(n))
#define NEW(type) NEWA(type,1)
#define MAX(a,b) ((a)>(b) ? (a) : (b))



/*******************************************************************
 * types 
 */

typedef struct {
  long k;
  float v;
} avl_payload_t;

typedef struct avl_node_t avl_node_t;

struct avl_node_t {
  avl_payload_t v;
  int depth;
  avl_node_t *left,*right;
};


struct avl_tree_t {
  fast_alloc_t *fa_nodes;
  avl_node_t *root;
};




/*******************************************************************
 * elementaty operations 
 */

#define DEPTH(n) ((n) ? (n)->depth : 0)

static void avl_swap_payloads(avl_node_t *a,avl_node_t*b) {
  avl_payload_t tmp=a->v;
  a->v=b->v;
  b->v=tmp;
}

static void avl_rebalance_big_left(avl_node_t *node) {
  /* make a subtree balanced if the left sub-tree has depth 2+ depth of right subtree */
  avl_node_t *node2=node->left;
  int ldepth=node2->depth;
  int lldepth=DEPTH(node2->left);
  if(lldepth==ldepth-1) { /* 2-node rotation */
    avl_swap_payloads(node,node2);
    node->left=node2->left;
    node2->left=node2->right;
    node2->right=node->right;
    node->right=node2;
    node2->depth=DEPTH(node2->left)+1;
    node->depth=node2->depth+1;
  } else { /* 3-node rotation */
    avl_node_t *node3=node2->right;
    avl_swap_payloads(node,node3);
    node2->right=node3->left;
    node3->left=node3->right;
    node3->right=node->right;
    node->right=node3;
    node2->depth=ldepth-1;
    node3->depth=ldepth-1;
    node->depth=ldepth;
  }
}


static void avl_rebalance_big_right(avl_node_t *node) {
  /* make a subtree balanced if the right sub-tree has depth 2+ depth of left subtree */
  avl_node_t *node2=node->right;
  int ldepth=node2->depth;
  int lldepth=DEPTH(node2->right);
  if(lldepth==ldepth-1) { /* 2-node rotation */
    avl_swap_payloads(node,node2);
    node->right=node2->right;
    node2->right=node2->left;
    node2->left=node->left;
    node->left=node2;
    node2->depth=DEPTH(node2->right)+1;
    node->depth=node2->depth+1;
  } else { /* 3-node rotation */
    avl_node_t *node3=node2->left;
    avl_swap_payloads(node,node3);
    node2->left=node3->right;
    node3->right=node3->left;
    node3->left=node->left;
    node->left=node3;
    node2->depth=ldepth-1;
    node3->depth=ldepth-1;
    node->depth=ldepth;
  }
}

/* >0 iff a>b */
static int avl_cmp_payloads(const avl_payload_t *a,const avl_payload_t *b) {
/*
  float epsilon=MAX(fabs(a->v),b->v)*1e-5;
*/
  float epsilon=0;

  if(a->v > b->v+epsilon) return 1;
  if(a->v+epsilon < b->v) return -1;

  if(a->k > b->k) return 1;
  if(a->k < b->k) return -1;
  return 0;
}

/*******************************************************************
 * constructor/destructor, scan 
 */

avl_tree_t* avl_tree_new(void) {
  avl_tree_t* t=NEW(avl_tree_t);
  t->root=NULL;
  t->fa_nodes=fa_alloc(sizeof(avl_node_t),1024);
  return t;
}

void avl_tree_delete(avl_tree_t*t) {
  fa_free(t->fa_nodes);
  free(t);
}

static void avl_check_rec(const avl_node_t *node) {
  if(!node) return;
  assert(!node->left || avl_cmp_payloads(&node->left->v,&node->v)<0);
  assert(!node->right || avl_cmp_payloads(&node->right->v,&node->v)>0);
  assert(abs(DEPTH(node->left)-DEPTH(node->right))<=1);
  assert(node->depth==1+MAX(DEPTH(node->left),DEPTH(node->right)));
  avl_check_rec(node->left);
  avl_check_rec(node->right);
}

void avl_tree_check(const avl_tree_t*t) {
  avl_check_rec(t->root);
}


static const avl_node_t *avl_max_rec(const avl_node_t*node) {
  while(node->right) node=node->right;
  return node;
}

void avl_tree_max(const avl_tree_t*t,long *k_out,float *v_out) {
  assert(t->root || !"max of empty tree");
  const avl_node_t *node=avl_max_rec(t->root);
  if(k_out) *k_out=node->v.k;
  if(v_out) *v_out=node->v.v;
}

static void avl_print_rec(int n_indent,const avl_node_t *node) {
  if(!node) return;
  avl_print_rec(n_indent+1,node->left);
  int i;
  for(i=0;i<n_indent;i++) printf("  ");
  printf("k=%ld v=%g depth=%d\n",node->v.k,node->v.v,node->depth);
  avl_print_rec(n_indent+1,node->right);
}


void avl_tree_print(const avl_tree_t*t) {
  avl_print_rec(0,t->root);
}



/*******************************************************************
 * add 
 */

static avl_node_t *avl_add_rec(avl_node_t *v,avl_node_t *node) {
  if(!node) return v;
  int cmp=avl_cmp_payloads(&v->v,&node->v);
  assert(cmp!=0 || !"value already in tree");
  if(cmp<0) {
    node->left=avl_add_rec(v,node->left);
    int rightdepth=DEPTH(node->right);
    if(node->left->depth==rightdepth+2)
      avl_rebalance_big_left(node);
    else {
      assert(abs(node->left->depth-rightdepth)<=1);
      node->depth=1+MAX(node->left->depth,rightdepth);
    }
  } else {
    node->right=avl_add_rec(v,node->right);
    int leftdepth=DEPTH(node->left);
    if(node->right->depth==leftdepth+2)
      avl_rebalance_big_right(node);
    else {
      assert(abs(node->right->depth-leftdepth)<=1);
      node->depth=1+MAX(node->right->depth,leftdepth);
    }
  }
  return node;
}



void avl_tree_add(avl_tree_t*t,long k,float v) {
  avl_node_t *newnode=fa_allocunit(t->fa_nodes);
  newnode->depth=1;
  newnode->left=NULL;
  newnode->right=NULL;
  newnode->v.k=k;
  newnode->v.v=v;
  t->root=avl_add_rec(newnode,t->root);
}


/*******************************************************************
 * remove */

static avl_node_t *avl_remove_rec(avl_tree_t *t,avl_payload_t *v,avl_node_t *node) {
  assert(node || !"node to remove not found");
  int cmp=avl_cmp_payloads(v,&node->v);

  if(cmp==0) {
    avl_node_t *replacement=node;
    if(!node->left) {
      if(!node->right)
        replacement=NULL;
      else
        replacement=node->right;
    } else {
      if(!node->right)
        replacement=node->left;
      else {
        avl_node_t *node2=(avl_node_t *)avl_max_rec(node->left);
        node->v=node2->v;
        /* now node2 has to be removed. Re-inject in algo... */
        v=&node->v;
        cmp=-1;
        goto newtarget;
      }
    }
    fa_freeunit(t->fa_nodes,node);
    return replacement;
  }

newtarget:
  if(cmp<0) {
    node->left=avl_remove_rec(t,v,node->left);
    int leftdepth=DEPTH(node->left);
    int rightdepth=DEPTH(node->right);
    if(leftdepth==rightdepth-2)
      avl_rebalance_big_right(node);
    else {
      node->depth=1+MAX(leftdepth,rightdepth);
      assert(abs(leftdepth-rightdepth)<=1);
    }
  } else {
    node->right=avl_remove_rec(t,v,node->right);
    int leftdepth=DEPTH(node->left);
    int rightdepth=DEPTH(node->right);
    if(leftdepth-2==rightdepth)
      avl_rebalance_big_left(node);
    else {
      node->depth=1+MAX(leftdepth,rightdepth);
      assert(abs(leftdepth-rightdepth)<=1);
    }
  }

  return node;
}


void avl_tree_remove(avl_tree_t*t,long k,float v) {
  avl_payload_t toremove={k,v};
  t->root=avl_remove_rec(t,&toremove,t->root);
}


