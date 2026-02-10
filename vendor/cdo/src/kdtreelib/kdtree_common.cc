/*!\file kdtree_common.c
 * \brief Routines common to the Cartesian and spherical kd-tree versions
 *
 */

#include "kdtree.h"
#include "pqueue.h"


extern int pmergesort(struct kd_point *base, size_t nmemb, int axis, int max_threads);


/* ********************************************************************

   general utility functions 

   ********************************************************************* */

void *
kd_malloc(size_t size, const char *msg)
{
  void *ptr;
  if ((ptr = malloc(size)) == NULL)
    perror(msg);
  return ptr;
}

kdata_t
kd_sqr(kdata_t x)
{
  return (x < 0 || x > 0) ? x * x : 0;
}

int
kd_isleaf(struct kdNode *n)
{
  return (n->left || n->right) ? 0 : 1;
}

/* end utility functions */

/* *******************************************************************
   
   Helper functions for debugging

   ******************************************************************* */

void
kd_printNode(struct kdNode *node)
{
  if (!node) {
    fprintf(stderr, "Node is empty.\n");
    return;
  }
  printf("Node %p at (%f, %f)\n", (void *) node, node->location[0],
         node->location[1]);

  printf("Split axis: %d\n", node->split);
  printf("Corners: (%f, %f)\t(%f, %f)\n", node->min[0], node->min[1],
         node->max[0], node->max[1]);
  printf("Children: %p\t%p\n", (void *) node->left, (void *) node->right);
  printf("Index: %zu\n", node->index);

  printf("\n");
}

void
kd_printTree(struct kdNode *node)
{
  if ( node == NULL ) return;

  kd_printTree(node->left);

  if ( kd_isleaf(node) )
    printf("%f\t%f\n", node->location[0], node->location[1]);

  kd_printTree(node->right);
}

/* End helper functions */

void kd_initArg(struct kd_thread_data *d, kdNodePool *nodepool, struct kd_point *points, size_t nPoints,
                kdata_t *min, kdata_t *max, int depth, int max_threads, int dim)
{
  d->nodepool = nodepool;
  d->points = points;
  d->nPoints = nPoints;
  memcpy(d->min, min, dim*sizeof(kdata_t));
  memcpy(d->max, max, dim*sizeof(kdata_t));
  d->depth = depth;
  d->max_threads = max_threads;
  d->dim = dim;
}

/* ******************************************************************
   
   Functions for building and destroying trees 

   ******************************************************************** */

/*!
 * \brief free the kd-tree data structure, 
 *
 * \param node the root node of the tree to be destroyed
 *
 * \param *destr a pointer to the destructor function for the data container.
 *
 * \return This function does not return a value
 */

static
void kd_freeNode(kdNode *node)
{
  if ( node ) free(node);
}


kdNodePool_t *kd_nodepool_new(size_t n)
{
  kdNodePool_t *nodepool = (kdNodePool_t *) malloc(sizeof(kdNodePool_t));
  nodepool->size = n;
  nodepool->imp = 0;
  nodepool->pool = (kdNode *) malloc(n*sizeof(kdNode));
  printf("kd_nodepool_new: size=%zu\n", n);
  return nodepool;
}


void kd_nodepool_delete(kdNodePool_t *nodepool)
{
  if ( nodepool == NULL ) return;

  if ( nodepool->pool ) free(nodepool->pool);
  free(nodepool);
}

struct kdNode *
kd_allocNode(struct kd_point *points, size_t pivot, kdata_t *min, kdata_t *max, int axis, int dim)
{
  struct kdNode *node;

  if ((node = (kdNode *)kd_malloc(sizeof(kdNode), "kd_allocNode (node): ")) == NULL)
    return NULL;

  node->split = axis;
  memcpy(node->location, points[pivot].point, dim * sizeof(kdata_t));
  memcpy(node->min, min, dim * sizeof(kdata_t));
  memcpy(node->max, max, dim * sizeof(kdata_t));
  node->left = node->right = NULL;
  node->index = 0;

  return node;
}

struct kdNode *
kd_allocNodeP(kdNodePool *nodepool, struct kd_point *points, size_t pivot, kdata_t *min, kdata_t *max, int axis, int dim)
{
  struct kdNode *node;

  if ( nodepool )
    {
      if ( nodepool->imp >= nodepool->size ) return NULL;
      node = &nodepool->pool[nodepool->imp++];
    }
  else
    {
      if ((node = (kdNode *)kd_malloc(sizeof(kdNode), "kd_allocNode (node): ")) == NULL)
        return NULL;
    }

  node->split = axis;
  memcpy(node->location, points[pivot].point, dim * sizeof(kdata_t));
  memcpy(node->min, min, dim * sizeof(kdata_t));
  memcpy(node->max, max, dim * sizeof(kdata_t));
  node->left = node->right = NULL;
  node->index = 0;

  return node;
}

void
kd_doDestroyTree(struct kdNode *node)
{
  if ( node == NULL ) return;

  kd_doDestroyTree(node->left);
  kd_doDestroyTree(node->right);
  kd_freeNode(node);
}

void
kd_destroyTree(kdTree_t *tree)
{
  if ( tree == NULL ) return;

  if ( tree->nodepool )
    kd_nodepool_delete(tree->nodepool);
  else
    kd_doDestroyTree(tree->node);

  free(tree);
}

void *kd_doBuildTree(void *threadarg)
{
  kdata_t tmpMaxLeft[KD_MAX_DIM], tmpMinRight[KD_MAX_DIM];
  struct kdNode *node;
  pthread_t threads[2];
  pthread_attr_t attr;
  struct kd_thread_data argleft, argright;
  struct kd_thread_data *my_data = (struct kd_thread_data *) threadarg;

  kdNodePool_t *nodepool = my_data->nodepool;
  struct kd_point *points = my_data->points;
  size_t nPoints = my_data->nPoints;
  kdata_t *min = my_data->min;
  kdata_t *max = my_data->max;
  int depth = my_data->depth;
  int max_threads = my_data->max_threads;
  int dim = my_data->dim;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  int sortaxis = depth % dim;

  if (nPoints == 1) {
    if ((node = kd_allocNodeP(nodepool, points, 0, min, max, sortaxis, dim)) == NULL) return NULL;
    node->index = points[0].index;
    return node;
  }

  // If this iteration is allowed to start more threads, we first use them to parallelize the sorting 
  pmergesort(points, nPoints, sortaxis, max_threads);

  size_t pivot = nPoints / 2;
  if ((node = kd_allocNodeP(nodepool, points, pivot, min, max, sortaxis, dim)) == NULL) return NULL;

  memcpy(tmpMaxLeft, max, dim * sizeof(kdata_t));
  tmpMaxLeft[sortaxis] = node->location[sortaxis];
  kd_initArg(&argleft, nodepool, points, pivot, min, tmpMaxLeft,
             depth + 1, max_threads / 2, dim);

  if (max_threads > 1) {
    pthread_create(&threads[0], &attr, kd_doBuildTree, (void *) &argleft);
  } else {
    node->left = (kdNode *)kd_doBuildTree((void *) &argleft);
    if (!node->left) {
      kd_doDestroyTree(node);
      return NULL;
    }
  }

  memcpy(tmpMinRight, min, dim * sizeof(kdata_t));
  tmpMinRight[sortaxis] = node->location[sortaxis];
  kd_initArg(&argright, nodepool, &points[pivot], nPoints - pivot, tmpMinRight, max,
             depth + 1, max_threads / 2, dim);

  if (max_threads > 1) {
    pthread_create(&threads[1], &attr, kd_doBuildTree, (void *) &argright);
  } else {
    node->right = (kdNode *)kd_doBuildTree((void *) &argright);
    if (!node->right) {
      kd_doDestroyTree(node);
      return NULL;
    }
  }

  if (max_threads > 1) {
    pthread_join(threads[0], (void **) (&node->left));
    pthread_join(threads[1], (void **) (&node->right));
    if (!node->left || !node->right) {
      kd_doDestroyTree(node);
      return NULL;
    }
  }

  return (void *) node;
}

/* end of tree construction and destruction */

/* Functions dealing with result heaps */

/* Insert the sub-tree starting at node into the result heap res */
int
kd_insertResTree(struct kdNode *node, struct pqueue *res)
{
  if ( node == NULL ) return 1;

  if ( !kd_insertResTree(node->left, res) ) return 0;

  if ( kd_isleaf(node) )
    {
      struct resItem *point;
      if ((point = (struct resItem *)kd_malloc(sizeof(struct resItem), "kd_insertResTree: "))
          == NULL)
        return 0;

      point->node = node;
      point->dist_sq = -1;
      pqinsert(res, point);
    }
  
  if ( !kd_insertResTree(node->right, res) ) return 0;

  return 1;
}
