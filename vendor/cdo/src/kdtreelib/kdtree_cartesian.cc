/*!\file kdtree_cartesian.c
 * \brief Routines specific to the n-dimensional cartesian kd-tree
 *
 */
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include "kdtree.h"
#include "pqueue.h"

/* ********************************************************************

   general utility functions 

   ********************************************************************* */

static constexpr
kdata_t square(const kdata_t x)  noexcept
{
  return x*x;
}

static constexpr
kdata_t kd_dist_sq(const kdata_t *restrict a, const kdata_t *restrict b)  noexcept
{
  return (float)(square((a[0]-b[0]))+square((a[1]-b[1]))+square((a[2]-b[2])));
}


kdata_t kd_min(kdata_t x, kdata_t y)
{
  return x < y ? x : y;
}


/* ******************************************************************
   
   Functions for building and destroying trees 

   ******************************************************************** */


/*! \brief build kd-tree structure
 *
 * \param points an array of kd_points (struct with position vector
 * and data container).
 *
 * \param nPoints the length of the points array.
 *
 * \param constr a pointer to a void *constructor() function to include
 * the data container in the tree; optional, can be NULL
 *
 * \param destr a pointer to a void destructor() function to free()
 * the data containers in the tree; optional, can be NULL, but should
 * be given if the constr argument is non-NULL.
 *
 * \param min a vector with the minimum positions of the corners of the
 * hyperrectangle containing the data.
 *
 * \param max a vector with the maximum positions of the corners of
 * the hyperrectangle containing the data.
 *
 * \param dim the dimensionality of the data.
 *
 * \param max_threads the maximal number of threads spawned for
 * construction of the tree. The threads will be unbalanced if this is
 * not a power of 2.
 *
 * \return root node of the tree
 */
kdTree_t *
kd_buildTree(struct kd_point *points, size_t nPoints,
             kdata_t *min, kdata_t *max, int dim, int max_threads)
{
  kdTree_t *tree = (kdTree_t *) calloc(1, sizeof(kdTree_t));
#ifdef MEMPOOL
  tree->nodepool = kd_nodepool_new(nPoints*2);
#endif
  struct kd_thread_data my_data;
  kd_initArg(&my_data, tree->nodepool, points, nPoints, min, max, 0, max_threads, dim);
  tree->node = (kdNode *)kd_doBuildTree(&my_data);
  if ( tree->node == NULL )
    {
      free(tree);
      tree = NULL;
    }
  
  return tree;
}


/* *******************************************************************

   Functions for range searches 
  
   *******************************************************************
*/

/* Returns 1 if node is a point in the hyperrectangle defined by
   minimum and maximum vectors min and max. */
int
kd_isPointInRect(struct kdNode *node, kdata_t *min, kdata_t *max, int dim)
{
    int i;

    if (node == NULL)
        return 0;

    for(i = 0; i < dim; i++) {
        if (node->location[i] < min[i] || node->location[i] > max[i])
            return 0;
    }
    return 1;
}

/* Returns 1 if the hyperrectangle of node is fully contained within
   the HR described by the minimum and maximum vectors min and
   max. Returns 0 otherwise. */
int
kd_isRectInRect(struct kdNode *node, kdata_t *min, kdata_t *max, int dim)
{
    int i;

    if (node == NULL)
        return 0;

    for(i = 0; i < dim; i++) {
        if (node->min[i] < min[i] || node->max[i] > max[i])
            return 0;
    }
    return 1;
}

/* Returns 1 if the hyperrectangle of node overlaps the HR described
   the minimum and maximum vectors min and max. Returns 0
   otherwise. */
int
kd_rectOverlapsRect(struct kdNode *node, kdata_t *min, kdata_t *max, int dim)
{
    int i;

    if (node == NULL)
        return 0;

    for(i = 0; i < dim; i++) {
        if ((node->min[i] > max[i] || node->max[i] < min[i]))
            return 0;
    }
    return 1;
}

/*!
 * \brief Perform orthogonal range search (get all points in a
 * hyperrectangle).
 *
 * \param node the root node of tree to be searched.
 *
 * \param min a vector with the minimum positions of the corners of the
 * hyperrectangle containing the data.
 *
 * \param max a vector with the maximum positions of the corners of
 * the hyperrectangle containing the data.
 *
 * \param dim the dimension of the data.
 *
 * \return  Pointer to a priority queue, NULL in case of problems.
*/
struct pqueue *
kd_ortRangeSearch(struct kdNode *node, kdata_t *min, kdata_t *max, int dim)
{
    struct pqueue *res;

    if ((res = pqinit(NULL, 1)) == NULL)
        return NULL;
    if (!kd_doOrtRangeSearch(node, min, max, dim, res)) {
        for(size_t i = 0; i < res->size; i++) {
            free(res->d[i]);
        }
        free(res->d);
        free(res);
        return NULL;
    }
    return res;
}

/* This is the orthogonal range search. Returns 1 if okay, 0 in case
   of problems. */
int
kd_doOrtRangeSearch(struct kdNode *node, kdata_t *min, kdata_t *max,
                    int dim, struct pqueue *res)
{

    if (node == NULL)
        return 1;

    if (kd_isleaf(node) && kd_isPointInRect(node, min, max, dim)) {
        return kd_insertResTree(node, res);
    } else {
        if (kd_isRectInRect(node->left, min, max, dim)) {
            if (!kd_insertResTree(node->left, res))
                return 0;
        } else {
            if (kd_rectOverlapsRect(node->left, min, max, dim))
                if (!kd_doOrtRangeSearch(node->left, min, max, dim, res))
                    return 0;
        }
        if (kd_isRectInRect(node->right, min, max, dim)) {
            if (!kd_insertResTree(node->right, res))
                return 0;
        } else {
            if (kd_rectOverlapsRect(node->right, min, max, dim))
                if (!kd_doOrtRangeSearch(node->right, min, max, dim, res))
                    return 0;
        }
    }
    return 1;
}

/*! 
 * \brief Find the nearest neighbor of a point.
 *
 * \param node the root node of the tree to be searched.
 *
 * \param p a vector to the point whose nearest neighbor is sought.
 *
 * \param max_dist_sq the square of the maximum distance to the
 * nearest neighbor.
 *
 * \param dim the dimension of the data.
 *
 * \return A pointer to node containing the nearest neighbor.
 * max_dist_sq is set to the square of the distance to the nearest
 * neigbor.
 */
struct kdNode *
kd_nearest(struct kdNode *node, kdata_t *p, kdata_t *max_dist_sq, int dim)
{
    if ( !node ) return NULL;

    struct kdNode *nearer, *further;
    if ( p[node->split] < node->location[node->split] ) {
        nearer  = node->left;
        further = node->right;
    } else {
        nearer  = node->right;
        further = node->left;
    }
    
    struct kdNode *tmp = kd_nearest(nearer, p, max_dist_sq, dim);
    struct kdNode *nearest = NULL;
    if ( tmp )
        nearest = tmp;
    else
        nearest = node;

    kdata_t dist_sq = kd_dist_sq(nearest->location, p);
    if (*max_dist_sq > dist_sq)
	*max_dist_sq = dist_sq;

    if (!further)
        return nearest;

    kdata_t dx = kd_min(KDATA_ABS(p[node->split] - further->min[node->split]),
                        KDATA_ABS(p[node->split] - further->max[node->split]));
    if ( *max_dist_sq > kd_sqr(dx) ) {
        /*
         * some part of the further hyper-rectangle is in the search
         * radius, search the further node 
         */
        struct kdNode *tmp_nearest;
        tmp = kd_nearest(further, p, max_dist_sq, dim);
        if ( tmp )
            tmp_nearest = tmp;
        else
            tmp_nearest = further;

        kdata_t tmp_dist_sq = kd_dist_sq(tmp_nearest->location, p);
        if ( tmp_dist_sq < dist_sq || (tmp_dist_sq <= dist_sq && tmp_nearest->index < nearest->index) ) {
            nearest = tmp_nearest;
            dist_sq = tmp_dist_sq;
            *max_dist_sq = kd_min(dist_sq, *max_dist_sq);
        }
    }
    return nearest;
}

/*! 
 * \brief Return the q nearest-neighbors to a point.
 *
 * \param node the root node of the tree to be searched.
 *
 * \param p a vector to the point whose nearest neighbors are sought.
 *
 * \param max_dist_sq the square of the maximum distance to the
 * nearest neighbors.
 *
 * \param q the maximum number of points to be retured.
 *
 * \param dim the dimension of the data.
 *
 * \return A pointer to a priority queue of the points found, or NULL
 * in case of problems.
 */
struct pqueue *
kd_qnearest(struct kdNode *node, kdata_t *p,
            kdata_t *max_dist_sq, size_t q, int dim)
{
    struct pqueue *res = pqinit(NULL, q + 2);
    if ( res == NULL) return NULL;
    
    if ( !kd_doQnearest(node, p, max_dist_sq, q + 1, dim, res) ) {
        for ( size_t i = 0; i < res->size; ++i ) free(res->d[i]);
        free(res->d);
        free(res);
        return NULL;
    }
    
    return res;
}

/* *
 * This is the q nearest-neighbor search.
 *
 * This is a modification of the range search in which the maximum
 * search radius is decreased to the maximum of the queue as soon as
 * the queue is filled.
 *
 * return 1 if okay, zero in case of problems
 */
// Uwe Schulzweida: extract kd_check_dist() from kd_doQnearest()
static bool
kd_check_dist(struct kdNode *node, kdata_t *p,
              kdata_t *max_dist_sq, size_t q, struct pqueue *res)
{
  kdata_t dist_sq = kd_dist_sq(node->location, p);
  if ( dist_sq < *max_dist_sq && kd_isleaf(node) )
    {
      struct resItem *point = (struct resItem *) kd_malloc(sizeof(struct resItem), "kd_doQnearest: ");
      if ( point == NULL) return false;
      point->node = node;
      point->dist_sq = dist_sq;
      pqinsert(res, point);
    }

  if ( res->size > q )
    {
      struct resItem *item;
      pqremove_max(res, &item);
      free(item);
      if ( res->size > 1 )
        {
          // Only inspect the queue if there are items left 
          pqpeek_max(res, &item);
          *max_dist_sq = item->dist_sq;
        }
      else
        {
          // Nothing was found within the max search radius 
          *max_dist_sq = 0;
        }
    }

  return true;
}

int
kd_doQnearest(struct kdNode *node, kdata_t *p,
              kdata_t *max_dist_sq, size_t q, int dim, struct pqueue *res)
{
  if ( !node ) return 1;

  if ( !kd_check_dist(node, p, max_dist_sq, q, res) ) return 0;

  struct kdNode *nearer, *further;
  if ( p[node->split] < node->location[node->split] ) {
    nearer  = node->left;
    further = node->right;
  } else {
    nearer  = node->right;
    further = node->left;
  }
  if ( !kd_doQnearest(nearer, p, max_dist_sq, q, dim, res) ) return 0;

  if ( !further ) return 1;

  kdata_t dx = kd_min(KDATA_ABS(p[node->split] - further->min[node->split]),
                      KDATA_ABS(p[node->split] - further->max[node->split]));
    
  if ( *max_dist_sq > kd_sqr(dx) ) {
    /*
     * some part of the further hyper-rectangle is in the search
     * radius, search the further node 
     */
    if (!kd_doQnearest(further, p, max_dist_sq, q, dim, res)) return 0;

    if (!kd_check_dist(node, p, max_dist_sq, q, res)) return 0;
  }

  return 1;
}


/*!
 * \brief Perform a range search around a point.
 *
 * \param node the root node of the tree to be searched.
 *
 * \param p the location of the point around which the search is carried out .
 *
 * \param max_dist_sq the square of the radius of the hypersphere.
 *
 * \param dim the dimension of the data.  \param ordered determines
 * whether the result list should be ordered in increasing distance
 * (KD_ORDERED) or unordered (KD_UNORDERED).
 *
 * \return A pointer to a priority queue containing the points found,
 * NULL in case of problems.
 */
struct pqueue *
kd_range(struct kdNode *node, kdata_t *p, kdata_t *max_dist_sq,
         int dim, int ordered)
{
    struct pqueue *res = pqinit(NULL, 1);
    if ( res == NULL ) return NULL;
    
    if ( !kd_doRange(node, p, max_dist_sq, dim, res, ordered) ) {
        for( size_t i = 0; i < res->size; ++i ) {
            free(res->d[i]);
        }
        free(res->d);
        free(res);
        return NULL;
    }
    return res;
}


/* This is the range search. Returns 1 if okay, 0 in case of problems */
int
kd_doRange(struct kdNode *node, kdata_t *p, kdata_t *max_dist_sq,
           int dim, struct pqueue *res, int ordered)
{

    struct kdNode *nearer, *further;
    struct resItem *point;
    kdata_t dist_sq, dx;

    if (!node)
        return 1;

    dist_sq = kd_dist_sq(node->location, p);
    if (dist_sq < *max_dist_sq && kd_isleaf(node)) {
        if ((point = (struct resItem *)kd_malloc(sizeof(struct resItem), "kd_doRange:"))
            == NULL)
            return 0;
        point->node = node;
        point->dist_sq = ordered ? dist_sq : -1;
        pqinsert(res, point);
    }

    if (p[node->split] < node->location[node->split]) {
        nearer = node->left;
        further = node->right;
    } else {
        nearer = node->right;
        further = node->left;
    }
    if (!kd_doRange(nearer, p, max_dist_sq, dim, res, ordered))
        return 0;

    if (!further)
        return 1;

    dx = kd_min(KDATA_ABS(p[node->split] - further->min[node->split]),
                KDATA_ABS(p[node->split] - further->max[node->split]));
    if (*max_dist_sq > kd_sqr(dx)) {
        /*
         * some part of the further hyper-rectangle is in the search
         * radius, search the further node 
         */
        if (!kd_doRange(further, p, max_dist_sq, dim, res, ordered))
            return 0;
        dist_sq = kd_dist_sq(node->location, p);

        if (dist_sq < *max_dist_sq && kd_isleaf(node)) {
            if ((point = (struct resItem *)kd_malloc(sizeof(struct resItem), "kd_doRange: "))
                == NULL)
                return 0;
            point->node = node;
            point->dist_sq = ordered ? dist_sq : -1;
            pqinsert(res, point);
        }

    }
    return 1;
}

/* End range searching functions */
