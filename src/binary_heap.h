#ifndef UCSF_MSG_BINARY_HEAP_H
#define UCSF_MSG_BINARY_HEAP_H

/*
 * A min heap has the property that if b is a child node of node a in the heap
 * then key(b) >= key(a).  A max heap has the property that if b is a child
 * node of node a in the heap then key(b) <= key(a).
 *
 * In this implementation, the nodes of the heap are stored as a complete
 * binary tree packed contiguously into a 1D array:  the children of the node
 * at index i are at the indices 2 * i + 1 and 2 * i + 2 (if 2 * i + 2 exceeds
 * the bounds of the array the node at i has at most one child; if 2 * i + 1
 * exceeds the bounds of the array the node at i has zero children (i.e. it is
 * a leaf)) and the parent of the node with index j is at the index
 * (j - 1) / 2.
 */

#include <stddef.h>  /* size_t */

/*
 * At entry, size must be positive, root must be in [0, size - 1), and the
 * subtrees rooted at 2 * root + 1 and 2 * root + 2 must satisfy the max heap
 * property (i.e. for j in [2 * root + 1, size), these conditions are true:
 *     (a) Either 2 * j + 1 >= size or data[j] >= data[2 * j + 1]
 *     (b) Either 2 * j + 2 >= size or data[j] >= data[2 * j + 2]
 * ).  At exit, the tree rooted at root satisfies the max heap property.
 */
void adjust_max_double_heap(size_t root, size_t size, double* data);

/*
 * At entry, size must be positive, root must be in [0, size - 1), and the
 * subtrees rooted at 2 * root + 1 and 2 * root + 2 must satisfy the min heap
 * property (i.e. for j in [2 * root + 1, size), these conditions are true:
 *     (a) Either 2 * j + 1 >= size or data[j] <= data[2 * j + 1]
 *     (b) Either 2 * j + 2 >= size or data[j] <= data[2 * j + 2]
 * ).  At exit, the tree rooted at root satisfies the min heap property.
 */
void adjust_min_double_heap(size_t root, size_t size, double* data);

/*
 * Returns one if the subtree rooted at root satisfies the max heap property;
 * otherwise returns zero.  size must be positive and root must be in
 * [0, size - 1).
 */
int is_max_double_heap(size_t root, size_t size, const double* data);

/*
 * Returns one if the subtree rooted at root satisfies the min heap property;
 * otherwise returns zero.  size must be positive and root must be in
 * [0, size - 1).
 */
int is_min_double_heap(size_t root, size_t size, const double* data);

#endif /* include guard */
