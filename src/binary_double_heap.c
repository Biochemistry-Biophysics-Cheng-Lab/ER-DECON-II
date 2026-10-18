#include "binary_heap.h"

void
adjust_max_double_heap(size_t root, size_t size, double* data)
{
    size_t child = root + root + 1;
    double root_value = data[root];

    while (child < size) {
	/* Compare the two children and take the larger. */
	if (child != size - 1 && data[child] < data[child + 1]) {
	    ++child;
	}
	if (root_value < data[child]) {
	    /*
	     * When the child is larger than the old root, descend the heap
	     * to find where the old root would go and shift the elements in
	     * the tree upward to make room.
	     */
	    data[(child - 1) / 2] = data[child];
	    child += child + 1;
	} else {
	    /* Found where the old root goes; break out. */
	    break;
	}
    }
    data[(child - 1) / 2] = root_value;
}


void
adjust_min_double_heap(size_t root, size_t size, double* data)
{
    size_t child = root + root + 1;
    double root_value = data[root];

    while (child < size) {
	/* Compare the two children and take the smaller. */
	if (child != size - 1 && data[child] > data[child + 1]) {
	    ++child;
	}
	if (root_value > data[child]) {
	    /*
	     * When the child is smaller than the old root, descend the heap
	     * to find where the old root would go and shift the elements in
	     * the tree upward to make room.
	     */
	    data[(child - 1) / 2] = data[child];
	    child += child + 1;
	} else {
	    /* Found where the old root goes; break out. */
	    break;
	}
    }
    data[(child - 1) / 2] = root_value;
}


int
is_max_double_heap(size_t root, size_t size, const double* data)
{
    size_t limit = (size - 1) / 2;
    size_t i;

    for (i = root; i < limit; ++i) {
	if (data[i] < data[i + i + 1]) {
	    return 0;
	}
	if (data[i] < data[i + i + 2]) {
	    return 0;
	}
    }

    return (root <= limit && limit + limit + 1 < size) ?
	data[limit] >= data[limit + limit + 1] : 1;
}


int
is_min_double_heap(size_t root, size_t size, const double* data)
{
    size_t limit = (size - 1) / 2;
    size_t i;

    for (i = root; i < limit; ++i) {
	if (data[i] > data[i + i + 1]) {
	    return 0;
	}
	if (data[i] > data[i + i + 2]) {
	    return 0;
	}
    }

    return (root <= limit && limit + limit + 1 < size) ?
	data[limit] <= data[limit + limit + 1] : 1;
}
