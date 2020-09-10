/*
Code written by Elliott Swanson of the Allen Institute for Immunology (elliott.swanson@alleninstitute.org). Free for academic use only.
See LICENSE for code reuse permissions.
*/

#ifndef BARCODES_H
#define BARCODES_H

#include <stdbool.h>

// set the length of 10X cell barcode
#define BC_LEN 16

// set the Q-score cutoff for low quality bases (PHRED - 33, ex. 53 == Q score of 20)
#define LOW_Q 53

// set the first position of barcodes in read1 sequences
#define BC_FIRST 0

// define bc_node struct for barcode trie. "counts" is a pointer to an array of int tag counts
typedef struct bc_node {
    bool exists;
    unsigned long int total;
    unsigned int *counts;
    struct bc_node* children[4];
}
bc_node;

// Loads a trie of barocdes of length BC_LEN into memory from plaintext whitelist file "input" using provided node* "root". Returns true if successful, else returns false.
// Allocates array of "t_count" unsigned short ints for cell tag counts.
bool load_bc_trie(const char *input, bc_node* root, int t_count);

// Loads a trie of barocdes of length BC_LEN into memory from a gzipped whitelist file "input" using provided node* "root". Returns true if successful, else returns false.
// Allocates array of "t_count" unsigned short ints for cell tag counts.
bool load_bc_trie_gzipped(const char *input, bc_node* root, int t_count);

// Returns a pointer to the bc_trie leaf bc_node corresponding to the input barcode. If the barocde doesn't exist in the trie, return NULL.
bc_node* get_bc_leaf(const char *seq, bc_node* root, int length);

// Unloads barcode trie from memory. Returns true if successful, else returns false.
bool unload_bc_trie(bc_node* root);

// helper function for unload_bc_trie
void unload_bc_helper(bc_node* unload_trav);



#endif // BARCODES_H