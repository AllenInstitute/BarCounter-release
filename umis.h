#ifndef UMIS_H
#define UMIS_H

#include <stdbool.h>

#include "barcodes.h"

// set the length of UMI
#define UMI_LEN 12

// set the first position of UMI in read1 sequences
#define UMI_FIRST 16

// define list_node struct for linked list of cell barcodes with UMI for per tag
typedef struct list_node {
    char barcode[BC_LEN + 1];
    struct list_node* next;
} list_node;

// define umi_node struct
typedef struct umi_node {
    bool exists;
    struct list_node** tag_lists;
    struct umi_node* children[4];
} umi_node;


// add a UMI sequene to a trie. Does not allow for 'N' bases.
// Each UMI leaf will track an array of pointers to linked lists of cell barcodes. This ensures that every combo of barcode/UMI/Tag is unique.
// If UMI is added: returns true. Else if UMI is NOT added, returns false.
bool add_umi(char *umi, umi_node* umi_root, int t_count, int t_index, char *cell);

// Unloads umi trie from memory. Returns true if successful, else returns false.
bool unload_umi_trie(umi_node* root, int t_count);

// helper function for unload_bc_trie
void unload_umi_helper(umi_node* unload_trav, int t_count);



#endif // UMIS_H