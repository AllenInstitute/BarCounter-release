/*
Code written by Elliott Swanson of the Allen Institute for Immunology (elliott.swanson@alleninstitute.org). Free for academic use only.
See LICENSE for code reuse permissions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "umis.h"
#include "barcodes.h"

// add a UMI sequene to a trie. Does not allow for 'N' bases.
// Each UMI leaf will track an array of pointers to linked lists of cell barcodes. This ensures that every combo of barcode/UMI/Tag is unique.
// If UMI is added: returns true. Else if UMI is NOT added, returns false.
bool add_umi(char *umi, umi_node* umi_root, int t_count, int t_index, char *cell)
{
    int i;
    // declare and initialize travelling node pointer to NULL
    umi_node* trav = umi_root;

    for (int c = 0; c < UMI_LEN; c++)
    {
        // convert current character to int i for selectng the correct child node
        switch(umi[c])
        {
            case 'A': i = 0; break;
            case 'C': i = 1; break;
            case 'G': i = 2; break;
            case 'T': i = 3; break;
            case 'N': return false;
            default: printf("Non DNA base included in UMI %s. Exiting...\n", umi);
                     exit(26);
        }
        //check the value at children[i]. If child doesn't exist, create child node move trav
        if (trav->children[i] == NULL)
        {
            trav->children[i] = calloc(1, sizeof(umi_node));
        }
        trav = trav->children[i];
    }
    trav->exists = true;

    // allocate array of t_count linked lists of cell barcodes
    if (trav->tag_lists == NULL)
    {
        trav->tag_lists = calloc(1, sizeof(list_node*) * t_count);
    }

    // add cell barcode to linked list for specific tag
    if (trav->tag_lists[t_index] == NULL)
    {
        trav->tag_lists[t_index] = calloc(1, sizeof(list_node));
    }
    // initialize travelling list pointer to beginning of correct tag list
    list_node* l_trav = trav->tag_lists[t_index];

    // check linked list for cell barcode. If cell barcode is already listed: return
    while (true)
    {
        // node has no "cell" barcode
        if (strlen(l_trav->barcode) == 0)
        {
            // update list node with cell barcode
            strcpy(l_trav->barcode, cell);
            return true;
        }
        // if "cell" is already in the list: no need to add, break loop
        else if (strcmp(l_trav->barcode, cell) == 0)
        {
            return false;
        }
        // if match is not found, proceed to next node. Allocate "next" if necessary.
        if (l_trav->next == NULL)
        {
            l_trav->next = calloc(1, sizeof(list_node));
        }
        // update l_trav to move to next node
        l_trav = l_trav->next;
    }
    return true;
}

// Unloads umi trie from memory. Returns true if successful, else returns false.
bool unload_umi_trie(umi_node *root, int t_count)
{
    umi_node *unload_trav = root;
    unload_umi_helper(unload_trav, t_count);
    return true;
}

// helper function for unload_bc_trie
void unload_umi_helper(umi_node* unload_trav, int t_count)
{
    // temp variable to track current list node
    list_node* temp = NULL;
    list_node* p_next = NULL;

    for (int i=0; i<4; i++)
    {
        if (unload_trav->children[i] != NULL)
        {
            unload_umi_helper(unload_trav->children[i], t_count);
        }
    }
    // if in leaf node: all linked lists that were allocated via add_umi
    if (unload_trav->tag_lists != NULL)
    {
        for (int j = 0; j < t_count; j++)
        {
            if (unload_trav->tag_lists[j] != NULL)
            {
                temp = unload_trav->tag_lists[j];
                while (true)
                {
                    // if not at end of list: track pointer to next node, free current node, move to next node
                    if (temp->next != NULL)
                    {
                        p_next = temp->next;
                        free(temp);
                        temp = p_next;
                    }
                    else {
                        free(temp);
                        break;
                    }
                }
            }
        }
        free(unload_trav->tag_lists);
    }
    free(unload_trav);
}