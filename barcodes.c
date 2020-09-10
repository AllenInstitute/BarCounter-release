/*
Code written by Elliott Swanson of the Allen Institute for Immunology (elliott.swanson@alleninstitute.org). Free for academic use only.
See LICENSE for code reuse permissions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <zlib.h>

#include "barcodes.h"

// Loads a trie of barocdes of length BC_LEN into memory from plaintext whitelist file "input" using provided node* "root". Returns true if successful, else returns false.
// Allocates array of "t_count" unsigned ints for cell tag counts.
bool load_bc_trie(const char *input, bc_node* root, int t_count)
{
    FILE *fp = fopen(input, "r");
    if (fp == NULL)
    {
        printf("%s could not be opened. Exiting...\n", input);
        exit(14);
    }

    // declare and initialize travelling node pointer to NULL
    bc_node* trav = NULL;

    int i;
    char barcode[20];

    while (fscanf(fp, "%s", barcode) != EOF){
        // reset trav pointer to root
        trav = root;
        // ensure barcodes in whitelist are the length of BC_LEN
        if (strlen(barcode) != BC_LEN)
        {
            printf("Barcode length of %li for %s is invalid. Length must be %i bases long.\n",strlen(barcode), barcode, BC_LEN);
            return false;
        }
        for (int c = 0; c < BC_LEN; c++){
            // convert current character to int i for selectng the correct child node
            switch(barcode[c])
            {
                case 'A': i = 0; break;
                case 'C': i = 1; break;
                case 'G': i = 2; break;
                case 'T': i = 3; break;
                default: printf("Non DNA base included in whitelist barcode %s. Exiting...\n", barcode);
                         exit(15);
            }
            //check the value at children[i]. If child doesn't exist, create child node move trav
            if (trav->children[i] == NULL)
            {
                trav->children[i] = calloc(1, sizeof(bc_node));
            }
            trav = trav->children[i];
        }
        trav->exists = true;

        // allocate counts array and initialize all values to 0.
        if (trav->counts == NULL) // ensure array is only allocated once for each barcode
        {
            trav->counts = malloc(sizeof(unsigned int) * t_count);
            for (int t = 0; t < t_count; t++)
            {
                trav->counts[t] = 0;
            }
            trav->total = 0;
        }
    }
    fclose(fp);
    printf("Barcode whitelist %s loaded successfully\n",input);
    return true;
}

// Loads a trie of barocdes of length BC_LEN into memory from a gzipped whitelist file "input" using provided node* "root". Returns true if successful, else returns false.
// Allocates array of "t_count" unsigned ints for cell tag counts.
bool load_bc_trie_gzipped(const char *input, bc_node* root, int t_count)
{
    gzFile fp = gzopen(input, "r");
    if (fp == NULL)
    {
        printf("%s could not be opened. Exiting...\n", input);
        exit(14);
    }

    // declare and initialize travelling node pointer to NULL
    bc_node* trav = NULL;

    int i;
    char barcode[20];

    while (true){

        // read line of whitelist into 'barcode'
        gzgets(fp, barcode, 20);

        // reset trav pointer to root
        trav = root;
        // ensure barcodes in whitelist are the length of BC_LEN
        if ((strlen(barcode) != BC_LEN) && (barcode[16] != '\n'))
        {
            printf("Barcode length of %li for %s is invalid. Length must be %i bases long.\n",strlen(barcode) - 1, barcode, BC_LEN);
            return false;
        }
        for (int c = 0; c < BC_LEN; c++){
            // convert current character to int i for selectng the correct child node
            switch(barcode[c])
            {
                case 'A': i = 0; break;
                case 'C': i = 1; break;
                case 'G': i = 2; break;
                case 'T': i = 3; break;
                default: printf("Non DNA base included in whitelist barcode %s. Exiting...\n", barcode);
                         exit(15);
            }
            //check the value at children[i]. If child doesn't exist, create child node move trav
            if (trav->children[i] == NULL)
            {
                trav->children[i] = calloc(1, sizeof(bc_node));
            }
            trav = trav->children[i];
        }
        trav->exists = true;

        // allocate counts array and initialize all values to 0.
        if (trav->counts == NULL) // ensure array is only allocated once for each barcode
        {
            trav->counts = malloc(sizeof(unsigned int) * t_count);
            for (int t = 0; t < t_count; t++)
            {
                trav->counts[t] = 0;
            }
            trav->total = 0;
        }

        // break when the end of the whitelist is reached
        if (gzeof(fp))
        {
            break;
        }
    }
    gzclose(fp);
    printf("Barcode whitelist %s loaded successfully\n",input);

    return true;
}

// Returns a pointer to the bc_trie leaf bc_node corresponding to the input barcode. If the barocde doesn't exist in the trie, return NULL.
bc_node* get_bc_leaf(const char *seq, bc_node* root, int length)
{
    bc_node* m_trav = root;
    int i;

    for (int c = 0; c < length; c++)
    {
        // convert current character to int i for selectng the correct child node
        switch(seq[c])
        {
            case 'A': i = 0; break;
            case 'C': i = 1; break;
            case 'G': i = 2; break;
            case 'T': i = 3; break;
            case 'N': return NULL;
            default: printf("Non DNA base included in barcode %s from input FastQ. Exiting...\n", seq);
                     exit(19);
        }
        //check the value at children[i]
        if (m_trav->children[i] == NULL)
        {
            return NULL;
        }
        m_trav = m_trav->children[i];
    }
    if (m_trav->exists == true)
    {
        return m_trav;
    }
    return NULL;
}

// Unloads barcode trie from memory. Returns true if successful, else returns false.
bool unload_bc_trie(bc_node *root)
{
    bc_node *unload_trav = root;
    unload_bc_helper(unload_trav);
    return true;
}

// helper function for unload_bc_trie
void unload_bc_helper(bc_node* unload_trav)
{
    for (int i=0; i<4; i++)
    {
        if (unload_trav->children[i] != NULL)
        {
            unload_bc_helper(unload_trav->children[i]);
        }
    }
    // if in leaf node: free counts array that was allocated in load barcode trie
    if (unload_trav->counts != NULL)
    {
        free(unload_trav->counts);
    }

    free(unload_trav);
}