#ifndef TAGS_H
#define TAGS_H

#include <stdbool.h>

#include "barcodes.h"

// set the maximum number of antibody tags
#define MAX_TAGS 300

// set the tag sequence length
#define TAG_LEN 15

// set the max allowable tag name length
#define NAME_LEN 50

// set the minimum acceptable hamming distance between two tag sequences
#define MIN_TAG_HDIST 3

// set the first position of antibody tag in read2 sequences
#define TAG_FIRST 0


// define tag_node struct
typedef struct tag_node {
    bool exists;
    int index;
    struct tag_node* children[5];
} tag_node;


// calculate the hamming distance of two strings
int hamming_distance(char *str1, char *str2);

// Load CSV taglist (15bp "Tag", tag name) into array. Verify a min hamming dist of MIN_TAG_HDIST for all combinations of tags.
int load_taglist(char *taglist, char tags[MAX_TAGS][TAG_LEN + 1], char names[MAX_TAGS][NAME_LEN + 1]);

/* Check tag array for tag string. Returns index of string if in array, otherwise returns -1.
    Array to be checked: "arr", string to check for: "s", length of array to check: "l" */
int in_tag_array(char tags[MAX_TAGS][TAG_LEN + 1], char *s, int l);

/* Check names array for tag name string. Returns index of string if in array, otherwise returns -1.
    Array to be checked: "arr", string to check for: "s", length of array to check: "l" */
int in_names_array(char names[MAX_TAGS][NAME_LEN + 1], char *s, int l);

// Check Taglist for hamming dist to ensure that each tag has a hamming dist of >= MIN_TAG_HDIST to every other tag
void check_tag_dist(char tags[MAX_TAGS][TAG_LEN + 1], int t_count);

// load a trie of every possible tag seq with 1 hamming distance from a tag in the taglist of length "t_count".
bool load_tag_trie(char tags[MAX_TAGS][TAG_LEN + 1], tag_node* tag_root, int t_count);

// add a tag sequene to a trie and update "index" to equal the tag index in "tags" and "names" array. Including 'N' to allow for a single low quality basecall during sequencing.
bool add_tag(char *tag, tag_node* tag_root, int t);

// Check tag_trie for tag seq. If present, returns tag index. Else, returns -1.
int get_tag_index(const char *tag, tag_node* root);

// Unloads tag trie from memory. Returns true if successful, else returns false.
bool unload_tag_trie(tag_node *root);

// helper function for unload_tag_trie
void unload_tag_helper(tag_node* unload_trav);


#endif // TAGS_H