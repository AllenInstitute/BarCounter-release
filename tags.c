/*
Code written by Elliott Swanson of the Allen Institute for Immunology (elliott.swanson@alleninstitute.org). Free for academic use only.
See LICENSE for code reuse permissions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "tags.h"
#include "barcodes.h"

// calculate the hamming distance of two strings
int hamming_distance(char *str1, char *str2)
{
    int dist = 0;
    int m = strlen(str1);
    int n = strlen(str2);

    // ensure string lengths are equal
    if (m != n)
    {
        printf("string legnths of %s and %s are not equal: %i vs %i\n", str1, str2, m, n);
        return -1;
    }

    // calculate hamming distance
    for (int i = 0; i < m; i++)
    {
        if (str1[i] != str2[i])
        {
            dist++;
        }
    }
    return dist;
}

// Load CSV taglist (15bp "Tag", tag name) into array. Returns # of tags "t_count".
int load_taglist(char *taglist, char tags[MAX_TAGS][TAG_LEN + 1], char names[MAX_TAGS][NAME_LEN + 1])
{
    // read in csv file of taglist
    FILE *ptagl = fopen(taglist, "r");
    if (ptagl == NULL)
    {
        printf("Taglist %s failed to open\n", taglist);
        exit(4);
    } else {
        printf("Taglist %s opened successfully\n",taglist);
        }

    char *tag = NULL;
    char *name = NULL;
    char line[100];

    // var t_count to track # of tags, t_index to track array index of tag in tags array
    int t_count = 0;
    int t_index = 0;
    while (fscanf(ptagl, "%s", line) != EOF)
    {
        tag = strtok(line, ",");
        name = strtok(NULL, ",");

        // ensure correct length of tag
        if (strlen(tag) != TAG_LEN)
        {
            printf("Tag %s has length %li. All tag lengths must be exacly %i. Exiting...\n", tag, strlen(tag), TAG_LEN);
            exit(10);
        }

        if (t_count >= MAX_TAGS)
        {
            printf("Maximum of %i tags exceded! Exiting...\n", MAX_TAGS);
            exit(11);
        }
        // check if tag is in tags array
        t_index = in_tag_array(tags, tag, t_count);
        if (t_index == -1)
        {
            strncpy(tags[t_count], tag, TAG_LEN + 1);
        } else {
            printf("tag seq %s is listed multiple times in the taglist. Exiting...\n", tag);
            exit(12);
        }
        // ensure the length of the tag name does not exceed the maximum allowable # of characters, NAME_LEN
        if (strlen(name) > NAME_LEN)
        {
            printf("Tag name %s has length %li. The maximum allowable tag name lengths is %i. Exiting...\n", name, strlen(name), NAME_LEN);
            exit(13);
        }
        // check if name is in names array
        t_index = in_names_array(names, name, t_count);
        if (t_index == -1)
        {
            strncpy(names[t_count], name, NAME_LEN + 1);

        } else {
            printf("tag name %s is listed multiple times in the taglist. Exiting...\n", name);
            exit(14);
        }
        // if tag and name are successfully added to tags and names arrays, increment tag count "t_count"
        t_count++;
    }
    fclose(ptagl);

    return t_count;
}

/* Check tag array for tag string. Returns index of string if in array, otherwise returns -1.
    Array to be checked: "arr", string to check for: "s", length of array to check: "l" */
int in_tag_array(char tags[MAX_TAGS][TAG_LEN + 1], char *s, int l)
{
    for (int i = 0; i < l; i++)
    {
        if (strcmp(tags[i], s) == 0)
        {
            return i;
        }
    }
    return -1;
}

/* Check names array for tag name string. Returns index of string if in array, otherwise returns -1.
    Array to be checked: "arr", string to check for: "s", length of array to check: "l" */
int in_names_array(char names[MAX_TAGS][NAME_LEN + 1], char *s, int l)
{
    for (int i = 0; i < l; i++)
    {
        if (strcmp(names[i], s) == 0)
        {
            return i;
        }
    }
    return -1;
}

// Check Taglist for hamming dist to ensure that each tag has a hamming dist of >= MIN_TAG_HDIST to every other tag
void check_tag_dist(char tags[MAX_TAGS][TAG_LEN + 1], int t_count)
{
    int dist;
    for (int i = 0; i < t_count - 1; i++)
    {
        for (int j = i+1; j < t_count; j++)
        {
            dist = hamming_distance(tags[i], tags[j]);
            if (dist < MIN_TAG_HDIST)
            {
                printf("Hamming distance between tags %s and %s is %i. The minimum allowed is %i.\nExiting...\n", tags[i], tags[j], dist, MIN_TAG_HDIST);
                exit(16);
            }
        }
    }
}

// load a trie of every possible tag seq with 1 hamming distance from a tag in the taglist of length "t_count".
bool load_tag_trie(char tags[MAX_TAGS][TAG_LEN + 1], tag_node* tag_root, int t_count)
{
    char tag[TAG_LEN + 1];

    // set variables for adding seqs + 1 hamming dist from tag. "next" is first base to test (index of bases), "temp" stored temp tag string. Includes 'N' to allow for a single 'N' in tag sequence.
    char bases[5] = "NACGT";
    char temp[TAG_LEN + 1];
    int next;

    // loop through all tags in taglist
    for (int t = 0; t < t_count; t++)
    {
        // update tag to check
        // tag = tags[t];
        strcpy(tag, tags[t]);

        if (!add_tag(tag, tag_root, t))
        {
            printf("failed to add tag %s\n", tag);
            return false;
        }

        // change base at every position and add to trie
        for (int m = 0; m < TAG_LEN; m++)
        {
            // reset temp to equal tag
            strcpy(temp, tag);

            // add tag string with each of the three remaining bases + 'N' substituted at postion m
            switch(temp[m])
            {
                case 'A': next = 2; break;
                case 'C': next = 3; break;
                case 'G': next = 4; break;
                case 'T': next = 5; break;
            }
            for (int r = 0; r < 4; r++)
            {
                temp[m] = bases[(next+r)%5];
                if (!add_tag(temp, tag_root, t))
                {
                    printf("failed to add tag %s\n", tag);
                    return false;
                }
            }
        }
    }
    return true;
}

// add a tag sequene to a trie and update "index" to equal the tag index in "tags" and "names" array. Including 'N' to allow for a single low quality basecall during sequencing.
bool add_tag(char *tag, tag_node* tag_root, int t)
{
    int i;
    // declare and initialize travelling node pointer to NULL
    tag_node* trav = tag_root;

    for (int c = 0; c < TAG_LEN; c++)
    {
        // convert current character to int i for selectng the correct child node
        switch(tag[c])
        {
            case 'A': i = 0; break;
            case 'C': i = 1; break;
            case 'G': i = 2; break;
            case 'T': i = 3; break;
            case 'N': i = 4; break;
            default: printf("Non DNA base included in taglist tag %s. Exiting...\n", tag);
                     exit(17);
        }
        //check the value at children[i]. If child doesn't exist, create child node move trav
        if (trav->children[i] == NULL)
        {
            trav->children[i] = calloc(1, sizeof(tag_node));
        }
        trav = trav->children[i];
    }
    trav->exists = true;
    trav->index = t;

    return true;
}

// Check tag_trie for tag seq. If present, returns tag index. Else, returns -1.
int get_tag_index(const char *tag, tag_node* root)
{
    tag_node *m_trav = root;
    int i;
    for (int c = 0; c < TAG_LEN; c++)
    {
        // convert current character to int i for selectng the correct child node
        switch(tag[c])
        {
            case 'A': i = 0; break;
            case 'C': i = 1; break;
            case 'G': i = 2; break;
            case 'T': i = 3; break;
            case 'N': i = 4; break;
            default: printf("Non DNA base included in tag %s from input FastQ. Exiting...\n", tag);
                     exit(25);
        }
        //check the value at children[i]
        if (m_trav->children[i] == NULL)
        {
            return -1;
        }
        m_trav = m_trav->children[i];
    }
    if (m_trav->exists == true)
    {
        return m_trav->index;
    }
    return -1;
}

// Unloads tag trie from memory. Returns true if successful, else returns false.
bool unload_tag_trie(tag_node *root)
{
    tag_node *unload_trav = root;
    unload_tag_helper(unload_trav);
    return true;
}

// helper function for unload_tag_trie
void unload_tag_helper(tag_node* unload_trav)
{
    for (int i=0; i<5; i++)
    {
        if (unload_trav->children[i] != NULL)
        {
            unload_tag_helper(unload_trav->children[i]);
        }
    }

    free(unload_trav);
}

