/*
Code written by Elliott Swanson of the Allen Institute for Immunology (elliott.swanson@alleninstitute.org). Free for academic use only.
See LICENSE for code reuse permissions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <libgen.h>

#include "barcodes.h"
#include "tags.h"
#include "umis.h"

#define MAX_FASTQ 100

// return string f_time with formatted current GMT (UTC)
char* get_datetime(char* f_time);

int main(int argc, char *argv[])
{
    // format usage string
    char *command = "./barcounter -w {barcode whitelist} -t {taglist} -1 {read1 fastqs} -2 {read2 fastqs} -o {output directory}";
    char *summary = "BarCounter counts the number of valid read2 antibody derived tags (ADTs) that match tags in the user provided taglist.\nTag counts are generated for each read1 cell barcode that is present in the user provided whitelist.\nTags will be counted once per read1 Unique Molecular Identifier (UMI).";
    char *description = "-w whitelist: list of valid cell barcodes (one per line) in .txt or .gz format\n-t taglist: list of valid ADTs and their names in .csv format (sequence,name)\n-1 read1: gzipped files in fastq format, comma separated file list with no spaces\n-2 read2: gzipped files in fastq format, comma separated file list with no spaces\n-o output directory: if the directory does not yet exist BarCounter will create it. All outputs will be created in this location.";
    char usage[1000];
    snprintf(usage, 1000, "%s\n\n%s\n\n%s\n", command, summary, description);

    // initialize two dimensional array for tag sequences
    char tags[MAX_TAGS][TAG_LEN + 1];
    // initialize two dimensional array for tag names
    char names[MAX_TAGS][NAME_LEN + 1];

    // track global variable for report
    unsigned long long int total_reads = 0;
    unsigned long long int valid_barcodes = 0;
    unsigned long long int corrected_barcodes = 0;
    unsigned long long int valid_tags = 0;

    // Verify command line arguments. If usage is incorrect print Usage and exit with code 1. Help option -h prints usage and exits the program.
    int a;
    char *read1 = NULL, *read2 = NULL, *whitelist = NULL, *taglist = NULL, *outdir = NULL;
    bool help = false;
    while ((a = getopt(argc, argv, "1:2:w:t:o:h")) != -1)
    {
        switch(a)
        {
            case '1': read1 = optarg; break;
            case '2': read2 = optarg; break;
            case 't': taglist = optarg; break;
            case 'w': whitelist = optarg; break;
            case 'o': outdir = optarg; break;
            case 'h': help = true; break;
        }
    }

    if (help == true)
    {
        printf("%s", usage);
        exit(0);
    }

    // format outpur dir to ensure it ends in trailing '/'. If option is omitted defailt to './'
    if (outdir == NULL || strlen(outdir) == 0)
    {
        printf("output directory must be provided using -o\n\n%s", usage);
        exit(1);
    }
    else if (outdir[strlen(outdir)-1] != '/')
    {
        int out_len = strlen(outdir);
        outdir[out_len] = '/';
        outdir[out_len+1] = '\0';
    }

    // ensure all required args are provided
    if (read1 == NULL || read2 == NULL || taglist == NULL || whitelist == NULL)
    {
        printf("Required argument is missing. Refer to Usage below:\n\n%s\n",usage);
        exit(1);
    }

    // process all read1 fastq paths
    // read each comma delimited path into a variable
    char** paths1 = malloc(sizeof(char *) * MAX_FASTQ);
    char* name_token = NULL;
    int read1_count = 0;

    name_token = strtok(read1, ",");
    // assign each path to pointer in char ** array "paths1"
    for (int v = 0; v < MAX_FASTQ; v++)
    {
        if (name_token == NULL)
        {
            break;
        }
        paths1[v] = name_token;
        read1_count++;
        name_token = strtok(NULL, ",");
    }
    // ensure that fastq read pair input is limited to MAX_FASTQ number
    if (strtok(NULL, ",") != NULL)
    {
        printf("Maximum number of fastq pairs %i exceeded. Exiting...\n", MAX_FASTQ);
        exit(2);
    }

    // process all read2 fastq paths
     // read each comma delimited path into a variable
    char** paths2 = malloc(sizeof(char *) * MAX_FASTQ);
    name_token = NULL;
    int read2_count = 0;
    name_token = strtok(read2, ",");
    // assign each path to pointer in char ** array "paths2"
    for (int o = 0; o < MAX_FASTQ; o++)
    {
        if (name_token == NULL)
        {
            break;
        }
        paths2[o] = name_token;
        read2_count++;
        name_token = strtok(NULL, ",");
    }
    // ensure that fastq read pair input is limited to MAX_FASTQ number
    if (strtok(NULL, ",") != NULL)
    {
        printf("Maximum number of fastq pairs %i exceeded. Exiting...\n", MAX_FASTQ);
        exit(2);
    }

    // verify that the same number of read1 and read2 fastq files were provided
    if (read1_count != read2_count)
    {
        printf("The number of read1 and read2 fastq files are not equal. %i read1 files and %i read2 files were provided. Exiting...\n", read1_count, read2_count);
        exit(3);
    }

    // verify that all fastq read pairs have the same sample name.
    // It is assumed that standard Illumina fastq naming convention (SampleName_S1_L001_R2_001.fastq.gz) is being adhered to and the first underscore delimited field is the sample_name.
    char r1_path[200];
    char r2_path[200];
    char *first_name = NULL;
    char *second_name = NULL;
    char *read1_num = NULL;
    char *read2_num = NULL;
    // extract sample name of first fastq file to check against
    char check[200];
    char *check_name;
    strcpy(check, basename(paths1[0]));
    check_name = strtok(check, "_");

    for (int t = 0; t < read1_count; t++)
    {
        strcpy(r1_path, basename(paths1[t]));
        strcpy(r2_path, basename(paths2[t]));

        first_name = strtok(r1_path, "_");
        for (int q = 0; q < 3; q++)
        {
            read1_num = strtok(NULL, "_");
        }
        second_name = strtok(r2_path, "_");
        for (int q = 0; q < 3; q++)
        {
            read2_num = strtok(NULL, "_");
        }

        // compare read pair sample names
        if (strcmp(first_name, check_name) != 0 || strcmp(second_name, check_name) != 0)
        {
            printf("Input fastqs must have the same sample name. Exiting...\n");
            exit(4);
        }
        // ensure each read1 filename contains "_R1_" and each read2 filename contains "_R2_", requres standard Illumina naming convention
        if (strcmp(read1_num, "R1") != 0)
        {
            printf("Read1 fastq file %s doesn not contain R1 label or is not in standard Illumina naming format. Exiting...\n", paths1[t]);
            exit(5);
        }
        if (strcmp(read2_num, "R2") != 0)
        {
            printf("Read2 fastq file %s doesn not contain R2 label or is not in standard Illumina naming format. Exiting...\n", paths2[t]);
            exit(6);
        }

    }

    // prepare log file
    // declare string to store formatted time
    char f_time[100];
    char log_file[500];
    snprintf(log_file, 500, "%s%s_BarCounter.log", outdir, first_name);
    char *user = getenv("USER");
    bool dir_exists;

    printf("\nBarCounter is being run by %s with the following arguments:\n", user);
    printf("\t-w %s (whitelist)\n", whitelist);
    printf("\t-t %s (taglist)\n", taglist);
    printf("\t-1 (read1 fastq)\n");
    for (int a = 0; a < read1_count; a++)
    {
        printf("\t\t%s\n", paths1[a]);
    }
    printf("\n\t-2 (read2 fastq)\n");
    for (int b = 0; b < read2_count; b++)
    {
        printf("\t\t%s\n", paths2[b]);
    }
    printf("\n\t-o %s (output directory)\n", outdir);
    printf("\n");

    // check fastq paths to ensure that each fastq file exists
    struct stat st;
    int is_file = -1;
    // check read 1 fastqs
    for (int f = 0; f < read1_count; f++)
    {
        is_file = stat(paths1[f], &st);
        // Check for file existence
        if (is_file != 0)
        {
            printf("Read 1 fastq path %s is invalid! Exiting...\n", paths1[f]);
            exit(7);
        }
    }
    // check read 2 fastqs
    for (int g = 0; g < read2_count; g++)
    {
        is_file = stat(paths2[g], &st);
        // Check for file existence
        if (is_file != 0)
        {
            printf("Read 2 fastq path %s is invalid! Exiting...\n", paths2[g]);
            exit(7);
        }
    }

    // Check if outdir exists and path is valid. If outdir doesn't exist create it
    if (stat(outdir, &st) == -1)
    {
        printf("Output directory %s doesn't exist. Creating %s\n", outdir,outdir);
        dir_exists = false;
        mkdir(outdir, 0777);
    } else {
        printf("Output will be written to existing directory %s\n", outdir);
        dir_exists = true;
        }

    // open log file for writing
    FILE *p_logfile = fopen(log_file, "w");

    fprintf(p_logfile, "%s\tBarCounter is being run by %s\n", get_datetime(f_time), user);
    fprintf(p_logfile, "%s\t-w %s (whitelist)\n", get_datetime(f_time), whitelist);
    fprintf(p_logfile, "%s\t-t %s (taglist)\n", get_datetime(f_time), taglist);
    fprintf(p_logfile, "%s\t-1 (read1 fastq)\n", get_datetime(f_time));
    for (int c = 0; c < read1_count; c++)
    {
        fprintf(p_logfile, "\t\t\t\t%s\n", paths1[c]);
    }
    fprintf(p_logfile, "%s\t-2 (read2 fastq)\n", get_datetime(f_time));
    for (int d = 0; d < read2_count; d++)
    {
        fprintf(p_logfile, "\t\t\t\t%s\n", paths2[d]);
    }
    fprintf(p_logfile, "%s\t-o %s (output directory)\n", get_datetime(f_time), outdir);
    if (dir_exists == false){
            fprintf(p_logfile, "%s\tOutput directory %s doesn't exist. Creating %s\n", get_datetime(f_time), outdir,outdir);
        } else {
            fprintf(p_logfile, "%s\tOutput will be written to existing directory %s\n", get_datetime(f_time), outdir);
        }

    // format output CSV file
    char counts_file[500];
    snprintf(counts_file, 500, "%s%s_Tag_Counts.csv", outdir, first_name);

    printf("Log file will be %s\n", log_file);
    printf("ADT counts will be written to %s\n\n", counts_file);
    fprintf(p_logfile, "%s\tADT counts will be written to %s\n", get_datetime(f_time), counts_file);

    // check if whitelist file is gzipped or plaintext format
    char *ext = NULL;
    ext = strrchr(whitelist, '.');
    bool gzipped;
    if (strcmp(ext, ".gz") == 0)
    {
        gzipped = true;
    }
    else if (strcmp(ext, ".txt") == 0)
    {
        gzipped = false;
    }
    else
    {
        printf("Unknown whitelist file extension %s. Exiting...\n", ext);
        fprintf(p_logfile, "%s\tUnknown whitelist file extension %s\n", get_datetime(f_time), ext);
        exit(8);
    }

    // load CSV taglist tags into array "tags" and names into array "names"
    int t_count = load_taglist(taglist, tags, names);
    // ensure taglist is not empty
    if (t_count == 0)
    {
        printf("Taglist is empty. Exiting...\n");
        fprintf(p_logfile, "%s\tTaglist is empty. Exiting...\n", get_datetime(f_time));
        exit(15);
    }
    // ensure adequate hamming distance between tags
    check_tag_dist(tags, t_count);

    // declare root nodes for each trie
    bc_node* bc_root = NULL;
    tag_node* tag_root = NULL;
    umi_node* umi_root = NULL;

    // create root node for each trie
    bc_root = calloc(1, sizeof(bc_node));
    tag_root = calloc(1, sizeof(tag_node));
    umi_root = calloc(1, sizeof(umi_node));

    // load taglist into tag trie
    if (!load_tag_trie(tags, tag_root, t_count))
    {
        printf("Failed to load all tags for processing. Exiting...\n");
        fprintf(p_logfile, "%s\tFailed to load all tags for processing. Exiting...\n", get_datetime(f_time));
        exit(18);
    }

    if (gzipped == true)
    {
        // load gzipped whitelist barcodes into barcode trie
        if (!load_bc_trie_gzipped(whitelist, bc_root, t_count))
        {
            printf("Failed to load barcodes for processing. Exiting...\n");
            fprintf(p_logfile, "%s\tFailed to load barcodes for processing. Exiting...\n", get_datetime(f_time));
            exit(21);
        }
    } else {
        // load plaintext whitelist barcodes into barcode trie
        if (!load_bc_trie(whitelist, bc_root, t_count))
        {
            printf("Failed to load barcodes for processing. Exiting...\n");
            fprintf(p_logfile, "%s\tFailed to load barcodes for processing. Exiting...\n", get_datetime(f_time));
            exit(21);
        }
    }

    // initilize variables for fastq parsing
    char curr_bc[BC_LEN + 1];
    char curr_umi[UMI_LEN + 1];
    char curr_tag[TAG_LEN + 1];
    int tag_index = -1;
    bc_node* p_bc = NULL;

    char R1_ID[200];
    char R1_seq[200];
    char R1_spacer[10];
    char R1_quals[200];

    char R2_ID[200];
    char R2_seq[200];
    char R2_spacer[10];
    char R2_quals[200];

    // set variable for mismatch testing: "tstart" is first base to test (index of bases), "test_n" is number of bases to test (4 only if base is N).
    char bases[4] = "ACGT";
    char temp[BC_LEN + 1];
    int tstart;
    int test_n;

    printf("\nBeginning fastq processing\n");
    fprintf(p_logfile, "%s\tBeginning fastq processing\n", get_datetime(f_time));

    // loop through each fastq read pair and process
    for (int x = 0; x < read1_count; x++)
    {

        // open input gzipped fastq files for reading
        gzFile pinR1 = gzopen(paths1[x], "r");
        gzFile pinR2 = gzopen(paths2[x], "r");

        // ensure all file pointers are valid
        if (pinR1 == NULL)
        {
            printf("Cannot open read1 fastq file %s\n", paths1[x]);
            fprintf(p_logfile, "%s\tCannot open read1 fastq file %s\n", get_datetime(f_time), paths1[x]);
            exit(22);
        }
        if (pinR2 == NULL)
        {
            printf("Cannot open read2 fastq file %s\n", paths2[x]);
            fprintf(p_logfile, "%s\tCannot open read2 fastq file %s\n", get_datetime(f_time), paths2[x]);
            exit(23);
        }
        printf("\nOpened input fastq files:\n%s\n%s\n\n",paths1[x],paths2[x]);
        fprintf(p_logfile, "%s\tOpened read1 fastq file %s\n", get_datetime(f_time), paths1[x]);
        fprintf(p_logfile, "%s\tOpened read2 fastq file %s\n", get_datetime(f_time), paths2[x]);


        // read input fastq files read by read and check if value is in trie
        while (true)
        {
            // assign each line of reach for each fastq file to variables
            // if any fastq files return NULL from gzgets the EOF has been reached, break the loop
            if (!gzgets(pinR1, R1_ID, 200))
            {
                break;
            }
            gzgets(pinR1, R1_seq, 200);
            gzgets(pinR1, R1_spacer, 10);
            gzgets(pinR1, R1_quals, 200);
            if (!gzgets(pinR2, R2_ID, 200))
            {
                break;
            }
            gzgets(pinR2, R2_seq, 200);
            gzgets(pinR2, R2_spacer, 10);
            gzgets(pinR2, R2_quals, 200);

            // update read count
            total_reads++;

            // parse read1 seq to assign cell barcode and UMI
            for (int b = 0; b < BC_LEN; b++)
            {
                curr_bc[b] = R1_seq[BC_FIRST + b];
                curr_bc[BC_LEN] = '\0';
            }
            for (int u = 0; u < UMI_LEN; u++)
            {
                curr_umi[u] = R1_seq[UMI_FIRST + u];
                curr_umi[UMI_LEN] = '\0';
            }
            // parse read2 seq to assign tag
            for (int t = 0; t < TAG_LEN; t++)
            {
                curr_tag[t] = R2_seq[TAG_FIRST + t];
                curr_tag[TAG_LEN] = '\0';
            }

            // ensure barcode is valid and in whitelist
            p_bc = get_bc_leaf(curr_bc, bc_root, BC_LEN);

            // allow for single mismatch at low quality basecall in barcode. Track correct barcode in bc_node pointer. */
            if (p_bc == NULL)
            {
                // copy barcode seq to temp and check mismatches at low quality bases
                strncpy(temp, curr_bc, BC_LEN + 1);
                test_n = 3;

                for (int m = 0; m < BC_LEN; m++)
                {
                    // check if basecall is below quality score of 20
                    if (R1_quals[m] < LOW_Q)
                    {
                        // check barcode string with each of the three remaining bases substituted at postion m
                        switch(temp[m])
                        {
                            case 'A': tstart = 1; break;
                            case 'C': tstart = 2; break;
                            case 'G': tstart = 3; break;
                            case 'T': tstart = 4; break;
                            // test all four bases in case of N
                            case 'N': test_n = 4; tstart = 0; break;
                        }
                        for (int r = 0; r < test_n; r++)
                        {
                            temp[m] = bases[(tstart+r)%4];

                            // if temp barcode with one mismatch corresponds to a sample, update out_s and break
                            p_bc = get_bc_leaf(temp, bc_root, BC_LEN);
                            if (p_bc != NULL)
                            {
                                break;
                            }
                        }
                    }
                    temp[m] = curr_bc[m];
                    // if match has been found break main loop
                    if (p_bc != NULL)
                    {
                        corrected_barcodes++;
                        break;
                    }
                }
            }
            if (p_bc != NULL)
            {
                // update valid barcode count
                valid_barcodes++;

                // ensure read2 seq is in the taglist
                tag_index = get_tag_index(curr_tag, tag_root);
                if (tag_index != -1)
                {
                    valid_tags++;

                    // if UMI added for barcode: update cell barcode tag counts
                    if (add_umi(curr_umi, umi_root, t_count, tag_index, curr_bc))
                    {
                        // update cell barcode tag count
                        p_bc->counts[tag_index]++;
                        p_bc->total++;
                    }
                }
            }
        }
        gzclose(pinR1);
        gzclose(pinR2);
    }

    // write tag counts to output CSV file
    FILE *out_counts = fopen(counts_file, "w");
    // test printing counts by re-reading whitelist file
    gzFile p_white = gzopen(whitelist, "r");
    int i;
    char barcode[20];
    p_bc = NULL;

    // write header
    fprintf(out_counts, "cell_barcode,total");
    for (int n = 0; n < t_count; n++)
    {
        fprintf(out_counts,",%s", names[n]);
    }
    fprintf(out_counts, "\n");

    while (true)
    {
        // read line of whitelist into 'barcode'
        gzgets(p_white, barcode, 20);
        barcode[16] = '\0';

        // break when the end of the whitelist is reached
        if (gzeof(p_white))
        {
            break;
        }

        p_bc = get_bc_leaf(barcode, bc_root, BC_LEN);
        if (p_bc != NULL)
        {
            // only write cell barcodes with counts
            if (p_bc->total != 0)
            {
                fprintf(out_counts, "%s,%li", barcode, p_bc->total);
                for (int fg = 0; fg < t_count; fg++)
                {
                    fprintf(out_counts, ",%i", p_bc->counts[fg]);
                }
                fprintf(out_counts, "\n");
            }
        }
    }
    gzclose(p_white);
    fclose(out_counts);

    // unload UMI trie
    if (!unload_umi_trie(umi_root, t_count))
    {
        printf("UMIs failed to unload\n");
        fprintf(p_logfile, "%s\tUMIs failed to unload\n", get_datetime(f_time));
    }
    // unload barcode trie
    if (!unload_bc_trie(bc_root))
    {
        printf("Barcodes failed to unload\n");
        fprintf(p_logfile, "%s\tBarcodes failed to unload\n", get_datetime(f_time));
    }
    // unload tag trie
    if (!unload_tag_trie(tag_root))
    {
        printf("Tags failed to unload\n");
        fprintf(p_logfile, "%s\tTags failed to unload\n", get_datetime(f_time));
    }

    printf("Processing complete\n");
    printf("Total reads processed: %lli\n", total_reads);
    printf("Uncorrected barcodes: %lli\n", valid_barcodes - corrected_barcodes);
    printf("Corrected barcodes: %lli\n", corrected_barcodes);
    printf("Total Valid barcodes: %lli\n", valid_barcodes);
    printf("Valid tags: %lli\n", valid_tags);
    printf("\nFINISHED\n");

    fprintf(p_logfile, "%s\tProcessing complete\n", get_datetime(f_time));
    fprintf(p_logfile, "%s\tTotal reads processed: %lli\n", get_datetime(f_time), total_reads);
    fprintf(p_logfile, "%s\tUncorrected barcodes: %lli\n", get_datetime(f_time), valid_barcodes - corrected_barcodes);
    fprintf(p_logfile, "%s\tCorrected barcodes: %lli\n", get_datetime(f_time), corrected_barcodes);
    fprintf(p_logfile, "%s\tTotal Valid barcodes: %lli\n", get_datetime(f_time), valid_barcodes);
    fprintf(p_logfile, "%s\tValid tags: %lli\n", get_datetime(f_time), valid_tags);
    fprintf(p_logfile, "%s\tFINISHED\n", get_datetime(f_time));

    fclose(p_logfile);

    free(paths1);
    free(paths2);

    return 0;
}

// return string f_time with formatted current GMT (UTC)
char* get_datetime(char* f_time)
{
    // declare raw time variable
    time_t current_time;
    // initialize with current time
    time(&current_time);
    // initialize time pointer to GMT
    struct tm *gmt = gmtime(&current_time);
    // format time string
    strftime(f_time, 100, "%d-%b-%Y %X", gmt);

    return f_time;
}