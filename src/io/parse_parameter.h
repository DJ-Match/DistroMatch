/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef PARSE_PARAMETER_H
#define PARSE_PARAMETER_H

#include <argtable3.h>

#include <algorithm>
#include <iostream>
#include <regex.h>
#include <sstream>
#include <string.h>

#include <config.h>

#define ALGOS "MRepLM|kMM|kCS|MRepS"
#define PARTITIONING "metis|umetis|uwmetis|ukahip|random|diststream|none"

bool parse_algorithm(const char *algo, Config &config) {
    if (strcmp("MRepLM", algo) == 0) {
        config.algo = MREP_LM;
    } else if (strcmp("kCS", algo) == 0) {
        config.algo = K_CS;

    } else if (strcmp("MRepS", algo) == 0) {
        config.algo = MREP_S;

    } else if (strcmp("kMM", algo) == 0) {
        config.algo = K_MM;

    } else {
        return false;
    }

    return true;
}

bool parse_partitioning(const char *partitioning, Config &config) {
    if (strcmp("metis", partitioning) == 0) {
        config.partitioning = METIS;
    } else if (strcmp("umetis", partitioning) == 0) {
        config.partitioning = METIS_UNWEIGHTED;
    } else if (strcmp("uwmetis", partitioning) == 0) {
        config.partitioning = METIS_NODE_WEIGHTS_UNWEIGHTED;
    } else if (strcmp("ukahip", partitioning) == 0) {
        config.partitioning = KAHIP_UNWEIGHTED;
    } else if (strcmp("random", partitioning) == 0) {
        config.partitioning = RANDOM_PARTITIONING;
    } else if (strcmp("diststream", partitioning) == 0) {
        config.partitioning = DIST_STREAM;
    } else {
        config.partitioning = NO_PARTITIONING;
    }
    return true;
}

int parse_parameters(int argn, char **argv, std::string &infile,
                     std::string &outfile, Config &config) {
    const char *programm = argv[0];

    // Setup argtable parameters.
    // For further details, see:
    // - https://www.argtable.org/tutorial/
    // - https://linux.die.net/man/3/argtable
    struct arg_lit *help = arg_lit0(NULL, "help", "Print help.");
    struct arg_file *in_file =
        arg_file1("i", "input_file", "<file>", "Path to graph file.");
    struct arg_file *out_file =
        arg_file0("o", "output_file", "<file>", "Path to .sol output file.");
    struct arg_rex *algo =
        arg_rex1(NULL, "algo", "^(" ALGOS ")$", "VARIANT", REG_EXTENDED,
                 "Define the algorithm to solve kEC: [" ALGOS "]");
    struct arg_int *k =
        arg_int1(NULL, "k", NULL, "The number of matchings to compute.");
    struct arg_end *end = arg_end(100);
    struct arg_rex *partitioning = arg_rex1(
        "p", "partitioning", "^(" PARTITIONING ")$", "VARIANT", REG_EXTENDED,
        "Define the partitioning used: [" PARTITIONING "]");
    struct arg_int *seed =
        arg_int1(NULL, "seed", NULL,
                 "Select a seed for the random fuction. [Default: 0]");
    // struct arg_lit *dp = arg_lit0(
    //     NULL, "dp",
    //     "Compute the double amount of matchings and merge two optimally.");
    struct arg_dbl *ccm_epsilon =
        arg_dbl0(NULL, "ccm_epsilon", NULL,
                 "In CCM all edges of a vertex v with weight "
                 "(1-ccm_epsilon)*max_weight_v are considered candidates. "
                 "Range: [0,1]. Default: 0.1");

#ifdef TESTING
    struct arg_lit *par1 = arg_lit0(NULL, "par1", "Test parameter 1.");
    struct arg_lit *par2 = arg_lit0(NULL, "par2", "Test parameter 2.");
#endif

    void *argtable[] = {help, in_file,      out_file, algo,
                        k,    partitioning, seed,     ccm_epsilon,
#ifdef TESTING
                        par1, par2,
#endif
                        end};

    // Parse arguments.
    arg_parse(argn, argv, argtable);

    // Catch case that help was requested.
    if (help->count > 0) {
        printf("Usage: %s", programm);
        arg_print_syntax(stdout, argtable, "\n");
        arg_print_glossary(stdout, argtable, "  %-40s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        exit(1);
        return 1;
    }

    if (in_file->count == 1) {
        infile = std::string(in_file->filename[0]);
    } else {
        printf("Define the graph file: ");
        std::cin >> infile;
    }

    if (out_file->count == 1) {
        outfile = std::string(out_file->filename[0]);
    }

    if (algo->count > 0) {
        std::string algo_s = algo->sval[0];
        while (!parse_algorithm(algo_s.c_str(), config)) {
            printf("Define algorithm (%s): ", ALGOS);
            std::cin >> algo_s;
        }
    } else {
        std::string algo;
        do {
            printf("Define algorithm (%s): ", ALGOS);
            std::cin >> algo;
        } while (!parse_algorithm(algo.c_str(), config));
    }

    if (k->count > 0) {
        config.k = k->ival[0];
    } else {
        int k = 0;
        do {
            printf("Define the number of matchings k to compute: ");
            std::cin >> k;
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } while (k == 0);

        config.k = k;
    }

    if (partitioning->count > 0) {
        std::string partition = partitioning->sval[0];
        parse_partitioning(partition.c_str(), config);
    } else {
        config.partitioning = NO_PARTITIONING;
    }

    if (seed->count > 0) {
        config.seed = seed->ival[0];
    } else {
        config.seed = 0;
    }

    if (ccm_epsilon->count > 0) {
        config.ccm_epsilon = ccm_epsilon->dval[0];
    }

#ifdef TESTING
    if (par1->count > 0) {
        config.test_param_1 = true;
    } else {
        config.test_param_1 = false;
    }

    if (par2->count > 0) {
        config.test_param_2 = true;
    } else {
        config.test_param_2 = false;
    }

#endif

    return 0;
}

#endif