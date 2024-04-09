#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <time.h>
#include "lib.h"

/* What the program does (docstring) */

int main (int argc, char **argv) {

	/* Variables */
	int N_spins = 10;
	int num_samples = 200;

	int** system = malloc(num_samples * sizeof(*system));
	for(int i = 0; i < num_samples; i++) {
		system[i] = malloc(N_spins * sizeof(*system[i]));
	}

	int n_mappings = 200;

	float** res = malloc( (N_spins + 1) * sizeof(*res));
	float** rel = malloc( (N_spins + 1) * sizeof(*rel));
	for(int i = 0; i < N_spins + 1; i++) {
		res[i] = malloc(n_mappings * sizeof(*res[i]));
		rel[i] = malloc(n_mappings * sizeof(*rel[i]));
	}

	/* Initialize the system */
	srand(time(NULL));

	int biased = 1;			/* 1 = biased; 0 = unbiased */

	if (biased) {
		/* Biased spins */
		for(int i = 0; i < num_samples; i++) {
			for (int j = 0; j < N_spins; j++){
				float  p;
				if (j < 5) {
					p = 1 - (float)j/10;
				} else {
					p = 0.5;
				}

				if ((float) rand()/RAND_MAX < p) {
					system[i][j] = 1;
				} else {
					system[i][j] = 0;
				}
			}
		}
	} else {
		/* Unbiased spins */
		for(int i = 0; i < num_samples; i++) {
			for (int j = 0; j < N_spins; j++){
				if ((float) rand()/RAND_MAX < 0.5) {
					system[i][j] = 1;
				} else {
					system[i][j] = 0;
				}
			}
		}
	}

	/* Status message */
	printf("\nParameters used for the simulation:\n");
	printf("\t-> N_spins = %d\n", N_spins);
	printf("\t-> biased = %d\n", biased);
	printf("\t-> samples = %d\n", num_samples);
	printf("\t-> mappings = %d\n", n_mappings);
	printf("\nStarting calculation ...\n");

	/* For each number of retained spins, I do n_mappings mappings and I evaluate the resolution and relevance */
	# pragma omp parallel
	{
		#pragma omp for
		for(int n = 0; n <= N_spins; n++) {
			for(int mapp = 0; mapp < n_mappings; mapp++) {
				eval_res_rel(system, num_samples, N_spins, n, mapp, res, rel);
			}
		}
	}
	#pragma omp barrier

	/* Status message */
	printf("\nCalculation finished.\n");

	/* Print on a file the results */
	int buffer_size = 256;
	char* output_file = malloc(buffer_size * sizeof(*output_file));
	snprintf(output_file, buffer_size, "./build/res_rel.txt");
	FILE* pf = fopen(output_file, "w");
	if (pf == NULL) {
		fprintf(stderr, "Error: could not open the file in the location\n%s\nPlease be sure that you are launching the program from the makefile directory or change the output file path.\n", output_file);
		exit(EXIT_FAILURE);
	}

	/* Status message */
	printf("\nWriting output to file: %s\n\n", output_file);

	for (int i = 0; i < n_mappings; i++) {
		for (int j = 0; j <= N_spins; j++) {
			fprintf(pf, "%f\t%f\t", res[j][i], rel[j][i]);
		}
		fprintf(pf, "\n");
	}

	fclose(pf);
	free(output_file);

	/* Status message */
	printf("Finished.\n\n");

	/* Free variables */
	for(int i = 0; i < num_samples; i++) {
		free(system[i]);
	}
	free(system);

	for(int i = 0; i < N_spins + 1; i++) {
		free(res[i]);
		free(rel[i]);
	}
	free(res);
	free(rel);

	return 0;
}

