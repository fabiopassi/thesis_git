#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include "lib.h"

/* Functions body */

void eval_res_rel(int** system, int num_samples, int N_spins, int n, int mapp, float** res, float** rel){

    /* Pick n survived spins out of N_spins, with no repetition */
    int survived_idx[n];
    for (int i = 0; i < n; i++) {survived_idx[i] = -1;}

    for (int i = 0; i < n; i++) {
        int stop = 0;
        while(! stop) {
            int idx = rand() % N_spins;
            int add = 1;
            for(int j = 0; j < n; j++) {
                if (survived_idx[j] == idx ) {add = 0; break;}
            }
            if (add == 1) {survived_idx[i] = idx; stop = 1;}
        }
    }

    /* Build a subsystem with just the survived spins */
    int** system_survived = malloc(num_samples * sizeof(*system_survived));
	for(int i = 0; i < num_samples; i++) {
		system_survived[i] = malloc(n * sizeof(*system_survived[i]));
	}

    for(int i = 0; i < num_samples; i++) {
        for(int j = 0; j < n; j++) {
            system_survived[i][j] = system[i][survived_idx[j]];
        }
    }

    /* Assign labels */
    int labels[num_samples];
    for (int i = 0; i < num_samples; i++) {labels[i] = -1;}

    for(int i = 0; i < num_samples; i++) {

        /* If the i-th state is equal to a previous state, I assign the same label */
        for(int k = 0; k < i; k++) {
            int same = 1;
            for(int j = 0; j < n; j++) {
                if (system_survived[i][j] != system_survived[k][j]) {
                    same = 0;
                    break;
                }
            }
            if (same == 1) {
                labels[i] = labels[k];
            }
        }

        /* If the i-th state is different from every previous state, I use i as label */
        if (labels[i] == -1) {
            labels[i] = i;
        }

    }

    /* Calculate resolution */
    int unique_labels[num_samples];
    int count_labels[num_samples];
    int len = 0;
    for (int i = 0; i < num_samples; i++) {unique_labels[i] = -1; count_labels[i] = 0;}

    for(int i = 0; i < num_samples; i++) {
        int add = 1;
        for(int j = 0; j < len; j++) {
            if (labels[i] == unique_labels[j]){
                add = 0;
                count_labels[j]++;
                break;
            }
        }
        if (add == 1) {
            unique_labels[len] = labels[i];
            count_labels[len]++;
            ++len;
        }
    }

    float resolution = 0;
    for(int i=0; i<len; i++) {
        if (count_labels[i] != 0) {
            float prob = (float) count_labels[i] / num_samples;
            resolution -= prob * log(prob);
        }
    }
    
    res[n][mapp] = resolution;

    /* Calculate relevance */
    int unique_counts[len];
    int count_counts[len];
    int num_diff_counts = 0;
    for (int i = 0; i < len; i++) {count_counts[i] = 0; unique_counts[i] = -1;}

    for(int i = 0; i < len; i++) {
        int add = 1;
        for(int j = 0; j < num_diff_counts; j++) {
            if (count_labels[i] == unique_counts[j]){
                add = 0;
                count_counts[j]++;
                break;
            }
        }
        if (add == 1) {
            unique_counts[num_diff_counts] = count_labels[i];
            count_counts[num_diff_counts]++;
            ++num_diff_counts;
        }
    }

    float relevance = 0;
    for(int i=0; i<num_diff_counts; i++) {
        if (count_counts[i] != 0) {
            float prob = (float) count_counts[i] * unique_counts[i] / num_samples;
            relevance -= prob * log(prob);
        }
    }
    
    rel[n][mapp] = relevance;

    /* Free stuff */
    for(int i = 0; i < num_samples; i++) {
		free(system_survived[i]);
	}
	free(system_survived);

    return;
}
