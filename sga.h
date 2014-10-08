#include <stdio.h>
#include <stdlib.h>
#include "stdint.h"
#include <math.h>
#include <time.h>

// Global variables
int generations;     // Number of generations to run.
int num_chromo;      // Size of Population; Number of chromosomes.
int num_genes;       // Number of 16-bit "genes" per chromosome.
double rate_elite;   // Rate of elite chromosomes per population.
double rate_cross;   // Rate of crossover per population.
double rate_mutate;  // Rate of mutation per population.
int num_elite;       // Number of elite chromosomes given rate_elite.
int num_crossed;     // Number of chromosomes crossed given rate_cross.
int num_remaining;   // Number of chromosomes that are drawn directly from old population.
int num_mutations;   // Number of bits mutated given rate_mutate;

char* file_name;
FILE *output;

// Data Structures
typedef struct {
	uint16_t gene[10];
	double fitness;
} chromosome;

// Printing Methods
void print_stats();
void print_help();

// Helper Methods
double gene2value(uint16_t gene);
uint16_t value2gene(double value);
int choose_chromosome(chromosome* pop);

// SGA Methods
void setup(int argc, char* argv[], int rank);
void setup_generated();
void setup_defaults();

void initialize(chromosome* pop);
double mpi_computation(chromosome* pop, int size);
void mpi_calculation(chromosome* pop, int generation, int rank, int size);
void sort(chromosome* pop);
void elite(chromosome* pop, chromosome* new_pop);
void cross_genes(uint16_t* gene1, uint16_t* gene2, int index);
void cross_chromosomes(chromosome* chromo1, chromosome* chromo2);
void crossover(chromosome* pop, chromosome* new_pop);
void move_remaining(chromosome* pop, chromosome* new_pop);
void mutate(chromosome* pop);
