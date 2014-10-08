#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "stdint.h"
#include <math.h>
#include <time.h>
#include "sga.h"

void print_stats(char* header) {
    printf("\nSimple Genetic Algorithm\n");
	printf("----- %s -----\n", header);
    printf("Generations: %d\n", generations);
	printf("Chromosomes: %d\n", num_chromo);
	printf("Genes per Chromosome: %d\n", num_genes);
    printf("Elite rate: %+4.2f\n", rate_elite);
    printf("Crossover Rate: %+4.2f\n", rate_cross);
    printf("Mutation rate: %+4.2f\n", rate_mutate);
    printf("Per Generation:\n");
	printf("  elites: %d\n", num_elite);
    printf("  crossed chromosomes: %d\n", num_crossed);
	printf("  remaining chromosomes: %d\n", num_remaining);
    printf("  mutations: %d genes\n", num_mutations);
	printf("\n");
    printf("Output printed to: %s\n", file_name);
    printf("\n"); 
}

void print_help() {
    printf("\nSimple Genetic Algorithm\n");
    printf("------------------------\n");
    printf("Usage: mpirun [OPTIONS] ./sga [OPTIONS]\n");
    printf("  See man mpirun for MPI usage.\n");
    printf("Options:\n");
    printf("  -f {file name}    output file name\n");
    printf("  -g {value}        number of generations   0 < value\n");
    printf("  -p {value}        size of population      0 < value\n");
    printf("  -c {value}        crossover rate:         0 < value < 1\n");
    printf("  -m {value}        mutation rate:          0 < value < 1\n");
    printf("  -e {value}        elite rate:             0 < value < 1\n");
    printf("\n");
    printf("  -d                display default values\n");
    printf("  -h                display help\n\n");
    
}

/**
 * Transforms a 16 bit integer to a value in the range of
 * -1024 <= igene <= 1023.  
 * Notes:
 * 0000 0000 0000 0000 : 16-bits
 *       000 0000 0000 : 12-bits used 2^11 = 2048
 * 0000 0              : 5-bits unused
 * 
 * (16-bits / 5 bits) - 1024 (signed shift)
 */
double gene2value(uint16_t gene) {
	return (double)((gene / 32.0) - 1024);
}

/**
 * Mutates bits in chromosomes.  The range of chromosomes for 
 * possible mutations excludes the elite chromosomes.  The number of 
 * mutations is determine by @global num_mutation.
 */
void mutate(chromosome* pop) {
	int ichromo;
	int igene;
	int ibit;
	
	int count;
	for (count = 0; count < num_mutations; count++) {
		// Get random chromosome, excluding elites
		ichromo = (rand() % (num_chromo - num_elite)) + num_elite;
		// Get random gene
		igene = rand() % num_genes;
		// Get random bit
		ibit = rand() % 16;
		// Mutate gene
		pop[ichromo].gene[igene] ^= 1 << ibit;
	}
}

/**
 * Fills in the remainder of the new population with chromosomes
 * that are non-elite, and had the lowest fitness in the pervious
 * population.
 */
void move_remaining(chromosome* pop, chromosome* new_pop) {
	int start = num_elite + num_crossed - 1;
	
    // For every remaining chromosome, move it to the new population.
	int index;
	for (index = 0; index <= num_remaining; index++) {
		new_pop[start + index] = pop[start + index];
	}
}

/**
 * A roulette-wheel algorithm to choose a chromosome for crossing.
 * Uses num_elite to prevent crossing with elite chromosomes.
 */
int choose_chromosome(chromosome* pop) {
	double sum_fitness = 0;
	
	// Compute the sum of the fitness values for the population, 
    // excluding the elites.
	int ichromo;
	for (ichromo = num_elite; ichromo < num_chromo; ichromo++) {
		sum_fitness += pop[ichromo].fitness;
	}
	
	// Get a random number from 0 to sum of the fitness values.
	int rand_num = rand() % (int)sum_fitness;
	
	// Slough off the size of each fitness from the sum of the fitness values
	// until we have reached 0.  Then we have the index for the proper chromosome.
	//   Note: excludes the elites
	for (ichromo = num_elite; ichromo < num_chromo && rand_num > 0; ichromo++) {
		rand_num = rand_num - pop[ichromo].fitness;
	}
	
	// Remove the last increment before returning the chromosome index.
	return ichromo - 1;  
}

/**
 * Swap gene at a given pivot bit.
 * if index = 11:
 * 	        GENE1                 GENE2 
 *   AAAA AAAA AAAA AAAA | BBBB BBBB BBBB BBBB =>
 *   BBBB BAAA AAAA AAAA | AAAA ABBB BBBB BBBB
 * Note: upper bits are swapped, so the swapping of the chromosome
 * 		 should take place on the lower indexes. 
 */
void cross_genes(uint16_t* gene1, uint16_t* gene2, int pivot) {
	uint16_t new_gene1 = 0;
	uint16_t new_gene2 = 0;
	
	int i;
    // Transfer the bottom of the genes up to the pivot.
	for (i = 0; i < pivot; i++) {
		if (*gene1 & (1 << i) ) {
			new_gene1 |= (1 << i);
		}
		if (*gene2 & (1 << i) ) {
			new_gene2 |= (1 << i);
		}
	}
	
    // Swap the top of the gene, starting at the pivot.
	for (i = pivot; i < 16; i++) {
		if (*gene1 & (1 << i) ) {
			new_gene2 |= (1 << i);
		}
		if (*gene2 & (1 << i) ) {
			new_gene1 |= (1 << i);
		}
	}
	
    // Assign the new gene values to the old gene values.
	*gene1 = new_gene1;
	*gene2 = new_gene2;
}

/**
 * Cross two chromosomes by determining a pivot point, keeping the
 * upper half the same, and swapping the lower half.
 */
void cross_chromosomes(chromosome* chromo1, chromosome* chromo2) {
	// Two temporary chromosomes
    chromosome new_chromo1;
	chromosome new_chromo2;

    // Determine the pivot locations.
	int cut = rand() % 160; 		// Determine bit on chromosome to cut
	int gene_to_cross = (cut / 16); // get the gene for that number
	int bit_to_cross = (cut % 16);	// get the bit on the gene

    // Pull out the genes that contain the cut and cross them.
	uint16_t cross_gene1 = chromo1->gene[gene_to_cross];
	uint16_t cross_gene2 = chromo2->gene[gene_to_cross];
	cross_genes(&cross_gene1, &cross_gene2, bit_to_cross);

	// Swap the genes up to the crossover gene.
	int igene;
	for (igene = 0; igene < gene_to_cross; igene++) {
		new_chromo1.gene[igene] = chromo2->gene[igene];
		new_chromo2.gene[igene] = chromo1->gene[igene];
	}

	// Place the crossover genes.
	new_chromo1.gene[igene] = cross_gene1;
	new_chromo2.gene[igene] = cross_gene2;

	// Transfer the rest of the genes to their chromosomes.
	for (igene++; igene < num_genes; igene++) {
		new_chromo1.gene[igene] = chromo1->gene[igene];
		new_chromo2.gene[igene] = chromo2->gene[igene];
	}

	// Assign the new chromosomes to their old position in the population.
	*chromo1 = new_chromo1;
	*chromo2 = new_chromo2;
}


/**
 * Perform the crossover operation on pairs of chromosomes in the
 * population and transfer the new chromosome pairs to the new population.
 */
void crossover(chromosome* pop, chromosome* new_pop) {
    // We perform crossover operations on half the amount of chromosomes
    // that need to be crossed.
	int num_crossover_steps = num_crossed / 2;

	// For every pair of chromosomes to be crossed...	
	int icrossed;
	for (icrossed = 0; icrossed < num_crossover_steps; icrossed++) {
		// get an index for the chromosomes to be crossed and...
		int chromo1 = 0;
		int chromo2 = 0;
		while( chromo1 == chromo2 ) {
			chromo1 = choose_chromosome(pop);
			chromo2 = choose_chromosome(pop);
		}
        
        // get the indexes for the new chromosomes on the new population,
		// skipping the elites...
		int new_chromo1 = num_elite + (icrossed * 2);
		int new_chromo2 = num_elite + (icrossed * 2) + 1;
		
		// move those chromosomes to the new population,
		// putting them in the left-most position available...
		new_pop[new_chromo1] = pop[chromo1];
		new_pop[new_chromo2] = pop[chromo2];
		
		// cross over the two chromosomes in their new locations.
		cross_chromosomes( &(new_pop[new_chromo1]), 
						   &(new_pop[new_chromo2]) );
	}
}

/**
 * Determine the elite chromosomes in the current population
 * and move them to the new population. Assumes population has
 * been sorted.
 * @require num_elite > 0
 */
void elite(chromosome* pop, chromosome* new_pop) {
    // Move the elite chromosomes from the current population
    // to the new population.
	int ichromo;
    for (ichromo = 0; ichromo < num_elite; ichromo++) {
		new_pop[ichromo] = pop[ichromo];
	}
}

/**
 * Sort the chromosomes based on their fitness.
 */
void sort(chromosome* pop) {
	// For every chromosome in the population...
    int ichromo;
	for (ichromo = 0; ichromo < num_chromo; ichromo++) {
        // store the chromosome...
		chromosome temp = pop[ichromo];
        // identify the location of a possible hole (if we move it)...
		int hole = ichromo;
        // as long as the fitness of the chromosome at the current
        // location is greater than the one at the possible hole...
		while (hole > 0 && temp.fitness > pop[hole-1].fitness) {
            // put the current one in the hole...
			pop[hole] = pop[hole-1];
            // and continue to shift...
			hole--;
		}
        // Now put the chromosome in the hole we created.
		pop[hole] = temp;
	}
}

/**
 * Compute the fitness of chromosomes in the process-specific
 * population.
 */
double mpi_computation(chromosome* local_pop, int local_size) {
	double best_fitness;  // Stores the best fitness value for this population.
    
    // For every chromosome in the local population...
    int ichromo;
	for (ichromo = 0; ichromo < local_size; ichromo++) {
		double sigma = 0;   // The value of the sigma portion of the equation.
        // For every gene...
		int igene;
        for (igene = 0; igene < num_genes; igene++) {
            // aggregate the sigma values.
			sigma += pow(gene2value(local_pop[ichromo].gene[igene]), 2);
		}
		sigma = sigma / 4000;

		double pi = 0;      // The value of the pi portion of the equation.
        // For every gene...
		for (igene = 0; igene < num_genes; igene++) {
            // aggregate the pi values.
			pi *= cos( gene2value(local_pop[ichromo].gene[igene]) / sqrt( (double)(igene+1) ) );
		}
		pi *= pi;
		
        // and finally, assign the fitness to the chromosome.
		local_pop[ichromo].fitness = 1 + sigma - pi;
        
        // ##### Below is for preserving the best fitness value for this population. #####
        // Initialize to first fitness 
        if (ichromo == 0) {
            best_fitness = local_pop[ichromo].fitness;
        }
        
        // Store the best fitness value for this process
        if (local_pop[ichromo].fitness > best_fitness) {
            best_fitness = local_pop[ichromo].fitness;
        }
	}
    return best_fitness;
}


/**
 * Uses MPI to calculate the fitness for the population.  First, the population
 * is divided amongst the processes.  Then each process calculates the 
 * fitness value for the chromosomes it was assigned.  The master process
 * then gathers the results.  
 * 
 * When the population cannot be evenly distributed
 * amongst the processes, it further divides the leftover chromosomes and 
 * distributes each to one of the processes.  No process will have to operate
 * on more than one extra chromosome.
 * 
 * This implimentation records additional statistics about each process 
 * for the first 10 generations.
 */
void mpi_calculation(chromosome* pop, int generation, int rank, int size){
    // Process Variables
	chromosome* local_pop;      // An array of chromosomes specific to this process.
	chromosome* extra_chromo;   // An additional chromosome for this process.
	int left_overs = num_chromo % size; // The amount of chromosomes that will not
                                        // be evenly distributed between processes.
    
    // Allocate memory for the process-specific chromosomes.
    local_pop = (chromosome*)malloc(sizeof(chromosome) * (num_chromo/size));
    extra_chromo = (chromosome*)malloc(sizeof(chromosome));
	
    // Setup MPI types
	MPI_Datatype gene_array_mpi_t;  // An MPI gene
	MPI_Type_contiguous(10, MPI_UNSIGNED_SHORT, &gene_array_mpi_t);
	MPI_Type_commit(&gene_array_mpi_t);
	
	MPI_Aint gene_array_extent;
	MPI_Type_extent(gene_array_mpi_t, &gene_array_extent);
	
	MPI_Datatype types[] = {gene_array_mpi_t, MPI_DOUBLE};
	MPI_Aint offsets[] = {0, gene_array_extent};
	int blocks[] = {1, 1};
	
	MPI_Datatype chromosome_mpi_t;	// An MPI chromosome
	MPI_Type_struct(2, blocks, offsets, types, &chromosome_mpi_t);
	MPI_Type_commit(&chromosome_mpi_t);
	
	MPI_Datatype chromosome_array_mpi_t; // An MPI chromosome array
	MPI_Type_contiguous(num_chromo/size, chromosome_mpi_t, &chromosome_array_mpi_t);
	MPI_Type_commit(&chromosome_array_mpi_t);
    
    // ##### Scatter #####
	MPI_Barrier(MPI_COMM_WORLD);
    // Have the master process Scatter() the population to multiple processes.
    int err = MPI_Scatter(pop, num_chromo/size, chromosome_mpi_t,
						  local_pop, 1, chromosome_array_mpi_t,
						  0, MPI_COMM_WORLD);
                          
    // If Scatter() is uneven, scatter the leftovers                      
    if (left_overs != 0) {
		int err = MPI_Scatter(&(pop[(num_chromo/size)*size]), 1, chromosome_mpi_t,
							  extra_chromo, 1, chromosome_mpi_t,
							  0, MPI_COMM_WORLD);
	}
    
    // ##### Computation #####
	MPI_Barrier(MPI_COMM_WORLD);
    double best_fitness ;  // Holds the best fitness calculated by this process.
    
    // Compute the fitness for the process's chromosomes.
    best_fitness = mpi_computation(local_pop, num_chromo/size);
    // If were chromosomes leftover, compute the one passed to this process
    // if one was given to it.
    if (left_overs != 0 && rank < left_overs) {
        double left_over_fitness;
		left_over_fitness = mpi_computation(extra_chromo, 1);
        // consider the leftover chromosome fpr the process's best fitness.
        if (left_over_fitness > best_fitness) {
            best_fitness = left_over_fitness;
        }
	}
    
    //##### Gather #####
	// Have the master process Gather() the results from the processes.
	MPI_Barrier(MPI_COMM_WORLD);
	int err2 = MPI_Gather(local_pop, 1, chromosome_array_mpi_t,
						  pop, num_chromo/size, chromosome_mpi_t,
						  0, MPI_COMM_WORLD);
    // If there were leftover chromosomes, have the master process gather
    // those as well.
    if (left_overs != 0) {
		int err2 = MPI_Gather(extra_chromo, 1, chromosome_mpi_t,
						      &(pop[(num_chromo/size)*size]), 1, chromosome_mpi_t,
						      0, MPI_COMM_WORLD);
	}
    
    //##### Record Process Data #####
    MPI_Barrier(MPI_COMM_WORLD);
    // We are only recording process-specific data for the first 10 generations.
    if (generation < 10) {
        double* fitness_array;  // Best Fitness value produced by each process.
        fitness_array = (double*)malloc(sizeof(double) * size);
        
        // Gather the value stored in each process that held the highest
        // fitness value it computed.
        MPI_Gather(&best_fitness, 1, MPI_DOUBLE, 
                   fitness_array, 1, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        // The master process records the data to a file.
        if (rank == 0) {
            int irank;
            for (irank = 0; irank < size; irank++) {
                printf("  Process %d, ", irank);
                printf("Range: %d - %d", (num_chromo/size)*irank, (num_chromo/size)*(irank+1)-1);
                if (left_overs != 0 && irank < left_overs) {
                    printf(", %d", (num_chromo/size)*size + irank);
                }
                printf(" Best Fitness %+4.2f\n", fitness_array[irank]);
                
                fprintf(output, "  Process %d, ", irank);
                fprintf(output, "Range: %d - %d", (num_chromo/size)*irank, (num_chromo/size)*(irank+1)-1);
                if (left_overs != 0 && irank < left_overs) {
                    fprintf(output, ", %d", (num_chromo/size)*size + irank);
                }
                fprintf(output, " Best Fitness %+4.2f\n", fitness_array[irank]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

/**
 * Initialize the first population of chromosomes to random values.
 */
void initialize(chromosome* pop) {
	
    //  For every gene in every chromosome in the population...
    int ichromo;
	int igene;
	for (ichromo = 0; ichromo < num_chromo; ichromo++) {
		for (igene = 0 ; igene < num_genes; igene++) {
            // assign a random value...
			pop[ichromo].gene[igene] = (uint16_t)rand();
		}
        // and clear the chromosome's fitness.
		pop[ichromo].fitness = 0;
	}
}

/**
 * Set up default values for the global variables.
 */
void setup_defaults() {
    file_name = "output.txt";
    generations = 2000;
	num_chromo = 200;
    num_genes = 10;
    rate_cross = .80;
    rate_mutate = .10;
    rate_elite = .05;
}

/**
 * Generate the default values for the rest of
 * the global variables.
 * @ensure num_elite > 0
 */
void setup_generated() {
    srand(time(0));
   	num_elite = num_chromo * rate_elite;
	if (num_elite < 1)
		num_elite = 1;
	num_crossed = num_chromo * rate_cross;
	num_remaining = num_chromo - num_elite - num_crossed;
	num_mutations = num_chromo * rate_mutate;
}

/**
 * Setup global variable default values or override the 
 * defaults with user provided values.
 */
void setup(int argc, char* argv[], int rank){
    setup_defaults();
    int i;
    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            if (argv[i][1] == 'f') {
                file_name = argv[++i];
            } else if (argv[i][1] == 'p') {
                num_chromo = atoi(argv[++i]);
            } else if (argv[i][1] == 'g') {
                generations = atoi(argv[++i]);
            } else if (argv[i][1] == 'c') {
                rate_cross = atof(argv[++i]);
            } else if (argv[i][1] == 'm') {
                rate_mutate = atof(argv[++i]);
            } else if (argv[i][1] == 'e') {
                rate_elite = atof(argv[++i]);
            } else if (argv[i][1] == 'h') {
                if (rank == 0) {
                    print_help();
                }
                MPI_Finalize();
                exit(0);
            } else if (argv[i][1] == 'd') {
                if (rank == 0){
                    setup_generated();
                    print_stats("Default Values");
                }
                MPI_Finalize();
                exit(0);
            } else {
                if (rank == 0) {
                    printf("Invalid flag: %c\n", argv[i][1]);
                }
                MPI_Finalize();
                exit(1);
            }
        } else {
            if (rank == 0) {
                printf("Invalid parameter: %s\n", argv[i]);
            }

            MPI_Finalize();
            exit(1);
        }
    }
    setup_generated();
}

/**
 * Simple Genetic Algorithm
 * 
 * Uses a genetic algorithm to estimate the results of an equation.
 * Each generation involves the preservation of the "elite" chromosomes,
 * i.e., chromosomes with the greatest fitness, and a crossover operation
 * on a percentage of the remaining chromosomes.  Random mutations will
 * occur on non-elite chromosomes.
 * 
 * Results are written to a file.  The first ten generations include 
 * the best fitness score computed by each process.  Every generation
 * has its highest fitness value recorded.
 * 
 * This program use Message Passing Interface (MPI) to distribute
 * computational work between multiple processes.
 */
int main (int argc, char* argv[]) {    
    // Setup MPI
    int rank;
	int size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Setup SGA
    setup(argc, argv, rank);
    if (rank == 0) {
        print_stats("Values");
    }
    
    // Initialize the two populations.
    chromosome* pop = (chromosome*) malloc(sizeof(chromosome) * num_chromo);
    chromosome* new_pop = (chromosome*) malloc(sizeof(chromosome) * num_chromo);
	
    // Setup Master Process
    if (rank == 0) {
        
        // Open Output File.
        if ( !(output = fopen(file_name, "w")) ) {
            printf("Unable to open output file\n");
            exit(1);
        }
        
        // Initialize the Population
    	initialize(pop);
	}
    // For every generation...
    int igeneration;
	for (igeneration = 0; igeneration < generations; igeneration++) {
        // have the master process record the generation...
        if (rank == 0) {
            printf("Generation: %d\n", igeneration);
            fprintf(output, "Generation: %d\n", igeneration);
        }
        
        // use MPI to calculate the fitness values for the generation...
        MPI_Barrier(MPI_COMM_WORLD);
        mpi_calculation(pop, igeneration, rank, size);
        MPI_Barrier(MPI_COMM_WORLD);
        
        // have the master process...
        if (rank == 0) {         
            // sort the population based on its fitness values...
            sort(pop);
            
            // and record the best fitness value for the generation...
            printf("  Best Fit: %+4.2f\n", pop[0].fitness);
            fprintf(output, "Best Fit: %+4.2f\n", pop[0].fitness);
            
            // as long as this is not the last generation...
            if (igeneration + 1 < generations) {
                // move the elite chromosomes to the new population...                
                elite(pop, new_pop);
                                    
                // perform crossover on some of the other chromosomes
                // and move them to the new population...
                crossover(pop, new_pop);
              
                // move the remaining chromosomes to finish off the
                // new population...
                move_remaining(pop, new_pop);

                // perform mutations on the new population...
                mutate(new_pop);	   
                                
                // and, finally, swap which name refers to the old 
                // population and which name refers to the new population.
                chromosome* temp;
                temp = pop;
                pop = new_pop;
                new_pop = temp;
                
            }            
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Clean up our two populations.
    free(pop);
    free(new_pop);
    
    // Clean up the master process.
    if (rank == 0) {
        fclose(output);
    }
	
    MPI_Finalize();
    return 0;
}
