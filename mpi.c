#ifndef __MPI_C
#define __MPI_C

#include <mpi.h>
#include "omp.c"

MPI_Datatype MPI_UNSIGNED_LONG_1K;

int mpi_init(int argc, char** args, int mpi_threading_level)
{	
	fprintf(LL_LOG, "MPI Init\n");

	int thread_provided;
	int mret = MPI_Init_thread(&argc, &args, mpi_threading_level, &thread_provided);
	char* val_names[] = {"MPI_THREAD_SINGLE", "MPI_THREAD_FUNNELED", "MPI_THREAD_SERIALIZED", "MPI_THREAD_MULTIPLE"};
	int vals[] = {MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE};
	for(int i = 0; i < sizeof(vals)/sizeof(unsigned int); i++)
		if(thread_provided == vals[i])
			fprintf(LL_LOG, "  Provided MPI threading level: %s\n", val_names[i]);

	mret = MPI_Comm_rank(MPI_COMM_WORLD, &mw_rank);
	assert(mret == MPI_SUCCESS);
	mret = MPI_Comm_size(MPI_COMM_WORLD, &mw_size);
	assert(mret == MPI_SUCCESS);
	(mw_rank == 0) && printf("MPI Init, threading level: %s.\n", val_names[thread_provided]);
	(mw_rank == 0) && printf("  MPI world size (mw_size): %'u\n", mw_size);
		
	char hostname[MPI_MAX_PROCESSOR_NAME + 1];
	hostname[MPI_MAX_PROCESSOR_NAME] = '\0';
	int length;
	mret = MPI_Get_processor_name (hostname, &length);
	assert(mret == MPI_SUCCESS);
	fprintf(LL_LOG, "  Process: %d/%d, started on : %s \n", mw_rank, mw_size, hostname);
	
	int version,subversion;
	mret = MPI_Get_version(&version, &subversion);
	assert(mret == MPI_SUCCESS);
	fprintf(LL_LOG, "  MPI Version: %d.%d\n", version, subversion);
	mw_rank == 0 && printf("  MPI Version: %d.%d\n", version, subversion);

	char cversion[MPI_MAX_LIBRARY_VERSION_STRING + 1];
	cversion[MPI_MAX_LIBRARY_VERSION_STRING] = '\0';
	int resultlen;
	mret = MPI_Get_library_version(cversion, &resultlen);
	assert(mret == MPI_SUCCESS);
	fprintf(LL_LOG, "  MPI Library Version: %s\n", cversion);
	mw_rank == 0 && printf("  MPI Library Version: %s\n", cversion);

	unsigned int* attr_val_pointer = NULL;
	int flag=0;
	mret = MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &attr_val_pointer, &flag);
	assert(mret == MPI_SUCCESS);
	assert(flag != 0);
	mt_upperbound = *attr_val_pointer;
	fprintf(LL_LOG, "  Upper bound for tag value: %'u.\n", *attr_val_pointer);
	mw_rank == 0 && printf("  MPI upper bound for tag value: %'u.\n", *attr_val_pointer);

	// Identifying the number of cores on all machines, used for partitioning
	{
		nodes_cores = calloc(sizeof(unsigned int), mw_size);
		assert(nodes_cores != NULL);

		unsigned int local_cores = pe->threads_count;
		mret = MPI_Allgather(&local_cores, 1, MPI_UNSIGNED, nodes_cores, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
		assert(mret == MPI_SUCCESS);

		fprintf(LL_LOG, "  Node cores: ");
		mw_rank == 0 && printf("  Node cores: ");

		for(int i = 0; i < mw_size; i++)
		{
			if(i < mw_rank)
				lower_ranks_cores += nodes_cores[i];
			total_cores += nodes_cores[i];
		
			fprintf(LL_LOG, "%u, ", nodes_cores[i]);
			mw_rank == 0 && printf("%u, ", nodes_cores[i]);
		}
		
		fprintf(LL_LOG, "\n  Total cores: %'u.\n", total_cores);
		fprintf(LL_LOG, "  Lower ranks cores: %'u.\n", lower_ranks_cores);
		
		mw_rank == 0 && printf("\n  Total cores: %u\n", total_cores);
	}

	ul_array = calloc(sizeof(unsigned long), mw_size);
	assert(ul_array != NULL);

	dl_array = calloc(sizeof(double), mw_size);
	assert(dl_array != NULL);

	// Total mem
	{
		unsigned long tm = get_total_mem();
		mret = MPI_Allgather(&tm, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		assert(mret == MPI_SUCCESS);
		total_mem = 0;
		for(int r = 0; r < mw_size; r++)
			total_mem += ul_array[r];

		(mw_rank == 0) && printf("  Total mem: %'.1f GB\n", 1.0 * total_mem / (1UL << 30));
		fprintf(LL_LOG, "  Total mem: %'.1f GB\n", 1.0 * total_mem / (1UL << 30));
	}

	fprintf(LL_LOG, "\n");
	fflush(LL_LOG);

	if(mw_rank == 0)
	{
		printf( "\n\033[2;32mEnv Vars:\033[0;37m\n");
		printf( "    LL_INPUT_GRAPH_PATH:         %s\n", LL_INPUT_GRAPH_PATH);
		printf( "    LL_INPUT_GRAPH_BASE_NAME:    %s\n", LL_INPUT_GRAPH_BASE_NAME);
		printf( "    LL_INPUT_GRAPH_TYPE:         %s\n", LL_INPUT_GRAPH_TYPE);
		printf( "    LL_INPUT_GRAPH_BATCH_ORDER:  %u\n", LL_INPUT_GRAPH_BATCH_ORDER);
		printf( "    LL_INPUT_GRAPH_IS_SYMMETRIC: %u\n", LL_INPUT_GRAPH_IS_SYMMETRIC);
		printf( "    LL_GRAPH_RA_BIN_FOLDER:      %s\n", LL_GRAPH_RA_BIN_FOLDER);
		printf( "    LL_OUTPUT_REPORT_PATH:       %s\n", LL_OUTPUT_REPORT_PATH);
		printf( "\n");
	}

	MPI_Type_contiguous(1024, MPI_UNSIGNED_LONG, &MPI_UNSIGNED_LONG_1K);
	MPI_Type_commit(&MPI_UNSIGNED_LONG_1K);
	
	return 0;
}



#endif