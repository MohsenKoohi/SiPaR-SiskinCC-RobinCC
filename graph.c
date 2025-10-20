#ifndef __GRAPH_C
#define __GRAPH_C

#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <libgen.h>
#include <mpi.h>

#include "paragrapher.h"
#include "mpi.c"
#include "omp.c"

struct mem_block
{
	void* mb_address;
	unsigned long mb_size;
	struct mem_block* mb_next;
};

struct stolen_partition
{
	unsigned long sp_partition_id;
	struct stolen_partition* sp_next;
};

/*
	Each machine contains a subset of `offsets` and a subset of `edges` in the graph

	Size: 8 Bytes IDs per vertex, without weights on edges or vertices.
	On each machine, we have two types of partitions:
		(1) Initial Partitions (IP) which are read from storage and stored as the first_mem_block
*/

struct dll_800_graph
{
	unsigned long vertices_count;
	unsigned long edges_count;

	unsigned long* whole_offsets_list;          // used by rank 0
	unsigned long* whole_offsets_list_raw;      // used by rank 0
	unsigned long whole_offsets_list_raw_size;  // used by rank 0

	unsigned long offsets_list_first_vertex_ID;  // Inclusive
	unsigned long offsets_list_last_vertex_ID;    // The vertex is excluded but offsets_list has its offset value
	unsigned long* offsets_list;
	unsigned long offsets_list_size;

	unsigned long partitions_count;                // Total number of partitions of graph (hosted and not hosted)
	unsigned long* partitions_start_vertex_ID;     // ID of the first vertex of the partition. -1 for not hosted partitions
	unsigned long** partition_edges;               // NULL for not hosted partitions
	
	unsigned long IP_start_id;    // inclusive
	unsigned long IP_end_id;      // exclusive
	
	struct mem_block* first_mem_block;
};


#include "partitioning.c"

void print_dll_800_graph(struct dll_800_graph* g)
{
	fprintf(LL_LOG, "\nprint_dll_800_gaph:\n  |V|: %'20lu\n  |E|: %'20lu\n", g->vertices_count, g->edges_count);
	fprintf(LL_LOG, "  First offsets: ");
	for(unsigned long v = 0; v < min(g->offsets_list_last_vertex_ID + 1, 20); v++)
		fprintf(LL_LOG, "%lu, ", g->offsets_list[v]);
	if(g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID > 20)
	{
		fprintf(LL_LOG, "...\n  Last offsets: ... ");
		for(unsigned long v = g->offsets_list_last_vertex_ID - 20; v <= g->offsets_list_last_vertex_ID; v++)
			fprintf(LL_LOG, ", %lu", g->offsets_list[v -  g->offsets_list_first_vertex_ID ]);
	}

	fprintf(LL_LOG, "\n  Initial Partitions (IP): %u-%u\n", g->IP_start_id, g->IP_end_id);

	if(g->first_mem_block != NULL)
	{
		fprintf(LL_LOG, "  IP first edges: ");
		unsigned long ec = g->first_mem_block->mb_size/sizeof(unsigned long);
		for(unsigned long e = 0; e < min(ec, 20); e++)
			fprintf(LL_LOG, "%lu, ", ((unsigned long*)g->first_mem_block->mb_address)[e]);
		if(ec > 20)
		{
			fprintf(LL_LOG, " ...\n  IP last edges: ... ");
			for(unsigned long e = ec - 20; e < ec; e++)
				fprintf(LL_LOG, "%lu, ", ((unsigned long*)g->first_mem_block->mb_address)[e]);
		}
		fprintf(LL_LOG, "\n  Total edges: %'lu\n", ec);
	}

	fprintf(LL_LOG, "\n\n");

	return;
}

void __dll_800_webgraph_callback(paragrapher_read_request* req, paragrapher_edge_block* eb, void* in_offsets, void* in_edges, void* buffer_id, void* in_args)
{
	void** args = (void**) in_args;
	unsigned long* completed_callbacks_count = (unsigned long*)args[0];
	unsigned long* graph_edges = (unsigned long*)args[1];
	unsigned long first_edge_offset = (unsigned long)args[2];

	unsigned long* offsets = (unsigned long*)in_offsets;
	unsigned long ec = offsets[eb->end_vertex] + eb->end_edge - offsets[eb->start_vertex] - eb->start_edge;
	unsigned long dest_off = offsets[eb->start_vertex] + eb->start_edge - first_edge_offset;
	unsigned long* ul_in_edges = (unsigned long*)in_edges;

	// No need to parallelize this loop as multiple instances of this callback are concurrently called by the ParaGrapher 
	for(unsigned long e = 0; e < ec; e++, dest_off++)
		graph_edges[dest_off] = ul_in_edges[e];

	paragrapher_csx_release_read_buffers(req, eb, buffer_id);

	__atomic_add_fetch(completed_callbacks_count, 1UL, __ATOMIC_RELAXED);

	return;
}

/*
	Tags used: 21, 22, 23 

	`flags`:	
		bit 0: do not collect stats on end
		bit 1: do not use ParaGrapher FUSE (pg_fuse)
		bit 2: keep whole_offsets_list on R0
		bit 3: do not call flush_os_cache

	`exec_info`:
		0: total time in nanoseconds
		1: total load imbalance avg in percentage multiplied by 1e4
		2: total load imbalance std dev in percentage multiplied by 1e4
		3: total net usage
		4: loading offsets time
		5: loading offsets load imbalance avg in percentage multiplied by 1e4
		6: loading offsets load imbalance std dev in percentage multiplied by 1e4
		7: loading offsets net usage
*/
struct dll_800_graph* get_dll_800_webgraph(
	char* file_name, char* type, unsigned int* flags, unsigned long partitions_per_core, unsigned long* exec_info
)
{	
	// Opening the graph
		unsigned long t1 = get_nano_time();
		unsigned long n1 = get_net_usage();

		fprintf(LL_LOG, "\n\033[3;34mget_dll_800_webgraph()\033[0;37m\nLoading %s, format: binary, flags: 0x%x \n", file_name, *flags);
		(mw_rank == 0) && printf("\n\033[3;34mget_dll_800_webgraph()\033[0;37m\nLoading %s, format: binary, flags: 0x%x\n", file_name,*flags);
			
		int ret = paragrapher_init();
		assert(ret == 0);

		paragrapher_graph_type pgt;
		if(!strcmp(type, "PARAGRAPHER_CSX_WG_400_AP") || !strcmp(type, "PARAGRAPHER_CSX_WG_800_AP"))
			pgt = PARAGRAPHER_CSX_WG_800_AP;
		else
		{
			assert(0 && "get_dll_800_webgraph() does not work for this type of graph.");
			return NULL;
		}

		// Check if the _offsets.bin file exists
		if(mw_rank == 0)
		{
			char bin_offsets_file[1024 + PATH_MAX];
			sprintf(bin_offsets_file,"%.*s_offsets.bin", PATH_MAX, file_name);

			if(access(bin_offsets_file, F_OK) != 0)
			{

				paragrapher_graph* graph = NULL;
				graph = paragrapher_open_graph(file_name, pgt, NULL, 0);
				assert(graph != NULL);

				unsigned long* offsets = (unsigned long*)paragrapher_csx_get_offsets(graph, NULL, 0, -1UL, NULL, 0);
				assert(offsets != NULL);

				unsigned long vertices_count = 0;
				unsigned long edges_count = 0;	
				void* op_args []= {&vertices_count, &edges_count};
				ret = paragrapher_get_set_options(graph, PARAGRAPHER_REQUEST_GET_VERTICES_COUNT, op_args, 1);
				assert (ret == 0);
				ret = paragrapher_get_set_options(graph, PARAGRAPHER_REQUEST_GET_EDGES_COUNT, op_args + 1, 1);
				assert (ret == 0);
		
				if(offsets[vertices_count] != edges_count)
				{
					printf("Offset file does not match.\n");
					exit(-1);
				}

				// Closing the graph
				paragrapher_csx_release_offsets_weights_arrays(graph, offsets);
				offsets = NULL;
				ret = paragrapher_release_graph(graph, NULL, 0);
				assert(ret == 0);				
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		paragrapher_graph* graph = NULL;
		{
			void** open_args = calloc(10, sizeof(void*));
			assert(open_args != NULL);
			int open_argc = 0;
			if(flags == NULL || ((*flags & 2U) == 0) )
				open_args[open_argc++] = "USE_PG_FUSE";
		
			graph = paragrapher_open_graph(file_name, pgt, open_args, open_argc);
			assert(graph != NULL);

			free(open_args);
			open_args = NULL;
		}
		
		unsigned long vertices_count = 0;
		unsigned long edges_count = 0;	
		{
			void* op_args []= {&vertices_count, &edges_count};

			ret = paragrapher_get_set_options(graph, PARAGRAPHER_REQUEST_GET_VERTICES_COUNT, op_args, 1);
			assert (ret == 0);
			ret = paragrapher_get_set_options(graph, PARAGRAPHER_REQUEST_GET_EDGES_COUNT, op_args + 1, 1);
			assert (ret == 0);
			if(mw_rank == 0)
			{
				printf("  Vertices: %'lu\n",vertices_count);
				printf("  Edges: %'lu\n",edges_count);
			}
			fprintf(LL_LOG,"  Vertices: %'lu\n",vertices_count);
			fprintf(LL_LOG,"  Edges: %'lu\n",edges_count);

			unsigned long min_mem = 8 * edges_count + 8 * vertices_count;
			if(min_mem > total_mem)
			{
				(mw_rank == 0) && printf("  Requires %lu GB which is not available.\n", min_mem/1e9);

				paragrapher_release_graph(graph, NULL, 0);
				graph = NULL;

				return NULL;
			}
	
			// val = 1UL << (unsigned int)(log(edges_count)/log(2) - 3);
			// op_args[0] = &val;
			// ret = paragrapher_get_set_options(graph, PARAGRAPHER_REQUEST_SET_BUFFER_SIZE, op_args, 1);
			// assert (ret == 0);
		}

	// Allocating memory
		struct dll_800_graph* g = calloc(sizeof(struct dll_800_graph), 1);
		assert(g != NULL);
		g->vertices_count = vertices_count;
		g->edges_count = edges_count;

	// Initiating partitions
		if(partitions_per_core == 0)
			partitions_per_core = 128;
		g->partitions_count = total_cores * partitions_per_core;
		g->partitions_start_vertex_ID = calloc(sizeof(unsigned long), g->partitions_count + 1);
		assert(g->partitions_start_vertex_ID != NULL);
		#pragma omp parall for
		for(unsigned long p = 0; p < g->partitions_count; p++)
			g->partitions_start_vertex_ID[p] = -1UL;

		g->IP_start_id = partitions_per_core * lower_ranks_cores;
		g->IP_end_id = partitions_per_core * (lower_ranks_cores + pe->threads_count);
		if(mw_rank == mw_size - 1)
			g->IP_end_id = g->partitions_count;

		g->partition_edges = calloc(sizeof(unsigned long*), g->partitions_count);
		assert(g->partition_edges != NULL);

		printf("  Rank %u, partition ranges: %lu - %lu / %lu.\n", mw_rank, g->IP_start_id, g->IP_end_id, g->partitions_count);
		fprintf(LL_LOG, "  Rank %u, partition ranges: %lu - %lu / %lu.\n", mw_rank, g->IP_start_id, g->IP_end_id, g->partitions_count);

	// Reading offsets & transferring
	{
		unsigned long tlo = - get_nano_time();
		unsigned long nlo = - get_net_usage();
		
		if(mw_rank == 0)
		{
			fprintf(LL_LOG, "  Loading offsets from storage & transferring.\n");
			printf("  Loading offsets from storage & transferring.\n");
		}

		// Reading by R0 and partitioning
		if(mw_rank == 0)
		{
			g->whole_offsets_list_raw_size = sizeof(unsigned long) * (1 + g->vertices_count);
			g->whole_offsets_list_raw = numa_alloc_interleaved(g->whole_offsets_list_raw_size);
			assert(g->whole_offsets_list_raw != NULL);
			g->whole_offsets_list = g->whole_offsets_list_raw;

			unsigned long* offsets = (unsigned long*)paragrapher_csx_get_offsets(graph, NULL, 0, -1UL, NULL, 0);
			assert(offsets != NULL);

			#pragma omp parallel for
			for(unsigned long v = 0; v <= vertices_count; v++)
				g->whole_offsets_list[v] = offsets[v];

			paragrapher_csx_release_offsets_weights_arrays(graph, offsets);
			offsets = NULL;

			// Partitioning 
			dll_800_parallel_edge_partitioning(g, g->partitions_start_vertex_ID, g->partitions_count);
			
			// Sending IP_start_id and IP_end_id IDs to nodes
			unsigned long used_partitions = nodes_cores[0] * partitions_per_core;
			for(int n = 1; n < mw_size; n++)
			{
				unsigned long start = used_partitions;
				unsigned long end = used_partitions + nodes_cores[n] * partitions_per_core;
				if(n == mw_size - 1)
					end = g->partitions_count;

				assert(end - start + 1 < (1UL<<31) );

				int mret = MPI_Send(g->partitions_start_vertex_ID + start, end - start + 1,  MPI_UNSIGNED_LONG, n, 21, MPI_COMM_WORLD);
				assert(mret == MPI_SUCCESS);

				used_partitions = end;
			}
		}
		else
		{
			MPI_Status status;
			int mret = MPI_Recv(
				g->partitions_start_vertex_ID + g->IP_start_id, 
				g->IP_end_id - g->IP_start_id + 1, 
				MPI_UNSIGNED_LONG, 0, 21, MPI_COMM_WORLD, &status
			);
			assert(mret == MPI_SUCCESS);

			int count = 0;
			mret = MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
			assert(mret == MPI_SUCCESS);
			assert(count == g->IP_end_id - g->IP_start_id + 1);
		}

		// Print
			printf("  R%3u, vertex IDs: %15lu - %15lu, |V_%u|: %.2f%%\n", mw_rank, 
				g->partitions_start_vertex_ID[g->IP_start_id], g->partitions_start_vertex_ID[g->IP_end_id], mw_rank,
				100.0 * (g->partitions_start_vertex_ID[g->IP_end_id] - g->partitions_start_vertex_ID[g->IP_start_id])/g->vertices_count 
			);
			fprintf(LL_LOG, "  R%3u, vertex IDs: %15lu - %15lu, |V_%u|: %.2f%%\n", mw_rank, 
				g->partitions_start_vertex_ID[g->IP_start_id], g->partitions_start_vertex_ID[g->IP_end_id], mw_rank,
				100.0 * (g->partitions_start_vertex_ID[g->IP_end_id] - g->partitions_start_vertex_ID[g->IP_start_id])/g->vertices_count 
			);

		// Transferring offsets
		g->offsets_list_first_vertex_ID = g->partitions_start_vertex_ID[g->IP_start_id];
		g->offsets_list_last_vertex_ID = g->partitions_start_vertex_ID[g->IP_end_id];

		if(mw_rank == 0)
		{
			unsigned long used_partitions = nodes_cores[0] * partitions_per_core;
			for(int n = 1; n < mw_size; n++)
			{
				unsigned long start_partition = used_partitions;
				unsigned long end_partition = used_partitions + nodes_cores[n] * partitions_per_core;
				if(n == mw_size - 1)
					end_partition = g->partitions_count;

				unsigned long start_vertex_ID = g->partitions_start_vertex_ID[start_partition];
				unsigned long end_vertex_ID = g->partitions_start_vertex_ID[end_partition];

				if(end_vertex_ID - start_vertex_ID + 1 < (1UL<<31))
				{
					int mret = MPI_Send(
						g->whole_offsets_list + start_vertex_ID, 
						end_vertex_ID - start_vertex_ID + 1, 
						MPI_UNSIGNED_LONG, n, 22, MPI_COMM_WORLD
					);
					assert(mret == MPI_SUCCESS);
				}
				else
				{
					assert(end_vertex_ID - start_vertex_ID + 1 < (1UL<<41));
					unsigned long length = (end_vertex_ID - start_vertex_ID + 1)/1024;
					if((end_vertex_ID - start_vertex_ID + 1) % 1024)
						length++;

					int mret = MPI_Send(
						g->whole_offsets_list + start_vertex_ID, 
						length, 
						MPI_UNSIGNED_LONG_1K, n, 23, MPI_COMM_WORLD
					);
					assert(mret == MPI_SUCCESS);
				}

				used_partitions = end_partition;
			}

			// Setting offsets for R0
				if(*flags & 4U)
				{
					// Keep the whole_offsets_list
					g->offsets_list = g->whole_offsets_list;
					g->offsets_list_size = 0;
				}
				else
				{
					// Copying offsets for hosted vertices and releasing mem of entire offsets
					g->offsets_list_size = sizeof(unsigned long) * (g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1);
					g->offsets_list = numa_alloc_interleaved(g->offsets_list_size);
					assert(g->offsets_list != NULL);
					
					#pragma omp parallel for num_threads(pe->threads_count)
					for(unsigned long v = g->offsets_list_first_vertex_ID; v < g->offsets_list_last_vertex_ID + 1; v++)
						g->offsets_list[v] = g->whole_offsets_list[v];
					
					numa_free(g->whole_offsets_list_raw, g->whole_offsets_list_raw_size);
					g->whole_offsets_list_raw_size = 0;
					g->whole_offsets_list_raw = NULL;
					g->whole_offsets_list = NULL;
				}
		}
		else
		{
			g->offsets_list_size = sizeof(unsigned long) * (g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1)
				+ 8 * 1024;
			g->offsets_list = numa_alloc_interleaved(g->offsets_list_size);
			assert(g->offsets_list != NULL);
			
			MPI_Status status;

			if(g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1 < (1UL<<31))
			{
				int mret = MPI_Recv(
					g->offsets_list, 
					g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1,
					MPI_UNSIGNED_LONG, 0, 22, MPI_COMM_WORLD, &status
				);
				assert(mret == MPI_SUCCESS);

				int count = 0;
				mret = MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
				assert(mret == MPI_SUCCESS);
				assert(count == g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1);
			}
			else
			{
				assert(g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1 < (1UL<<41));
				unsigned long length = (g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1)/1024;
				if((g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1) % 1024)
					length++;

				int mret = MPI_Recv(
					g->offsets_list, 
					length,
					MPI_UNSIGNED_LONG_1K, 0, 23, MPI_COMM_WORLD, &status
				);
				assert(mret == MPI_SUCCESS);

				int count = 0;
				mret = MPI_Get_count(&status, MPI_UNSIGNED_LONG_1K, &count);
				assert(mret == MPI_SUCCESS);
				assert(count == length);
			}
		}

		tlo += get_nano_time();
		nlo += get_net_usage();
		
		// Writing stats
		{
			int mret = MPI_Gather(&tlo, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			unsigned long max_t = 0;
			if(mw_rank == 0)
			{
				double load_imb = get_idle_percentage_std_dev_avg_max(ul_array, mw_size, NULL, NULL, &max_t);
				
				if(exec_info != NULL)
				{
					exec_info[4] = max_t;
					exec_info[5] = load_imb * 1e4;
					exec_info[6] = 0;
				}

				printf("  Loading offsets: Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
				fprintf(LL_LOG, "  Loading offsets: Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
			}

			mret = MPI_Gather(&nlo, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
			{
				assert(ul_array[0] == nlo);

				unsigned long sum = 0;
				for(int r = 0; r < mw_size; r++)
					sum += ul_array[r];
				
				if(exec_info != NULL)
				{
					exec_info[7] = sum;
				}

				printf("  Loading offsets: network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
				fprintf(LL_LOG, "  Loading offsets: network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Allocating memory for the first mem blocks including IPs
		g->first_mem_block = calloc(sizeof(struct mem_block), 1);
		assert(g->first_mem_block != NULL);
		g->first_mem_block->mb_size = sizeof(unsigned long) * 
			( g->offsets_list[g->partitions_start_vertex_ID[g->IP_end_id] - g->offsets_list_first_vertex_ID] -
			g->offsets_list[0] );
		g->first_mem_block->mb_address = numa_alloc_interleaved(g->first_mem_block->mb_size);
		assert(g->first_mem_block->mb_address != NULL);
		g->first_mem_block->mb_next = NULL;

		// printf("Rank: %u; Start reading from %'.3f (seconds)\n", mw_rank,  (get_nano_time() - t1)/1e9);

	// Reading edges
	unsigned long read_edges = 0;
	{
		unsigned long completed_callbacks_count = 0;
		unsigned long first_edge_offset = g->offsets_list[0];
		void* callback_args[] = {
			(void*)&completed_callbacks_count, 
			(void*)g->first_mem_block->mb_address, 
			(void*)first_edge_offset
		};
		paragrapher_edge_block eb;
		eb.start_vertex = g->partitions_start_vertex_ID[g->IP_start_id];
		eb.start_edge=0;
		eb.end_vertex =  g->partitions_start_vertex_ID[g->IP_end_id];
		eb.end_edge= 0;

		paragrapher_read_request* req= paragrapher_csx_get_subgraph(graph, &eb, NULL, NULL, __dll_800_webgraph_callback, (void*)callback_args, NULL, 0);
		assert(req != NULL);

		struct timespec ts = {0, 200 * 1000 * 1000};
		long status = 0;
		unsigned long callbacks_count = 0;
		void* op0_args []= {req, &status};
		void* op1_args []= {req, &read_edges};
		void* op2_args []= {req, &callbacks_count};
		unsigned long next_edge_limit_print = 0;
		do
		{
			nanosleep(&ts, NULL);
			
			ret = paragrapher_get_set_options(graph, PARAGRAPHER_REQUEST_READ_STATUS, op0_args, 2);
			assert (ret == 0);
			ret = paragrapher_get_set_options(graph, PARAGRAPHER_REQUEST_READ_EDGES, op1_args, 2);
			assert (ret == 0);
			if(callbacks_count == 0)
			{
				ret = paragrapher_get_set_options(graph, PARAGRAPHER_REQUEST_READ_TOTAL_CALLBACKS, op2_args, 2);
				assert (ret == 0);
			}

			if(read_edges >= next_edge_limit_print)
			{
				fprintf(LL_LOG,
					"  Loading ..., time: %'.1f (s), status: %'ld, read_edges: %'15lu (%.1f%%), completed callbacks: %'lu/%'lu .\n", 
					(get_nano_time() - t1)/1e9, status, read_edges, 100.0 *read_edges/edges_count, 
					completed_callbacks_count, callbacks_count
				);

				next_edge_limit_print = read_edges + 0.05 * edges_count;
			}
		}
		while(status == 0);

		// printf("  Reading graph finished, status: %'ld, read_edges: %'lu, completed callbacks: %'u/%'u .\n", status, read_edges, completed_callbacks_count, callbacks_count);

		// Waiting for all buffers to be processed
		while(completed_callbacks_count < callbacks_count)
		{
			nanosleep(&ts, NULL);
			// printf("  Waiting for callbacks ..., completed callbacks: %'u/%'u .\n", completed_callbacks_count, callbacks_count);
		}
		// Releasing the req
		paragrapher_csx_release_read_request(req);
		req = NULL;
	}

	// Releasing the paragrapher graph
		ret = paragrapher_release_graph(graph, NULL, 0);
		assert(ret == 0);
		graph = NULL;

	// Setting Partitions
		#pragma omp paralle for
		for(unsigned long p = 0; p < g->partitions_count; p++)
		{
			g->partition_edges[p] = NULL;
			if(p < g->IP_start_id || p >= g->IP_end_id)
				continue;

			g->partition_edges[p] = ((unsigned long*)g->first_mem_block->mb_address) + 
				g->offsets_list[g->partitions_start_vertex_ID[p] - g->offsets_list_first_vertex_ID] -
				g->offsets_list[0];
		}

	// Finalizing
		t1 = get_nano_time() - t1;
		n1 = get_net_usage() - n1;

		fprintf(LL_LOG, "  Rank: %u; Reading %'lu edges (%.2f%%); Completed in %'.3f (seconds), net usage: %'.3f (GB)\n", 
			mw_rank,  read_edges, 100.0 * read_edges / g->edges_count, t1/1e9, n1 /1e9
		);
		printf("  Rank: %u; Reading %'lu edges (%.2f%%); Completed in %'.3f (seconds), net usage: %'.3f (GB)\n", 
			mw_rank,  read_edges, 100.0 * read_edges / g->edges_count, t1/1e9, n1 /1e9
		);

		if(flags == NULL || ((*flags & 1U) == 0) )
		{
			unsigned long total_edges = 0;
			int mret = MPI_Gather(&read_edges, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
			{
				for(int r = 0; r < mw_size; r++)
					total_edges += ul_array[r];
				
				printf("  Total read edges: %'lu / %'lu.\n", total_edges, g->edges_count);
				assert(g->edges_count == total_edges);
			}

			mret = MPI_Gather(&t1, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			unsigned long max_t = 0;
			if(mw_rank == 0)
			{
				double load_imb = get_idle_percentage_std_dev_avg_max(ul_array, mw_size, NULL, NULL, &max_t);
				
				if(exec_info != NULL)
				{
					exec_info[0] = max_t;
					exec_info[1] = load_imb * 1e4;
					exec_info[2] = 0;
				}

				printf("  Graph loading, Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
				fprintf(LL_LOG, "  Graph loading, Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
			}

			mret = MPI_Gather(&n1, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
			{
				assert(ul_array[0] == n1);

				unsigned long sum = 0;
				for(int r = 0; r < mw_size; r++)
					sum += ul_array[r];
				
				if(exec_info != NULL)
				{
					exec_info[3] = sum;
				}

				printf("  Graph loading: network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
				fprintf(LL_LOG, "  Graph loading: network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
			}
		}

	// Printing the first vals in the read graph
		print_dll_800_graph(g);

	// Flush the OS cache
		if((flags == NULL) || ((*flags & 8U) == 0))
			flush_os_cache();	
		
		(mw_rank == 0) && printf("\n");
		fprintf(LL_LOG, "\n");
		fflush(stdout);
		fflush(LL_LOG);

	return g;	
}

void __dll_800_binary_callback(char* in_buffer, unsigned long start_offset, unsigned long length, void** arg)
{
	unsigned char* buffer = (unsigned char*)in_buffer;
	unsigned long* mem = arg[0];
	int bytes_per_vertex_ID = (int)(unsigned long)arg[1];
	unsigned long l_offset = (unsigned long)arg[2];

	// if(mw_rank == 0)
	// {
	// 	printf("rank: %u, bc, start_offset:%lu, length: %lu, bpvid: %u\n", mw_rank, start_offset, length, bytes_per_vertex_ID);
	// 	for(int i = 0; i < 10; i++)
	// 		printf("%u ", (unsigned int)buffer[i]);
	// 	printf("\n\n");
	// 	fflush(stdout);
	// }

	assert(start_offset % bytes_per_vertex_ID == 0);
	assert(length % bytes_per_vertex_ID == 0);

	unsigned long val = 0;
	unsigned long mul = 1;
	unsigned long e = start_offset / bytes_per_vertex_ID;

	for(unsigned long b = start_offset; b < start_offset + length; b++)
	{
		val += (*buffer) * mul;
		buffer++;
		mul *= 256;

		if(b % bytes_per_vertex_ID == (bytes_per_vertex_ID - 1))
		{
			mem[e++ - l_offset] = val;

			val = 0;
			mul = 1;
		}
	}
}

/*
	Tags used: 18, 19, 20 

	`flags`:	
		bit 0: do not collect stats on end
		bit 1: apply average for vertex- and edge-balanced partitions 
		bit 2: keep whole_offsets_list on R0
		bit 3: allow caching (disable O_DIRECT) and do not call flush_os_cache

	`exec_info`:
		0: total time in nanoseconds
		1: total load imbalance avg in percentage multiplied by 1e4
		2: total load imbalance std dev in percentage multiplied by 1e4
		3: total net usage
		4: loading offsets time
		5: loading offsets load imbalance avg in percentage multiplied by 1e4
		6: loading offsets load imbalance std dev in percentage multiplied by 1e4
		7: loading offsets net usage
*/
struct dll_800_graph* get_dll_800_binary(
	char* file_name, unsigned int* flags, unsigned long partitions_per_core, unsigned long* exec_info
)
{	
	// Reading vertices/edges count
		char temp[4096 + PATH_MAX];
		char res[128];
		unsigned long t1 = get_nano_time();
		unsigned long n1 = get_net_usage();
		
		fprintf(LL_LOG, "\n\033[3;34mget_dll_800_binary()\033[0;37m\nLoading %s, format: binary, flags: 0x%x \n", file_name, *flags);
		(mw_rank == 0) && printf("\n\033[3;34mget_dll_800_binary()\033[0;37m\nLoading %s, format: binary, flags: 0x%x\n", file_name,*flags);

		unsigned long vertices_count = 0;
		unsigned long edges_count = 0;	
		unsigned int bytes_per_vertex_ID = 0;
		{
			sprintf(temp, "cat %s_props.txt | grep vertices-count | cut -f2 -d:", file_name);
			long ret = run_command(temp, res, 128);
			assert((int)ret  == 0);
			vertices_count = atol(res);
			assert(vertices_count < (1UL << 42)  && "Limiting |V| to 4 trillion.");

			sprintf(temp, "cat %s_props.txt | grep edges-count | cut -f2 -d:", file_name);
			ret = run_command(temp, res, 128);
			assert((int)ret  == 0);
			edges_count = atol(res);

			sprintf(temp, "cat %s_props.txt | grep bytes-per-vert | cut -f2 -d:", file_name);
			ret = run_command(temp, res, 128);
			assert((int)ret  == 0);
			bytes_per_vertex_ID = atoi(res);

			if(mw_rank == 0)
			{
				printf("  Vertices: %'lu\n",vertices_count);
				printf("  Edges: %'lu\n",edges_count);
				printf("  Bytes per vertex ID: %u\n",bytes_per_vertex_ID);
			}
			fprintf(LL_LOG,"  Vertices: %'lu\n",vertices_count);
			fprintf(LL_LOG,"  Edges: %'lu\n",edges_count);
			fprintf(LL_LOG,"  Bytes per vertex ID: %u\n",bytes_per_vertex_ID);


			unsigned long min_mem = 8 * edges_count + 8 * vertices_count;
			if(min_mem > total_mem)
			{
				(mw_rank == 0) && printf("  Requires %lu GB which is not available.\n", min_mem/1e9);
				return NULL;
			}
		}

	// Allocating memory
		struct dll_800_graph* g = calloc(sizeof(struct dll_800_graph), 1);
		assert(g != NULL);
		g->vertices_count = vertices_count;
		g->edges_count = edges_count;

	// Initiating partitions
		if(partitions_per_core == 0)
			partitions_per_core = 128;
		g->partitions_count = total_cores * partitions_per_core;
		g->partitions_start_vertex_ID = calloc(sizeof(unsigned long), g->partitions_count + 1);
		assert(g->partitions_start_vertex_ID != NULL);
		#pragma omp parall for
		for(unsigned long p = 0; p < g->partitions_count; p++)
			g->partitions_start_vertex_ID[p] = -1UL;

		g->IP_start_id = partitions_per_core * lower_ranks_cores;
		g->IP_end_id = partitions_per_core * (lower_ranks_cores + pe->threads_count);
		if(mw_rank == mw_size - 1)
			g->IP_end_id = g->partitions_count;

		g->partition_edges = calloc(sizeof(unsigned long*), g->partitions_count);
		assert(g->partition_edges != NULL);

		printf("  Rank %u, partition ranges: %lu - %lu / %lu.\n", mw_rank, g->IP_start_id, g->IP_end_id, g->partitions_count);
		fprintf(LL_LOG, "  Rank %u, partition ranges: %lu - %lu / %lu.\n", mw_rank, g->IP_start_id, g->IP_end_id, g->partitions_count);

	// Reading offsets & transferring
	{
		unsigned long tlo = - get_nano_time();
		unsigned long nlo = - get_net_usage();
		
		if(mw_rank == 0)
		{
			fprintf(LL_LOG, "  Loading offsets from storage & transferring.\n");
			printf("  Loading offsets from storage & transferring.\n");
		}

		// Reading by R0 and partitioning
		if(mw_rank == 0)
		{
			g->whole_offsets_list_raw_size = 2 * 4096 + 8 * 1024 + sizeof(unsigned long) * (1 + g->vertices_count);
			g->whole_offsets_list_raw = numa_alloc_interleaved(g->whole_offsets_list_raw_size);
			assert(g->whole_offsets_list_raw != NULL);
		
			// Retrieving offsets from storage
			{
				sprintf(temp, "%s_offsets.bin", file_name);
				unsigned int prf_flags = 0;
				if((flags != NULL) && (*flags & 8U))
					prf_flags = 1U;
				g->whole_offsets_list = (unsigned long*)par_read_file_flags(
					temp, 0, sizeof(unsigned long) * (1 + g->vertices_count), (char*)g->whole_offsets_list_raw, prf_flags
				);
			}

			// Partitioning 
			dll_800_parallel_edge_partitioning(g, g->partitions_start_vertex_ID, g->partitions_count);
			
			// Sending IP_start_id and IP_end_id IDs to nodes
			if((flags != NULL) && (*flags & 2U) )
			{
				printf("  Applying average for vertex- and edge-balanced partitioning.\n");
				for(unsigned long p = 1; p < g->partitions_count; p++)
				{
					unsigned long nv = p * (g->vertices_count / g->partitions_count);
					unsigned long sv0 = g->partitions_start_vertex_ID[p];
					
					unsigned long sv1 = (nv + sv0) / 2;
					assert(sv1 > g->partitions_start_vertex_ID[p-1]);
					g->partitions_start_vertex_ID[p] = sv1;
				}
				__sync_synchronize();
			}

			unsigned long used_partitions = nodes_cores[0] * partitions_per_core;
			for(int n = 1; n < mw_size; n++)
			{
				unsigned long start = used_partitions;
				unsigned long end = used_partitions + nodes_cores[n] * partitions_per_core;
				if(n == mw_size - 1)
					end = g->partitions_count;

				assert(end - start + 1 < (1UL<<31) );

				int mret = MPI_Send(g->partitions_start_vertex_ID + start, end - start + 1,  MPI_UNSIGNED_LONG, n, 18, MPI_COMM_WORLD);
				assert(mret == MPI_SUCCESS);

				used_partitions = end;
			}
		}
		else
		{
			MPI_Status status;
			int mret = MPI_Recv(
				g->partitions_start_vertex_ID + g->IP_start_id, 
				g->IP_end_id - g->IP_start_id + 1, 
				MPI_UNSIGNED_LONG, 0, 18, MPI_COMM_WORLD, &status
			);
			assert(mret == MPI_SUCCESS);

			int count = 0;
			mret = MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
			assert(mret == MPI_SUCCESS);
			assert(count == g->IP_end_id - g->IP_start_id + 1);
		}

		// Print
			printf("  R%3u, vertex IDs: %15lu - %15lu, |V_%u|: %.2f%%\n", mw_rank, 
				g->partitions_start_vertex_ID[g->IP_start_id], g->partitions_start_vertex_ID[g->IP_end_id], mw_rank,
				100.0 * (g->partitions_start_vertex_ID[g->IP_end_id] - g->partitions_start_vertex_ID[g->IP_start_id])/g->vertices_count 
			);
			fprintf(LL_LOG, "  R%3u, vertex IDs: %15lu - %15lu, |V_%u|: %.2f%%\n", mw_rank, 
				g->partitions_start_vertex_ID[g->IP_start_id], g->partitions_start_vertex_ID[g->IP_end_id], mw_rank,
				100.0 * (g->partitions_start_vertex_ID[g->IP_end_id] - g->partitions_start_vertex_ID[g->IP_start_id])/g->vertices_count 
			);

		// Transferring offsets
		g->offsets_list_first_vertex_ID = g->partitions_start_vertex_ID[g->IP_start_id];
		g->offsets_list_last_vertex_ID = g->partitions_start_vertex_ID[g->IP_end_id];

		if(mw_rank == 0)
		{
			unsigned long used_partitions = nodes_cores[0] * partitions_per_core;
			for(int n = 1; n < mw_size; n++)
			{
				unsigned long start_partition = used_partitions;
				unsigned long end_partition = used_partitions + nodes_cores[n] * partitions_per_core;
				if(n == mw_size - 1)
					end_partition = g->partitions_count;

				unsigned long start_vertex_ID = g->partitions_start_vertex_ID[start_partition];
				unsigned long end_vertex_ID = g->partitions_start_vertex_ID[end_partition];

				if(end_vertex_ID - start_vertex_ID + 1 < (1UL<<31))
				{
					int mret = MPI_Send(
						g->whole_offsets_list + start_vertex_ID, 
						end_vertex_ID - start_vertex_ID + 1, 
						MPI_UNSIGNED_LONG, n, 19, MPI_COMM_WORLD
					);
					assert(mret == MPI_SUCCESS);
				}
				else
				{
					assert(end_vertex_ID - start_vertex_ID + 1 < (1UL<<41));
					unsigned long length = (end_vertex_ID - start_vertex_ID + 1)/1024;
					if((end_vertex_ID - start_vertex_ID + 1) % 1024)
						length++;

					int mret = MPI_Send(
						g->whole_offsets_list + start_vertex_ID, 
						length, 
						MPI_UNSIGNED_LONG_1K, n, 20, MPI_COMM_WORLD
					);
					assert(mret == MPI_SUCCESS);
				}

				used_partitions = end_partition;
			}

			// Setting offsets for R0
				if(*flags & 4U)
				{
					// Keep the whole_offsets_list
					g->offsets_list = g->whole_offsets_list;
					g->offsets_list_size = 0;
				}
				else
				{
					// Copying offsets for hosted vertices and releasing mem of entire offsets
					g->offsets_list_size = sizeof(unsigned long) * (g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1);
					g->offsets_list = numa_alloc_interleaved(g->offsets_list_size);
					assert(g->offsets_list != NULL);
					
					#pragma omp parallel for num_threads(pe->threads_count)
					for(unsigned long v = g->offsets_list_first_vertex_ID; v < g->offsets_list_last_vertex_ID + 1; v++)
						g->offsets_list[v] = g->whole_offsets_list[v];
					
					numa_free(g->whole_offsets_list_raw, g->whole_offsets_list_raw_size);
					g->whole_offsets_list_raw_size = 0;
					g->whole_offsets_list_raw = NULL;
					g->whole_offsets_list = NULL;
				}
		}
		else
		{
			g->offsets_list_size = sizeof(unsigned long) * (g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1)
				+ 8 * 1024;
			g->offsets_list = numa_alloc_interleaved(g->offsets_list_size);
			assert(g->offsets_list != NULL);
			
			MPI_Status status;

			if(g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1 < (1UL<<31))
			{
				int mret = MPI_Recv(
					g->offsets_list, 
					g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1,
					MPI_UNSIGNED_LONG, 0, 19, MPI_COMM_WORLD, &status
				);
				assert(mret == MPI_SUCCESS);

				int count = 0;
				mret = MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
				assert(mret == MPI_SUCCESS);
				assert(count == g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1);
			}
			else
			{
				assert(g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1 < (1UL<<41));
				unsigned long length = (g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1)/1024;
				if((g->offsets_list_last_vertex_ID - g->offsets_list_first_vertex_ID + 1) % 1024)
					length++;

				int mret = MPI_Recv(
					g->offsets_list, 
					length,
					MPI_UNSIGNED_LONG_1K, 0, 20, MPI_COMM_WORLD, &status
				);
				assert(mret == MPI_SUCCESS);

				int count = 0;
				mret = MPI_Get_count(&status, MPI_UNSIGNED_LONG_1K, &count);
				assert(mret == MPI_SUCCESS);
				assert(count == length);
			}
		}

		tlo += get_nano_time();
		nlo += get_net_usage();
		
		// Writing stats
		{
			int mret = MPI_Gather(&tlo, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			unsigned long max_t = 0;
			if(mw_rank == 0)
			{
				double load_imb = get_idle_percentage_std_dev_avg_max(ul_array, mw_size, NULL, NULL, &max_t);
				
				if(exec_info != NULL)
				{
					exec_info[4] = max_t;
					exec_info[5] = load_imb * 1e4;
					exec_info[6] = 0;
				}

				printf("  Loading offsets: Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
				fprintf(LL_LOG, "  Loading offsets: Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
			}

			mret = MPI_Gather(&nlo, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
			{
				assert(ul_array[0] == nlo);

				unsigned long sum = 0;
				for(int r = 0; r < mw_size; r++)
					sum += ul_array[r];
				
				if(exec_info != NULL)
				{
					exec_info[7] = sum;
				}

				printf("  Loading offsets: network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
				fprintf(LL_LOG, "  Loading offsets: network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	// Allocating memory for the first mem blocks including IPs
		g->first_mem_block = calloc(sizeof(struct mem_block), 1);
		assert(g->first_mem_block != NULL);
		g->first_mem_block->mb_size = sizeof(unsigned long) * 
			( g->offsets_list[g->partitions_start_vertex_ID[g->IP_end_id] - g->offsets_list_first_vertex_ID] -
			g->offsets_list[0] );
		g->first_mem_block->mb_address = numa_alloc_interleaved(g->first_mem_block->mb_size);
		assert(g->first_mem_block->mb_address != NULL);
		g->first_mem_block->mb_next = NULL;

		// printf("Rank: %u; Start reading from %'.3f (seconds)\n", mw_rank,  (get_nano_time() - t1)/1e9);

	// Reading edges
		unsigned long read_edges;
		{
			void* args[] = {
				(void*)g->first_mem_block->mb_address, 
				(void*)(unsigned long)bytes_per_vertex_ID,
				(void*)g->offsets_list[0]
			};
				
			sprintf(temp, "%s_edges.bin", file_name);
			unsigned long start_offset = bytes_per_vertex_ID * g->offsets_list[0];
			unsigned long end_offset =  bytes_per_vertex_ID * 
				g->offsets_list[g->partitions_start_vertex_ID[g->IP_end_id] - g->offsets_list_first_vertex_ID];
				
			unsigned int prf_flags = 0;
			if((flags != NULL) && (*flags & 8U))
				prf_flags = 1U;
			par_read_file_invoke_flags(temp, start_offset, end_offset, bytes_per_vertex_ID, __dll_800_binary_callback, args, prf_flags);
			
			read_edges = (end_offset - start_offset)/bytes_per_vertex_ID;
		}

	// Setting Partitions
		#pragma omp paralle for
		for(unsigned long p = 0; p < g->partitions_count; p++)
		{
			g->partition_edges[p] = NULL;
			if(p < g->IP_start_id || p >= g->IP_end_id)
				continue;

			g->partition_edges[p] = ((unsigned long*)g->first_mem_block->mb_address) + 
				g->offsets_list[g->partitions_start_vertex_ID[p] - g->offsets_list_first_vertex_ID] -
				g->offsets_list[0];
		}

	// Finalizing
		t1 = get_nano_time() - t1;
		n1 = get_net_usage() - n1;

		fprintf(LL_LOG, "  Rank: %u; Reading %'lu edges (%.2f%%); Completed in %'.3f (seconds), net usage: %'.3f (GB)\n", 
			mw_rank,  read_edges, 100.0 * read_edges / g->edges_count, t1/1e9, n1 /1e9
		);
		printf("  Rank: %u; Reading %'lu edges (%.2f%%); Completed in %'.3f (seconds), net usage: %'.3f (GB)\n", 
			mw_rank,  read_edges, 100.0 * read_edges / g->edges_count, t1/1e9, n1 /1e9
		);

		if(flags == NULL || ((*flags & 1U) == 0) )
		{
			unsigned long total_edges = 0;
			int mret = MPI_Gather(&read_edges, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
			{
				for(int r = 0; r < mw_size; r++)
					total_edges += ul_array[r];
				
				printf("  Total read edges: %'lu / %'lu.\n", total_edges, g->edges_count);
				assert(g->edges_count == total_edges);
			}

			mret = MPI_Gather(&t1, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			unsigned long max_t = 0;
			if(mw_rank == 0)
			{
				double load_imb = get_idle_percentage_std_dev_avg_max(ul_array, mw_size, NULL, NULL, &max_t);
				
				if(exec_info != NULL)
				{
					exec_info[0] = max_t;
					exec_info[1] = load_imb * 1e4;
					exec_info[2] = 0;
				}

				printf("  Graph loading, Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
				fprintf(LL_LOG, "  Graph loading, Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
			}

			mret = MPI_Gather(&n1, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
			{
				assert(ul_array[0] == n1);

				unsigned long sum = 0;
				for(int r = 0; r < mw_size; r++)
					sum += ul_array[r];
				
				if(exec_info != NULL)
				{
					exec_info[3] = sum;
				}

				printf("  Graph loading: network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
				fprintf(LL_LOG, "  Graph loading: network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
			}
		}

	// Printing the first vals in the read graph
		print_dll_800_graph(g);

	// Flush the OS cache
		if((flags == NULL) || ((*flags & 8U) == 0))
			flush_os_cache();	
		
		(mw_rank == 0) && printf("\n");
		fprintf(LL_LOG, "\n");
		fflush(stdout);
		fflush(LL_LOG);

	return g;
}

void release_numa_interleaved_dll_800_graph(struct dll_800_graph* g)
{
	assert(g != NULL);

	if(g->whole_offsets_list_raw != NULL)
	{
		numa_free(g->whole_offsets_list_raw, g->whole_offsets_list_raw_size);
		g->whole_offsets_list_raw_size = 0;
		g->whole_offsets_list_raw = NULL;
		g->whole_offsets_list = NULL;
	}
	else
	{
		numa_free(g->offsets_list, g->offsets_list_size);
		g->offsets_list_size = 0;
	}
	g->offsets_list = NULL;

	// Releasing mem blocks
	struct mem_block* mb = g->first_mem_block;
	while(mb != NULL)
	{
		struct mem_block* next_mb = mb->mb_next;
		if(mb->mb_address != NULL && mb->mb_size != 0)
		{
			numa_free(mb->mb_address, mb->mb_size);
			mb->mb_address = NULL;
		}

		free(mb);

		mb = next_mb;
	}

	free(g->partitions_start_vertex_ID);
	g->partitions_start_vertex_ID = NULL;

	free(g->partition_edges);
	g->partition_edges = NULL;

	free(g);
	g = NULL;

	return;
}

#endif
