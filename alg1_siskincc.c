/*
	SiskinCC
*/

#include "aux.c"
#include "graph.c"
#include "dist_cc.c"

void* siskincc_distributed_cc(
	struct dll_800_graph* g, unsigned long* exec_info, unsigned long* ccs_p, int graph_is_symmetric, 
	int* process_complete_p, int binom_1_star_0
);

#ifndef __alg1_dont_include_main
int main(int argc, char** args)
{
	// Locale initialization
		setlocale(LC_NUMERIC, "");
		setbuf(stdout, NULL);
		setbuf(stderr, NULL);
		read_env_vars();
		printf("\n");

	// Initializing omp
		initialize_omp_par_env();

	// MPI Init
		mpi_init(argc, args, MPI_THREAD_FUNNELED);
		
	// Loading the grpah		
		unsigned long exec_info[50] = {0}; 
		
		unsigned int flags = 0U;// 1U: do not collect stats; 8U: do not flush_os_cache
		struct dll_800_graph* graph = NULL;
		if(strcmp(LL_INPUT_GRAPH_TYPE,"binary") == 0)
			graph = get_dll_800_binary(LL_INPUT_GRAPH_PATH, &flags, 64, exec_info);
		else if(
				strcmp(LL_INPUT_GRAPH_TYPE,"PARAGRAPHER_CSX_WG_800_AP") == 0  || 
				strcmp(LL_INPUT_GRAPH_TYPE,"PARAGRAPHER_CSX_WG_404_AP") == 0
			)
			graph = get_dll_800_webgraph(LL_INPUT_GRAPH_PATH, LL_INPUT_GRAPH_TYPE, &flags, 64, exec_info);
		else
		{
			 	if(mw_rank == 0)
					printf("Input format should be either binary or PARAGRAPHER_CSX_WG_800_AP.\n");

				MPI_Finalize();
				return -1;
		}
		if(graph == NULL)
			return -2;

		MPI_Barrier(MPI_COMM_WORLD);

		unsigned long ccs_count = 0;
		int process_complete = 0;
	
	// SiskinCC, star
		void* cc = siskincc_distributed_cc(graph, exec_info + 10, &ccs_count, 
			LL_INPUT_GRAPH_IS_SYMMETRIC, &process_complete, 0);  // last arg: binom_1_star_0
		if(process_complete == 0)
			goto report;

		unsigned int* cc_ui = NULL;
		unsigned long* cc_ul = NULL;		
		if(graph->vertices_count < (1UL<<32))
			cc_ui = (unsigned int*)cc;
		else
			cc_ul = (unsigned long*)cc;
		
		MPI_Barrier(MPI_COMM_WORLD);
			
		// Verification
		int validated = 0;
		if(graph->edges_count < 3e10) 
		{
			// Collect CC results if it is distributed

			// Validation
			{
				unsigned long cenol_ccs_count = 0;
				// Alg_6
				void* cenol_cc = dist_cc_jt_cent_non_olap_async(graph, exec_info + 20, &cenol_ccs_count, LL_INPUT_GRAPH_IS_SYMMETRIC);
				unsigned int* cenol_cc_ui = NULL;
				unsigned long* cenol_cc_ul = NULL;
				if(graph->vertices_count < (1UL<<32))
					cenol_cc_ui = (unsigned int*)cenol_cc;
				else
					cenol_cc_ul = (unsigned long*)cenol_cc;

				if(mw_rank == 0 && cenol_cc != NULL)
				{
					if(cenol_cc_ui != NULL)
						#pragma omp parallel for
						for(unsigned long v = 0; v < graph->vertices_count; v++)
						{
							assert(cc_ui[cenol_cc_ui[v]] == cc_ui[v]);
							assert(cenol_cc_ui[cc_ui[v]] == cenol_cc_ui[v]);
						}
					else
						#pragma omp parallel for
						for(unsigned long v = 0; v < graph->vertices_count; v++)
						{
							assert(cc_ul[cenol_cc_ul[v]] == cc_ul[v]);
							assert(cenol_cc_ul[cc_ul[v]] == cenol_cc_ul[v]);
						}

					assert(cenol_ccs_count == ccs_count);

					printf("\n\nValidation completed.\n");
					fflush(stdout);

					if(cenol_cc_ui != NULL)
					{
						numa_free(cenol_cc_ui, sizeof(unsigned int) * graph->vertices_count);
						cenol_cc_ui = NULL;
					}
					else
					{
						numa_free(cenol_cc_ul, sizeof(unsigned long) * graph->vertices_count);
						cenol_cc_ul = NULL;
					}

					validated = 1;
				}
			}
		}

		//Release mem 
		if(cc_ui != NULL)
		{
			numa_free(cc_ui, sizeof(unsigned int) * graph->vertices_count);
			cc_ui = NULL;
		}
		if(cc_ul != NULL)
		{
			numa_free(cc_ul, sizeof(unsigned long) * graph->vertices_count);
			cc_ul = NULL;
		}
		cc = NULL;

	// Writing to the report
		report:
		if(LL_OUTPUT_REPORT_PATH != NULL && mw_rank == 0)
		{
			FILE* out = fopen(LL_OUTPUT_REPORT_PATH, "a");
			if(out == NULL)
				printf("Could not open the report file on %s.\n", LL_OUTPUT_REPORT_PATH);
			else
			{
				if(LL_INPUT_GRAPH_BATCH_ORDER == 0)
				{
					fprintf(out, "%-20s; %-8s; %-8s;", "Dataset", "|V|", "|E|");
					fprintf(out, " %-6s; %-6s; %-9s;", "#Nodes", "#Cores", "LD(s)");
					fprintf(out, " %-10s; %-6s; %-9s;", "Time(s)", "LdImb%", "Net (GB)");
					fprintf(out, " %-12s; %-1s;", "|CCs|", "V");
					fprintf(out, "\n");
				}

				char temp1 [16];
				char temp2 [16];			

				fprintf(out, "%-20s; %-8s; %-8s;", // "Dataset", "|V|", "|E|" 
					LL_INPUT_GRAPH_BASE_NAME, ul2s(graph->vertices_count, temp1), ul2s(graph->edges_count, temp2)
				);
				fprintf(out, " %6u; %6u; %9.1f;",  // "#Nodes", "#Cores" "LD(s)"
					mw_size, total_cores, exec_info[0] / 1e9
				);

				unsigned long* cei = exec_info + 10;
				fprintf(out, " %10.1f; %6.1f; %9.1f;", //  "Binom.sec", "LdImb%", "Net (GB)"
					cei[0] / 1e9, cei[1]/1e4, cei[3]/1e9
				);
				fprintf(out, " %12lu;", //  "|CCs|"
					ccs_count
				);
				fprintf(out, " %1d;", validated); // "Vld"

				fprintf(out, "\n");

				fclose(out);
				out = NULL;
			}
		}

	// Releasing memory
		release_numa_interleaved_dll_800_graph(graph);
		graph = NULL;

	// Finalizing
		int	mret = MPI_Finalize();
		assert(mret == MPI_SUCCESS);

	if(mw_rank == 0)
		printf("\n\n");
	
	return 0;
}
#endif

/*
	MPI tags: 
		50 (became non-root), 
		51 (#buffers for 50 vals), 
*/
struct __siskincc_buffer
{
	unsigned long lock; // protecting writes to any of `allocated_count`, `access_count`, and `effective_blocks_count` fields
	unsigned long allocated_count;
	unsigned long access_count;
	
	unsigned long sender_rank; 

	unsigned long effective_blocks_count; 
	unsigned long last_block_size;
	unsigned long recv_type;   // 50
};

/*
	Send buffer
		Used to send 50 tags
		Statuses: 
			(1) allocated_count < blocks_per_buffer
			(2) allocated_count == blocks_per_buffer && access_count == 0: buffer is going to be sent, wait

	Recv buffer 
		Used to receive 50 tags
		Statuses:
			(1) allocated_count == effective_blocks_count and access_count == 0: buffer does not contain any thing, waiting for receiving
			(2) allocated_count < effective_blocks_count: each thread can allocate a block and process it
*/
struct __siskincc_buffer_group
{
	struct __siskincc_buffer* buffers;
	MPI_Request* requests;
	struct __siskincc_buffer_group* next_group;

	unsigned char* mem;
};

struct __siskincc_buffers
{
	unsigned long send_0_recv_1;

	unsigned long bytes_per_vertex;
	unsigned long vertex_pairs_per_block;
	unsigned long block_size;
	unsigned long blocks_per_buffer;
	unsigned long buffer_size;
	unsigned long buffers_per_group;
	
	unsigned long group_adding_lock;
	unsigned long groups_count;
	
	// For recv
	unsigned long* processed_50_buffers; 

	struct __siskincc_buffer_group* first_group;
};

struct __siskincc_block_ID
{
	struct __siskincc_buffer_group* group;
	struct __siskincc_buffer* buff; 
	unsigned long buff_index;
	unsigned long block_index;
	unsigned long block_size;
};

int __siskincc_group_allocate(struct __siskincc_buffers* buffs)
{
	assert(buffs != NULL);

	if(!__sync_bool_compare_and_swap(&buffs->group_adding_lock, 0UL, 1UL))
		return -1;
	
	struct __siskincc_buffer_group** gp = &buffs->first_group;
	while((*gp) != NULL)
		gp = &((*gp)->next_group);

	struct __siskincc_buffer_group* g = calloc(sizeof(struct __siskincc_buffer_group), 1);
	assert(g != NULL);

	g->buffers = calloc(sizeof(struct __siskincc_buffer), buffs->buffers_per_group);
	assert(g->buffers != NULL);
	g->requests = calloc(sizeof(MPI_Request), buffs->buffers_per_group);
	assert(g->requests != NULL);
	for(unsigned long r = 0; r < buffs->buffers_per_group; r++)
		g->requests[r] = MPI_REQUEST_NULL;
	g->next_group = NULL;
	
	g->mem = numa_alloc_interleaved(buffs->buffer_size * buffs->buffers_per_group);
	assert(g->mem != NULL);

	for(unsigned long b = 0; b < buffs->buffers_per_group; b++)
	{
		g->buffers[b].effective_blocks_count = buffs->blocks_per_buffer;
		g->buffers[b].access_count = 0;
		g->buffers[b].sender_rank = -1UL;
		g->buffers[b].last_block_size = buffs->block_size;
		g->buffers[b].recv_type = -1UL;
		g->buffers[b].lock = 0;

		if(buffs->send_0_recv_1 == 1)
			g->buffers[b].allocated_count = buffs->blocks_per_buffer;  // Prevent accessing by comp thread before receiving something
		else
			g->buffers[b].allocated_count = 0;   // allow access
	}

	*gp = g;
	buffs->groups_count++;
	__sync_synchronize();
	buffs->group_adding_lock = 0;

	return 0;
}

void __siskincc_cont_releae_block(struct __siskincc_buffers* buffs, struct __siskincc_block_ID* prev_block)
{
	assert(buffs != NULL);
	assert(buffs->first_group != NULL);
	assert(prev_block != NULL);
	assert(buffs->send_0_recv_1 == 0);

	// Lock the buffer
	while(1)
	{
		if(__sync_bool_compare_and_swap(&prev_block->buff->lock, 0, 1))
			break;
		usleep(10);
	}

	prev_block->buff->access_count--;

	// Unlock the buffer
	__sync_synchronize();
	prev_block->buff->lock = 0;

	prev_block->group = NULL;
	prev_block->buff = NULL;
	prev_block->buff_index = -1UL;
	prev_block->block_index = -1UL;

	return;
}

unsigned char* __siskincc_get_a_block(struct __siskincc_buffers* buffs, struct __siskincc_block_ID* prev_block)
{
	assert(buffs != NULL);
	assert(buffs->first_group != NULL);
	assert(prev_block != NULL);

	if(prev_block->group == NULL)
	{
		// First access
		prev_block->group = buffs->first_group;
		prev_block->buff = prev_block->group->buffers;
		prev_block->buff_index = 0UL;
		prev_block->block_index = -1UL;
	}
	else
	{
		// Release the block
		// Lock the buffer
		while(1)
		{
			if(__sync_bool_compare_and_swap(&prev_block->buff->lock, 0, 1))
				break;
			usleep(10);
		}

		prev_block->buff->access_count--;

		if(buffs->send_0_recv_1 == 1)
			if(prev_block->buff->access_count == 0 && prev_block->buff->allocated_count == prev_block->buff->effective_blocks_count)
			{
				// A received buffer has been completely processed.
				assert(buffs->processed_50_buffers != NULL);
				assert(prev_block->buff->sender_rank < mw_size);
				assert(prev_block->buff->recv_type == 50);
				__sync_add_and_fetch(buffs->processed_50_buffers + prev_block->buff->sender_rank, 1UL);
			}

		// Unlock the buffer
			__sync_synchronize();
			prev_block->buff->lock = 0;
	}

	unsigned long tries = 1;
	struct __siskincc_buffer_group* g = prev_block->group;
	while(1)
	{
		for(unsigned long rbf = 0; rbf < buffs->buffers_per_group; rbf++)
		{
			unsigned long bfi;
			if(g == prev_block->group)
				bfi = (rbf + prev_block->buff_index) % buffs->buffers_per_group;
			else
				bfi = rbf;

			struct __siskincc_buffer* bf = &g->buffers[bfi];
			
			if(bf->allocated_count == bf->effective_blocks_count)
				continue;

			// Lock the buffer
			while(1)
			{
				if(__sync_bool_compare_and_swap(&bf->lock, 0, 1))
					break;
				usleep(10);
			}

			if(bf->allocated_count >= bf->effective_blocks_count)
			{
				bf->lock = 0;
				continue;
			}

			unsigned long bli = bf->allocated_count;
			bf->allocated_count++;
			bf->access_count++;

			// Unlock the buffer
			__sync_synchronize();
			bf->lock = 0;

			// Update the prev_block
			prev_block->group = g;
			prev_block->buff = bf;
			prev_block->buff_index = bfi;
			prev_block->block_index = bli;
			if(bli + 1 == bf->effective_blocks_count)
				prev_block->block_size = bf->last_block_size;
			else
				prev_block->block_size = buffs->block_size;

			__sync_synchronize();
			return g->mem + bfi * buffs->buffer_size + bli * buffs->block_size;
		}

		struct __siskincc_buffer_group* ng = g->next_group;
		if(ng == NULL) 
			ng = buffs->first_group;
		if(ng != prev_block->group)
		{
			g = ng;
			continue;
		}

		if(buffs->send_0_recv_1 == 1)
		{
			prev_block->group = NULL;
			prev_block->buff = NULL;
			prev_block->buff_index = -1UL;
			prev_block->block_index = -1UL;

			return NULL;
		}

		// Wait for a buffer to be released and try again
		usleep(100);
		g = ng;
		tries++;
		if(tries % 1000 == 0)
			printf("    R%u; Could not assign send buff; tries: %u.\n", mw_rank, tries);
		assert(tries < 1e6 );
	}

	return NULL;
}

void __siskincc_buffer_free(struct __siskincc_buffers* buffs)
{
	assert(buffs != NULL);
	
	if(buffs->processed_50_buffers != NULL)
	{
		free(buffs->processed_50_buffers);
		buffs->processed_50_buffers = NULL;
	}

	struct __siskincc_buffer_group* g = buffs->first_group;
	while(g != NULL)
	{
		struct __siskincc_buffer_group* ng = g->next_group;

		free(g->buffers);
		g->buffers = NULL;
		
		free(g->requests);
		g->requests = NULL;
		
		numa_free(g->mem, buffs->buffer_size * buffs->buffers_per_group);
		g->mem = NULL;

		free(g);

		g = ng;
	}

	free(buffs);
	buffs = NULL;

	return;
}

void* __siskincc_comp_thread(void* in_args)
{
	// Init
		assert(in_args != NULL);
		int min_paritions_per_iteration = 16;
		unsigned long local_t0 = - get_nano_time();
		unsigned int comp_threads_count;
		if(mw_size == 0)
			comp_threads_count = pe->threads_count;
		else
			comp_threads_count = pe->threads_count - 4;

	// Thread pinning
		#pragma omp parallel num_threads(comp_threads_count)
		{
			unsigned tid = omp_get_thread_num();
			assert(tid < comp_threads_count);
			
			cpu_set_t cs;
			CPU_ZERO(&cs);
			CPU_SET(pe->thread2cpu[tid], &cs);
			
			int ret = sched_setaffinity(0, sizeof(cpu_set_t), &cs);
			assert(ret == 0);
		}

		#pragma omp parallel num_threads(comp_threads_count)
		{
			unsigned tid = omp_get_thread_num();
			assert(tid < comp_threads_count);
			
			cpu_set_t cs;
			CPU_ZERO(&cs);
			int ret=sched_getaffinity(0, sizeof(cpu_set_t), &cs);
			if(ret)
			{
				fprintf(LL_LOG, "Can't get the affinity, %d, %s\n", errno, strerror(errno));
				assert(ret == 0);
			}

			for(unsigned int i = 0; i < CPU_SETSIZE; i++)
				if(CPU_ISSET(i, &cs))
					assert(i == pe->thread2cpu[tid]);
		}

	// Args	
		void** args = (void**)in_args;
		struct dll_800_graph* g = (struct dll_800_graph*) args[0];
		struct __siskincc_buffers* sbuffs = (struct __siskincc_buffers*) args[1];
		struct __siskincc_buffers* rbuffs = (struct __siskincc_buffers*) args[2];
		unsigned int* cc_ui = NULL;
		unsigned long* cc_ul = NULL;
		if(g->vertices_count < (1UL<<32))
			cc_ui = (unsigned int*)args[3];
		else
			cc_ul = (unsigned long*)args[3];
		unsigned long graph_is_symmetric = (unsigned long)args[4];
		unsigned long* total_recv_p = (unsigned long*)args[5];
		unsigned long t0 = (unsigned long)args[6];
		unsigned long* comp_concluded = (unsigned long*)args[7];
		unsigned long* nbsc = (unsigned long*)args[8];
		int comm_children_count = (int)(unsigned long)(unsigned long*)args[9];
		int* comm_children = (int*)args[10];
		// args[11];
		int comm_parent = (int)(unsigned long)(unsigned long*)args[12];
		// args[13];
		//args[14];
		// args[15];
		// args[16];
		// args[17];
		unsigned long* r_counters = (unsigned long*)args[18];
		//args[19];
		unsigned long* steps_info = (unsigned long*)args[20];

	// Variables and memory allocation
		unsigned int bytes_per_vertex = rbuffs->bytes_per_vertex;
		unsigned long partitions_count = g->IP_end_id - g->IP_start_id;
		struct ul_dynamic_partitioning* sdp = ul_dynamic_partitioning_initialize(pe, partitions_count);
			
		unsigned long* ttimes = calloc(sizeof(unsigned long), pe->threads_count);
		assert(ttimes != NULL);

	// Traversing edges & processing receieved buffers
		unsigned long last_report_timestamp = get_nano_time();
		unsigned long mt = - last_report_timestamp;
		unsigned long completed_parts = 0;
		unsigned long total_changed_count = 0;
		
		#pragma omp parallel num_threads(comp_threads_count) reduction(+:total_changed_count)
		{
			// Thread vars
				unsigned int tid = omp_get_thread_num();
				ttimes[tid] = - get_nano_time();
				unsigned long t_partition = -1UL;		
				unsigned long thread_changed_count = 0;

				struct __siskincc_block_ID* rbid = calloc(sizeof(struct __siskincc_block_ID), 1);
				assert(rbid != NULL);
				rbid->group = NULL;
				
				struct __siskincc_block_ID* sbid = NULL;
				unsigned char* sb_buff = NULL;
				unsigned long sb_written_bytes = -1UL;
				if(mw_rank != 0)
				{
					sbid = calloc(sizeof(struct __siskincc_block_ID), 1);
					assert(sbid != NULL);
					sbid->group = NULL;
				}

			while(1)
			{
				int doing_nothing = 1;
				int op = 0;

				// Report status
				// if(0)
				{
					unsigned long rt = get_nano_time();
					unsigned long lrt = last_report_timestamp;
					if((rt - lrt >= 60e9) && __sync_bool_compare_and_swap(&last_report_timestamp, lrt, rt) )
					{
						char temp[1024] = {0};
						if(rbuffs->processed_50_buffers != NULL)
							for(int r=0; r<mw_size; r++)
								if(rbuffs->processed_50_buffers[r])
									sprintf(temp+strlen(temp), "R%d: %ld/%lu/%lu, ", r, nbsc[r], r_counters[r], rbuffs->processed_50_buffers[r]);

						char temp2[1024] = {0};
						{
							unsigned long total_buffs = 0;
							unsigned long null_buffs = 0;
							unsigned long incomplete_rbuffs = 0;
							unsigned long accessing_rbuffs = 0;

							struct __siskincc_buffer * lb = NULL;

							struct __siskincc_buffer_group* g = rbuffs->first_group;
							while(g != NULL)
							{
								for(unsigned long bfi = 0; bfi < rbuffs->buffers_per_group; bfi++)
								{
									total_buffs++;

									if(g->requests[bfi] != MPI_REQUEST_NULL)
										null_buffs++;
									else
									{
										if(g->buffers[bfi].allocated_count < g->buffers[bfi].effective_blocks_count)
										{
											incomplete_rbuffs++;
											lb = &g->buffers[bfi];
										}

										if(g->buffers[bfi].access_count > 0)
											accessing_rbuffs++;
									}
								}

								g = g->next_group;		
							}

							sprintf(temp2, "tb: %lu; nb: %lu; ib: %lu; ab: %lu", total_buffs, null_buffs, incomplete_rbuffs, accessing_rbuffs);
							if(lb != NULL)
								sprintf(temp2 + strlen(temp2), "|%lu,%lu,%lu,%lu,%lu", lb->allocated_count, lb->access_count, lb->sender_rank, 
									lb->effective_blocks_count, lb->last_block_size);
						}
						
						printf("  \033[3;%um R%u,rep@%.1f\033[0;37m, tid:%3u,  sdp_r:%lu/%lu, children: %s; buffers: %s\n", 
							mw_rank + 31, mw_rank, (rt + t0)/1e9, tid, 
							sdp->partitions_remained, sdp->partitions_count, temp, temp2
						);
					}
				}

				// Process own edges
					while(completed_parts < partitions_count && op++ < min_paritions_per_iteration)
					{
						t_partition = ul_dynamic_partitioning_get_next_partition(sdp, tid, t_partition);
						if(t_partition == -1UL)
							break;
					
						// Processing the partition
						{	
							unsigned long partition = g->IP_start_id + t_partition;

							unsigned long edges_offset = 
								g->offsets_list[g->partitions_start_vertex_ID[partition] - g->offsets_list_first_vertex_ID];

							for(unsigned long v = g->partitions_start_vertex_ID[partition]; v < g->partitions_start_vertex_ID[partition + 1]; v++)
							{
								unsigned long offset_s = g->offsets_list[v - g->offsets_list_first_vertex_ID];
								unsigned long offset_e = g->offsets_list[v + 1 - g->offsets_list_first_vertex_ID];

								for(unsigned long e = offset_s; e < offset_e; e++)
								{
									unsigned long neighbor = g->partition_edges[partition][e - edges_offset];
									
									if(graph_is_symmetric)
										if(neighbor >= v)
											break;

									unsigned long changed_from = -1UL;
									unsigned long changed_to = -1UL;

									if(cc_ui != NULL)
									{
										unsigned int x = v;
										unsigned int y = neighbor;

										while(1)
										{
											while(cc_ui[x] != x)
												x = cc_ui[x];

											while(cc_ui[y] != y)
												y = cc_ui[y];
											
											if(x == y)
												break;

											if(x < y)
											{
												if(__sync_bool_compare_and_swap(&cc_ui[y], y, x))
												{
													changed_from = y;
													changed_to = x;
													break;
												}
											}
											else
											{
												if(__sync_bool_compare_and_swap(&cc_ui[x], x, y))
												{
													changed_from = x;
													changed_to = y;
													break;
												}
											}
										}
									}
									else
									{
										unsigned long x = v;
										unsigned long y = neighbor;
										
										while(1)
										{
											while(cc_ul[x] != x)
												x = cc_ul[x];

											while(cc_ul[y] != y)
												y = cc_ul[y];

											if(x == y)
												break;

											if(x < y)
											{
												if(__sync_bool_compare_and_swap(&cc_ul[y], y, x))
												{
													changed_from = y;
													changed_to = x;
													break;
												}
											}
											else
											{
												if(__sync_bool_compare_and_swap(&cc_ul[x], x, y))
												{
													changed_from = x;
													changed_to = y;
													break;
												}
											}
										}
									}

									// Write the changed vals to be sent to the parent
									if(changed_from != -1UL && sbid != NULL)
									{
										thread_changed_count++;

										if(sb_buff == NULL || sb_written_bytes == sbid->block_size)
										{
											sb_buff = __siskincc_get_a_block(sbuffs, sbid);
											assert(sb_buff != NULL);

											sb_written_bytes = 0;
										}

										unsigned long val = changed_from;	
										for(int b = 0; b < sbuffs->bytes_per_vertex; b++)
										{
											sb_buff[sb_written_bytes++] = (unsigned char)(val & 255);
											val = (val >> 8);
										}

										val = changed_to;	
										for(int b = 0; b < sbuffs->bytes_per_vertex; b++)
										{
											sb_buff[sb_written_bytes++] = (unsigned char)(val & 255);
											val = (val >> 8);
										}
									}
								}
							}

							doing_nothing = 0;
						}

						// Update stats
						{
							unsigned long total = __sync_add_and_fetch(&completed_parts, 1UL);
							assert(total <= partitions_count);
							if(total == partitions_count)
							{
								unsigned long gnt = t0 + get_nano_time();
								fprintf(LL_LOG, "  R%u finished processing own edges, @%.3f seconds.\n", 
									mw_rank, gnt/1e9);
								printf("  R%u finished processing own edges, @%.3f seconds.\n", 
									mw_rank, gnt/1e9);

								steps_info[0] = gnt;
							}
						}
					}

				// Process received buffers
					while(rbuffs->first_group != NULL)
					{
						unsigned char* buff = __siskincc_get_a_block(rbuffs, rbid);
						if(buff == NULL)
							break;

						assert(rbid->buff != NULL);
						assert(rbid->buff->recv_type == 50);
						assert(rbid->block_size % (2 * rbuffs->bytes_per_vertex) == 0);

						for(unsigned long c = 0; c < rbid->block_size; c += 2 * rbuffs->bytes_per_vertex)
						{
							unsigned long ul_x = 0;
							int all_255 = 1;															
							for(int i = 0; i < rbuffs->bytes_per_vertex; i++)
							{
								if(buff[c + i] != 255)
									all_255 = 0;
								ul_x += (((unsigned long)buff[c + i]) << (8 * i));
							}
							if(all_255)
								continue;

							unsigned long ul_y = 0;
							for(int i = 0; i < rbuffs->bytes_per_vertex; i++)
								ul_y += (((unsigned long)buff[c + rbuffs->bytes_per_vertex + i]) << (8 * i));

							assert(ul_x < g->vertices_count);
							assert(ul_y < ul_x);

							unsigned long changed_from = -1UL;
							unsigned long changed_to = -1UL;

							if(cc_ui != NULL)
							{
								unsigned int x = (unsigned int)ul_x;
								unsigned int y = (unsigned int)ul_y;

								while(1)
								{
									while(cc_ui[x] != x)
										x = cc_ui[x];

									while(cc_ui[y] != y)
										y = cc_ui[y];
									
									if(x == y)
										break;

									if(x < y)
									{
										if(__sync_bool_compare_and_swap(&cc_ui[y], y, x))
										{
											changed_from = y;
											changed_to = x;
											break;
										}
									}
									else
									{
										if(__sync_bool_compare_and_swap(&cc_ui[x], x, y))
										{
											changed_from = x;
											changed_to = y;
											break;
										}
									}
								}
							}
							else
							{
								unsigned long x = ul_x;
								unsigned long y = ul_y;
								
								while(1)
								{
									while(cc_ul[x] != x)
										x = cc_ul[x];

									while(cc_ul[y] != y)
										y = cc_ul[y];

									if(x == y)
										break;

									if(x < y)
									{
										if(__sync_bool_compare_and_swap(&cc_ul[y], y, x))
										{
											changed_from = y;
											changed_to = x;
											break;
										}
									}
									else
									{
										if(__sync_bool_compare_and_swap(&cc_ul[x], x, y))
										{
											changed_from = x;
											changed_to = y;
											break;
										}
									}
								}
							}

							// Write the changed vals to be sent to the parent
							if(changed_from != -1UL && sbid != NULL)
							{
								thread_changed_count++;

								if(sb_buff == NULL || sb_written_bytes == sbid->block_size)
								{
									sb_buff = __siskincc_get_a_block(sbuffs, sbid);
									assert(sb_buff != NULL);

									sb_written_bytes = 0;
								}

								unsigned long val = changed_from;	
								for(int b = 0; b < sbuffs->bytes_per_vertex; b++)
								{
									sb_buff[sb_written_bytes++] = (unsigned char)(val & 255);
									val = (val >> 8);
								}

								val = changed_to;	
								for(int b = 0; b < sbuffs->bytes_per_vertex; b++)
								{
									sb_buff[sb_written_bytes++] = (unsigned char)(val & 255);
									val = (val >> 8);
								}
							}
						}
						
						doing_nothing = 0;
					}

				// The Conclusion Procedure
					// All children processed?
						int processed_all_children = 1;
						
						if(comm_children_count)
							for(int c = 0; c < comm_children_count; c++)
							{
								int cr = comm_children[c];
								assert(cr < mw_size);

								if(nbsc[cr] == -1UL) // The child has not declared it has finished
								{
									processed_all_children = 0;
									break;
								}

								if(r_counters[cr] < nbsc[cr])   // We have not received all messages
								{
									processed_all_children = 0;
									break;
								}

								if(rbuffs->processed_50_buffers[cr] < nbsc[cr])  // There exists some unprocessed buffers
								{
									processed_all_children = 0;
									break;
								}
							}

					if(processed_all_children == 1 && completed_parts == partitions_count)
						break;

					if(doing_nothing == 1)
						usleep(100);
			}
			
			ttimes[tid] += get_nano_time();
			total_changed_count += thread_changed_count;

			free(rbid);
			rbid = NULL;

			if(sbid != NULL)
			{
				// Fill the rest of buffer with ff
				if(sb_buff != NULL)
				{
					for(; sb_written_bytes < sbid->block_size; sb_written_bytes++)
						sb_buff[sb_written_bytes] = (unsigned char)255;

					__siskincc_cont_releae_block(sbuffs, sbid);
				}

				free(sbid);
				sbid = NULL;
			}
		}
		mt += get_nano_time();
		PTIP("  Comp, Traversing edges");
		printf("  R%u; finished processing edges and input buffers @%.3f seconds, |sent CC vals|: %'lu\n", 
			mw_rank, (get_nano_time() + t0)/1e9, total_changed_count);

		if(comm_children_count)
		{
			printf("    R%u; received buffers (rank:sent|received|processed): ", mw_rank);
			for(int c = 0; c < comm_children_count; c++)
			{
				int cr = comm_children[c];
				printf("%d:%lu|%lu|%lu; ", cr, nbsc[cr], r_counters[cr], rbuffs->processed_50_buffers[cr]);
			}
			printf("\n");
		}
		ul_dynamic_partitioning_reset(sdp);
		steps_info[1] = get_nano_time() + t0;
		{
			unsigned long total_recv_buffs = 0;
			for(int c = 0; c < comm_children_count; c++)
				total_recv_buffs += nbsc[comm_children[c]];

			steps_info[2] = total_recv_buffs * rbuffs->buffer_size;
		}

	// Pointer Jumping	
		if(mw_rank == 0)
		{
			mt = - get_nano_time();
			#pragma omp parallel num_threads(pe->threads_count) 
			{
				unsigned int tid = omp_get_thread_num();
				ttimes[tid] = - get_nano_time();
					
				unsigned long partition = -1UL;
				while(1)
				{
					partition = ul_dynamic_partitioning_get_next_partition(sdp, tid, partition);
					if(partition == -1UL)
						break; 
					
					unsigned long v_start = (g->vertices_count / partitions_count) * partition;
					unsigned long v_end = (g->vertices_count / partitions_count) * (1 + partition);
					if(partition == partitions_count - 1)
						v_end = g->vertices_count;
					
					for(unsigned long v = v_start; v < v_end; v++)
					{
						if(cc_ui != NULL)
						{
							while(cc_ui[cc_ui[v]] != cc_ui[v])
								cc_ui[v] = cc_ui[cc_ui[v]];
						}
						else
						{
							while(cc_ul[cc_ul[v]] != cc_ul[v])
								cc_ul[v] = cc_ul[cc_ul[v]];
						}
					}
				}
				ttimes[tid] += get_nano_time();	
			}
			mt += get_nano_time();
			PTIP("  Comp, Pointer jumping");
		}
		steps_info[5] = get_nano_time() + t0;

	// Finalzing
		printf("  R%u; comp threads concluded @%'.3f seconds.\n", mw_rank, (get_nano_time() + t0)/1e9);
				
		__sync_synchronize();
		*comp_concluded = 1;

	// Releasing mem
		ul_dynamic_partitioning_release(sdp);
		sdp = NULL;
		
		free(ttimes);
		ttimes = NULL;

	return NULL;
}

void* siskincc_distributed_cc(
	struct dll_800_graph* g, unsigned long* exec_info, unsigned long* ccs_p, int graph_is_symmetric,
	int* process_complete_p, int binom_1_star_0
)
{
	// Initial checks
		unsigned long t0 = - get_nano_time();
		unsigned long n0 = - get_net_usage();
		assert(g != NULL);
		assert(process_complete_p != NULL);
		*process_complete_p = 0;

		fprintf(LL_LOG, "\n\033[3;33msiskincc_distributed_cc\033[0;37m, binom_1_star_0: %d\n", binom_1_star_0);
		(mw_rank == 0) && printf("\n\033[3;33msiskincc_distributed_cc\033[0;37m, binom_1_star_0: %d\n", binom_1_star_0);

	// Allow this thread use all cores
		{
			cpu_set_t cs;
			CPU_ZERO(&cs);
			for(unsigned int i = 0; i < CPU_SETSIZE; i++) 
				CPU_SET(i, &cs);
			int ret = sched_setaffinity(0, sizeof(cpu_set_t), &cs);
			assert(ret == 0);
		}

	// Allocating memory for CC vals
		volatile unsigned long* cc_ul =  NULL;
		volatile unsigned int* cc_ui = NULL;
		if(g->vertices_count < (1UL<<32))
			cc_ui = numa_alloc_interleaved(sizeof(unsigned int) * g->vertices_count);
		else
			cc_ul = numa_alloc_interleaved(sizeof(unsigned long) * g->vertices_count);
		assert(cc_ui != NULL || cc_ul != NULL);
		unsigned long* ttimes = calloc(sizeof(unsigned long), pe->threads_count);
		assert(ttimes != NULL);

		unsigned long partitions_count = g->IP_end_id - g->IP_start_id;
		struct ul_dynamic_partitioning* dp = ul_dynamic_partitioning_initialize(pe, partitions_count);

	// Creation of the binomial tree / start topology
		int comm_parent = -1;
		int* comm_children = calloc(sizeof(int), mw_size);
		assert(comm_children != NULL);
		unsigned int comm_children_count = 0;

		if(binom_1_star_0 == 1)
		{
			int bits_count = ceil(log2(mw_size));

			for(int n = 0; n < bits_count; n++)
			{
				unsigned int bit_val = (1U<<n) & mw_rank;
				
				if(!bit_val)
				{
					int child = mw_rank | (1U<<n);
					if(child < mw_size)
						comm_children[comm_children_count++] = child;
				}

				if(bit_val)
				{
					comm_parent = mw_rank & (-1U << (n+1));
					break;
				}
			}
			
			assert(mw_rank == 0 || comm_parent >= 0);
			assert(comm_parent < mw_rank);

			unsigned int total_child = 0;
			int mret = MPI_Reduce(&comm_children_count, &total_child, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
				assert(total_child == mw_size - 1);
		}
		else
		{
			if(mw_rank == 0)
			{
				for(int r = 1; r < mw_size; r++)
					comm_children[comm_children_count++] = r;
			}
			else
				comm_parent = 0;
		}

		// Printing topo
		{
			(mw_rank == 0) && printf("  Communication topology:\n");
			fprintf(LL_LOG, "  Communication topology:\n");

			char* temp = calloc(1, 128 + (2 + ceil(log10(mw_size))) * mw_size);
			assert(temp != NULL);

			sprintf(temp, "    R%d, parent: %d, children: ", mw_rank, comm_parent);
			for(int c = 0; c < comm_children_count; c++)
				sprintf(temp + strlen(temp), "%d,", comm_children[c]);
			sprintf(temp + strlen(temp), "\n");

			printf(temp);
			fprintf(LL_LOG, temp);

			free(temp);
			temp = NULL;
		}

	// Initializing CC
		unsigned long mt = - get_nano_time();
		#pragma omp parallel num_threads(pe->threads_count)  
		{
			unsigned tid = omp_get_thread_num();
			ttimes[tid] = - get_nano_time();
			unsigned long partition = -1UL;		
			
			while(1)
			{
				partition = ul_dynamic_partitioning_get_next_partition(dp, tid, partition);
				if(partition == -1UL)
					break; 

				unsigned long v_start = (g->vertices_count / partitions_count) * partition;
				unsigned long v_end = (g->vertices_count / partitions_count) * (1 + partition);
				if(partition == partitions_count - 1)
					v_end = g->vertices_count;
				
				if(cc_ui != NULL)
					for(unsigned long v = v_start; v < v_end; v++)
						cc_ui[v] = v;
				else
					for(unsigned long v = v_start; v < v_end; v++)
						cc_ul[v] = v;
			}

			ttimes[tid] += get_nano_time();
		}
		mt += get_nano_time();
		ul_dynamic_partitioning_reset(dp);
		PTIP("  Initializing CC vals");

	// Initializing recv & send buffers
		struct __siskincc_buffers* rbuffs = NULL;
		struct __siskincc_buffers* sbuffs = NULL;
		{
			rbuffs = calloc(sizeof(struct __siskincc_buffers), 1);
			assert(rbuffs != NULL);
			rbuffs->send_0_recv_1 = 1;
			rbuffs->bytes_per_vertex = ceil(log2(g->vertices_count)/8);
			rbuffs->vertex_pairs_per_block = 8UL * 1024;
			rbuffs->block_size = 2 * rbuffs->bytes_per_vertex * rbuffs->vertex_pairs_per_block;
			rbuffs->blocks_per_buffer = 1UL * 1024;
			rbuffs->buffer_size = rbuffs->blocks_per_buffer * rbuffs->block_size;
			assert(rbuffs->buffer_size < (1UL<<31));
			rbuffs->group_adding_lock = 0;
			rbuffs->groups_count = 0;
			rbuffs->first_group = NULL;
			
			unsigned long fm = get_free_mem() - 20*1024UL*1024; // keep 20GB 
			unsigned long vertex_pair_per_buff = rbuffs->vertex_pairs_per_block * rbuffs->blocks_per_buffer;
			unsigned long ideal_recv_buffs_count = comm_children_count * (1 + g->vertices_count / vertex_pair_per_buff);
			unsigned long ideal_mem_for_recv_buff = rbuffs->buffer_size * ideal_recv_buffs_count;
			unsigned long min_mem_for_recv_buffers = min(128 * rbuffs->buffer_size, ideal_mem_for_recv_buff );
			if(comm_children_count == 0)
			{
				min_mem_for_recv_buffers = 0;
				ideal_mem_for_recv_buff = 0;
			}
			unsigned long mem_for_send_buffers = g->vertices_count * 2UL * rbuffs->bytes_per_vertex;
			if(comm_parent == -1)
				mem_for_send_buffers = 0;

			unsigned long can_run = 0;
			if(min_mem_for_recv_buffers + mem_for_send_buffers < fm)
				can_run = 1;
			fprintf(LL_LOG, "  R%u, required min mem: %'lu / %'lu.\n", mw_rank, (min_mem_for_recv_buffers + mem_for_send_buffers)/1024/1024, fm/1024/1024);
			
			int mret = MPI_Allgather(&can_run, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			{
				int cont = 1;
				for(int r = 0; r < mw_size; r++)
				{
					if(ul_array[r] == 0)
					{
						(mw_rank == 0) && printf("Not enough memory on R%u.\n", mw_rank);
						cont = 0;
					}
					if(cont == 0)
						return NULL;
				}
			}
			
			if(comm_children_count > 0)
			{
				rbuffs->buffers_per_group = min(ideal_mem_for_recv_buff, fm - mem_for_send_buffers) / rbuffs->buffer_size;
				assert(rbuffs->buffers_per_group > 0);
				
				int ret = __siskincc_group_allocate(rbuffs);
				assert(ret == 0);
			}
			else
				rbuffs->buffers_per_group = 0;

			rbuffs->processed_50_buffers = calloc(sizeof(unsigned long), mw_size);
			assert(rbuffs->processed_50_buffers != NULL);

			printf(
				"  R%u, Free mem: %'.1f GB; Recv buffs, buffer size: %'lu Bytes; #blocks/buffer: %u, #buffs/group: %u, group size: %'lu MB\n", 
				mw_rank, fm/1024/1024/1024., rbuffs->buffer_size, rbuffs->blocks_per_buffer, rbuffs->buffers_per_group, 
				rbuffs->buffers_per_group * rbuffs->buffer_size / 1024 / 1024
			);

			if(comm_parent != -1)
			{
				sbuffs = calloc(sizeof(struct __siskincc_buffers), 1);
				assert(sbuffs != NULL);
				sbuffs->send_0_recv_1 = 0;
				sbuffs->bytes_per_vertex = rbuffs->bytes_per_vertex;
				sbuffs->vertex_pairs_per_block = rbuffs->vertex_pairs_per_block;
				sbuffs->block_size = rbuffs->block_size;
				sbuffs->blocks_per_buffer = rbuffs->blocks_per_buffer;
				sbuffs->buffer_size = rbuffs->buffer_size;
				assert(sbuffs->buffer_size < (1UL<<31));
				sbuffs->group_adding_lock = 0;
				sbuffs->groups_count = 0;
				sbuffs->first_group = NULL;

				sbuffs->buffers_per_group =  1 + (g->vertices_count * 2UL * sbuffs->bytes_per_vertex) / sbuffs->buffer_size;
				assert(sbuffs->buffers_per_group > 0);

				int ret = __siskincc_group_allocate(sbuffs);
				assert(ret == 0);

				printf("  R%u, Send buffs, #buffs/group: %u, group size: %'lu MB\n", mw_rank, sbuffs->buffers_per_group, 
					sbuffs->buffers_per_group * sbuffs->buffer_size / 1024 / 1024);
			}
		}

	// Creating the computation thread
		volatile unsigned long* r_counters = calloc(sizeof(unsigned long), mw_size);
		assert(r_counters != NULL);
			
		volatile unsigned long comp_concluded = 0;
		unsigned long* nbsc = calloc(sizeof(unsigned long), (unsigned int)mw_size); // Number of buffers sent by children
		assert(nbsc != NULL);
		for(int r = 0; r < mw_size; r++)
			nbsc[r] = -1UL;

		pthread_t comp_tid;
		volatile unsigned long total_recv = 0;
		volatile unsigned long steps_info[16] = {0};
		steps_info[12]= 1000.0 * (g->partitions_start_vertex_ID[g->IP_end_id] - g->partitions_start_vertex_ID[g->IP_start_id])/g->vertices_count;
		steps_info[13] = 1000.0 * (g->offsets_list[g->partitions_start_vertex_ID[g->IP_end_id] - g->offsets_list_first_vertex_ID] - g->offsets_list[0]) / g->edges_count; 
		volatile void* args[32] = {NULL};
		{			
			args[0] = g;
			args[1] = sbuffs;
			args[2] = rbuffs;
			if(cc_ui != NULL)
				args[3] = cc_ui;
			else
				args[3] = cc_ul;
			args[4] = (unsigned long*)((unsigned long)graph_is_symmetric);
			args[5] = &total_recv;
			args[6] = (unsigned long*)t0;
			args[7] = &comp_concluded;
			args[8] = nbsc;
			args[9] = (unsigned long*)((unsigned long)comm_children_count);
			args[10] = comm_children;
			// args[11] 
			args[12] = (unsigned long*)(unsigned long)comm_parent;
			// args[13]
			// args[14]
			// args[15] 
			// args[16] 
			// args[17] 
			args[18] = r_counters;
			// args[19]
			args[20] = steps_info;

			__sync_synchronize();
			int ret = pthread_create(&comp_tid, NULL, __siskincc_comp_thread, (void*)args);
			assert(ret == 0);
		}

	// Send/receive buffers to/from other nodes
		// Vars
			unsigned long send_counter = 0;
			unsigned long send_completed = 0;
			if(comm_parent == -1)
				send_completed = 1;

			MPI_Request* finish_reqs = NULL;
			unsigned long all_children_finished = 0;

		// Initate finish_reqs
			if(comm_children_count)
			{
				finish_reqs = calloc(sizeof(MPI_Request), comm_children_count);
				assert(finish_reqs != NULL);

				for(int c = 0; c < comm_children_count; c++)
				{
					int child_rank = comm_children[c];
					int mret = MPI_Irecv((void*)(&nbsc[child_rank]), 1, MPI_UNSIGNED_LONG, child_rank, 51, MPI_COMM_WORLD, &finish_reqs[c]);
					assert(mret == MPI_SUCCESS);
					assert(finish_reqs[c] != MPI_REQUEST_NULL);
				}
			}
			else
				all_children_finished = 1;

		while(1)
		{
			int doing_nothing = 1;

			// Send CC vals (tag: 50)
			if(send_completed == 0)
			{
				int sc = 1;
				if(comp_concluded == 0)
					sc = 0;

				struct __siskincc_buffer_group* g = sbuffs->first_group;
				while(g != NULL)
				{
					for(unsigned long bfi = 0; bfi < sbuffs->buffers_per_group; bfi++)
					{
						if(g->requests[bfi] != MPI_REQUEST_NULL)
						{
							int flag;
							int mret = MPI_Test(&g->requests[bfi], &flag, MPI_STATUS_IGNORE);
							assert(mret == MPI_SUCCESS);

							if(flag)
							{
								assert(g->requests[bfi] == MPI_REQUEST_NULL);
								
								int l0 = __sync_bool_compare_and_swap(&g->buffers[bfi].lock, 0, 1);
								assert(l0);
								
								g->buffers[bfi].allocated_count = 0;
								assert(g->buffers[bfi].effective_blocks_count == sbuffs->blocks_per_buffer);

								__sync_synchronize();
								g->buffers[bfi].lock = 0;
								doing_nothing = 0;
							}
							else
								sc = 0;
						}
						else if(							
							(g->buffers[bfi].lock == 0) &&
							(comp_concluded == 1 || g->buffers[bfi].allocated_count == g->buffers[bfi].effective_blocks_count) && 
							(g->buffers[bfi].access_count == 0) && (g->buffers[bfi].allocated_count > 0)
						)
						{
							void* b = g->mem + bfi * sbuffs->buffer_size;
							unsigned int len = g->buffers[bfi].allocated_count * sbuffs->block_size;
														
							int mret = MPI_Isend(b, len, MPI_CHAR, comm_parent, 50, MPI_COMM_WORLD, &g->requests[bfi]);
							assert(mret == MPI_SUCCESS);

							send_counter++;
							sc = 0;
							doing_nothing = 0;
						}
					}

					g = g->next_group;
				}

				if(sc == 1)
					send_completed = 1;
			}

			// Receive buffers (tag: 50)
			{
				struct __siskincc_buffer_group* g = rbuffs->first_group;
				while(g != NULL)
				{
					for(unsigned long bfi = 0; bfi < rbuffs->buffers_per_group; bfi++)
					{
						// Check the previously requests
						if(g->requests[bfi] != MPI_REQUEST_NULL)
						{
							MPI_Status status;
							int flag;
							int mret = MPI_Test(&g->requests[bfi], &flag, &status);
							assert(mret == MPI_SUCCESS);

							if(flag)
							{			
								assert(g->requests[bfi] == MPI_REQUEST_NULL);

								int count = 0;
								mret = MPI_Get_count(&status, MPI_CHAR, &count);
								assert(mret == MPI_SUCCESS);
								assert(status.MPI_TAG == 50);
								assert(count <= rbuffs->buffer_size);
								
								g->buffers[bfi].recv_type = status.MPI_TAG;
								g->buffers[bfi].sender_rank = status.MPI_SOURCE;
								
								assert(status.MPI_SOURCE < mw_size);									
								int source_found = 0;
								for(int c = 0; c < comm_children_count; c++)
									if(comm_children[c] == status.MPI_SOURCE)
									{
										source_found = 1;
										break;
									}
								assert(source_found == 1);

								r_counters[status.MPI_SOURCE]++;

								// Lock the buffer
								int l0 = __sync_bool_compare_and_swap(&g->buffers[bfi].lock, 0, 1);
								assert(l0);
								
								g->buffers[bfi].effective_blocks_count = count / rbuffs->block_size;
								if(count % rbuffs->block_size)
								{
									g->buffers[bfi].effective_blocks_count++;
									g->buffers[bfi].last_block_size = count % rbuffs->block_size;
								}
								else
									g->buffers[bfi].last_block_size = rbuffs->block_size;

								g->buffers[bfi].allocated_count = 0;

								__sync_synchronize();
								g->buffers[bfi].lock = 0;

								doing_nothing = 0;
							}
						}
						else if(
							(all_children_finished == 0) &&   // Do not start Irecv if we are in Conclusion Procedure
							(g->buffers[bfi].access_count == 0) && 
							(g->buffers[bfi].lock == 0) &&
							(g->buffers[bfi].allocated_count == g->buffers[bfi].effective_blocks_count)
						)
						{
							void* b = g->mem + bfi * rbuffs->buffer_size;
							g->buffers[bfi].sender_rank = -1UL;
							g->buffers[bfi].recv_type = -1UL;

							int mret = MPI_Irecv(b, rbuffs->buffer_size, MPI_CHAR, MPI_ANY_SOURCE, 50, MPI_COMM_WORLD, &g->requests[bfi]);
							assert(mret == MPI_SUCCESS);
							assert(g->requests[bfi] != MPI_REQUEST_NULL);

							// printf("  R%u, g:%lx, buff: %lu initiated.\n", mw_rank, g, bfi);

							doing_nothing = 0;
						}
					}

					g = g->next_group;		
				}
			}

			// Receive finish announcements (tag:51)
			if(all_children_finished == 0)
			{
				int all_finished = 1;

				for(int c = 0; c < comm_children_count; c++)
				{
					int child_rank = comm_children[c];

					if(finish_reqs[c] != MPI_REQUEST_NULL)
					{
						MPI_Status status;
						int flag = 0;
						int mret = MPI_Test(&finish_reqs[c], &flag, &status);
						assert(mret == MPI_SUCCESS);

						if(flag)
						{				
							assert(finish_reqs[c] == MPI_REQUEST_NULL);

							int count = 0;
							mret = MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
							assert(mret == MPI_SUCCESS);
							assert(count == 1);
							assert(status.MPI_SOURCE == child_rank);
							
							fprintf(LL_LOG, "    R%u; child completion R%u; counter: %'lu; local counter: %'lu @%'.3f seconds.\n", 
								mw_rank, child_rank, nbsc[child_rank], r_counters[child_rank], (get_nano_time()+t0)/1e9);
							printf("    R%u; child completion R%u; counter: %'lu; local counter: %'lu, @%'.3f seconds.\n", 
								mw_rank, child_rank, nbsc[child_rank], r_counters[child_rank], (get_nano_time()+t0)/1e9);

							assert(nbsc[child_rank] >= r_counters[child_rank]);

							doing_nothing = 0;
						}
						else
							all_finished = 0;
					}
				
					if(nbsc[child_rank] != -1UL && nbsc[child_rank] > r_counters[child_rank])
						all_finished = 0;
				}

				if(all_finished)
				{
					all_children_finished = 1;
					printf("  R%u, all children finished, @%'.3f seconds.\n", mw_rank, (get_nano_time()+t0)/1e9);		
					fprintf(LL_LOG, "  R%u, all children finished, @%'.3f seconds.\n", mw_rank, (get_nano_time()+t0)/1e9);		
				}
			} 

			// The Conclusion Procedure
			if(send_completed == 1 && all_children_finished == 1)
			{
				// Inform the parent
					if(mw_rank != 0)
					{
						MPI_Request req = MPI_REQUEST_NULL;
						int mret = MPI_Isend(&send_counter, 1, MPI_UNSIGNED_LONG, comm_parent, 51, MPI_COMM_WORLD, &req);
						assert(mret == MPI_SUCCESS);
						steps_info[4] = send_counter * sbuffs->buffer_size;
						steps_info[3] = t0 + get_nano_time();

						unsigned long ti =  - get_nano_time();
						while(1)
						{
							int flag;
							int mret = MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
							assert(mret == MPI_SUCCESS);

							if(flag)
								break;

							usleep(100);
						}

						unsigned long tt = get_nano_time();
						fprintf(LL_LOG, "  R%u; completion informed to R%u; buffers: %'lu @%'.3f seconds (wait: %.3f s).\n", 
							mw_rank, comm_parent, send_counter, (tt+t0)/1e9, (tt+ti)/1e9);
						printf("  R%u; completion informed to R%u; buffers: %'lu @%'.3f seconds  (wait: %.3f s).\n", 
							mw_rank, comm_parent, send_counter, (tt+t0)/1e9, (tt+ti)/1e9);
					}

				// Cancel MPI recv requests 
					struct __siskincc_buffer_group* g = rbuffs->first_group;
					while(g != NULL)
					{
						for(unsigned long bfi = 0; bfi < rbuffs->buffers_per_group; bfi++)
							if(g->requests[bfi] != MPI_REQUEST_NULL)   // Some buffers may have been released in the current iteration 
							{
								// Ensure all blocks have been processed and no one is accessing them
								assert(g->buffers[bfi].allocated_count == g->buffers[bfi].effective_blocks_count);
								assert(g->buffers[bfi].access_count == 0);

								int mret = MPI_Cancel(&g->requests[bfi]);
								assert(mret == MPI_SUCCESS);

								mret = MPI_Request_free(&g->requests[bfi]);
								assert(mret == MPI_SUCCESS);
								assert(g->requests[bfi] == MPI_REQUEST_NULL);

								doing_nothing = 0;
							}
					
						g = g->next_group;
					}

				break;
			}

			if(doing_nothing)
				usleep(1000);
		}

		pthread_join(comp_tid, NULL);
		
	// Stats
		t0 += get_nano_time();
		n0 += get_net_usage();

		fprintf(LL_LOG, "  R%u; Completed in %'.3f (seconds), net usage: %'.3f (GB)\n", mw_rank,  t0/1e9, n0 /1e9);
		printf("  R%u; Completed in %'.3f (seconds), net usage: %'.3f (GB)\n", mw_rank,  t0/1e9, n0 /1e9);
		{
			int mret = MPI_Gather(&t0, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
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

				printf("\n  Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
				fprintf(LL_LOG, "  Time: %.3f seconds; Load imbalance: %.1f %%\n", max_t/1e9, load_imb);
			}

			mret = MPI_Gather(&n0, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
			{
				assert(ul_array[0] == n0);

				unsigned long sum = 0;
				for(int r = 0; r < mw_size; r++)
					sum += ul_array[r];
				
				sum /= 2;
				if(exec_info != NULL)
				{
					exec_info[3] = sum;
				}

				printf("  Network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
				fprintf(LL_LOG, "  Network usage: %'.3f GB, rate: %'.3f Gbps\n", sum/1e9, 8. * sum / max_t );
			}
		}
		fflush(LL_LOG);

	// Counting number of CCs
		if(mw_rank == 0 && ccs_p != NULL)
		{
			unsigned long ccs = 0;

			#pragma omp parallel for num_threads(pe->threads_count) reduction(+:ccs)
			for(unsigned long v = 0; v < g->vertices_count; v++)
			{
				if(cc_ui != NULL)
				{
					if(cc_ui[v] == v)
						ccs++;
				}
				else
				{
					if(cc_ul[v] == v)
						ccs++;
				}
			}

			printf("  |CCs|: %'lu\n",ccs);
			*ccs_p = ccs;
		}

	// Steps info
		if(mw_rank == 0)
		{
			unsigned long* si = calloc(sizeof(unsigned long), 16 * (1 + mw_size));
			assert(si != NULL);

			int mret = MPI_Gather((void*)steps_info, 16, MPI_UNSIGNED_LONG, si, 16, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			
			printf("\nSteps timing \n");
			char temp[64];

			// Row 0
				printf("%-10s; ", "Rank");
				for(int r = 0; r < mw_size; r++)
				{
					sprintf(temp, "R%u", r);
					printf("%-10s; ", temp);
				}
				printf("\n");

			// Row 1: |V| and |E| 
				printf("%-10s; ", "|V|-|E|");
				for(int r = 0; r < mw_size; r++)
				{
					sprintf(temp, "%.1f-%.1f", si[16 * r + 12]/10.0, si[16 * r + 13]/10.0);
					printf("%-10s; ", temp);
				}
				printf("\n");

			// Rows 2: time finished processing own edges 
				printf("%-10s; ", "OE");
				for(int r = 0; r < mw_size; r++)
				{
					sprintf(temp, "%.1f", si[16 * r + 0]/1e9);
					printf("%-10s; ", temp);
				}
				printf("\n");

			// Rows 3: time finished processing received buffers and its size
				printf("%-10s; ", "RBuff");
				for(int r = 0; r < mw_size; r++)
				{
					sprintf(temp, "%.1f/%.1f", si[16 * r + 1]/1e9, 1.*si[16 * r + 2]/1024/1024/1024);
					printf("%-10s; ", temp);
				}
				printf("\n");

			// Rows 4: time finished sending and its size
				printf("%-10s; ", "SBuff");
				for(int r = 0; r < mw_size; r++)
				{
					sprintf(temp, "%.1f/%.1f", si[16 * r + 3]/1e9, 1.*si[16 * r + 4]/1024/1024/1024);
					printf("%-10s; ", temp);
				}
				printf("\n");

			// Rows 4: time finished pointer jumping
				printf("%-10s; ", "PJ");
				for(int r = 0; r < mw_size; r++)
				{
					sprintf(temp, "%.1f", si[16 * r + 5]/1e9);
					printf("%-10s; ", temp);
				}
				printf("\n");

			free(si);
			si = NULL;
		}
		else
		{
			int mret = MPI_Gather((void*)steps_info, 16, MPI_UNSIGNED_LONG, NULL, 16, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
		}

	// Releasing memory
		free(ttimes);
		ttimes = NULL;
	
		free(comm_children);
		comm_children = NULL;

		free((void*)r_counters);
		r_counters = NULL;

		ul_dynamic_partitioning_release(dp);
		dp = NULL;

		free(nbsc);
		nbsc = NULL;

		if(finish_reqs != NULL)
		{
			free(finish_reqs);
			finish_reqs = NULL;
		}

		if(rbuffs != NULL)
		{
			__siskincc_buffer_free(rbuffs);
			rbuffs = NULL;
		}	

		if(sbuffs != NULL)
		{
			__siskincc_buffer_free(sbuffs);
			sbuffs = NULL;
		}	

	// Return
		*process_complete_p = 1;
		if(mw_rank == 0)
		{
			printf("\n");

			if(cc_ui != NULL)
				return (void*)cc_ui;
			else
				return (void*)cc_ul;
		}
		else
		{
			if(cc_ui != NULL)
				numa_free((void*)cc_ui, sizeof(unsigned int) * g->vertices_count);
			else
				numa_free((void*)cc_ul, sizeof(unsigned long) * g->vertices_count);
			cc_ui = NULL;
			cc_ul = NULL;

			return NULL;
		}
}