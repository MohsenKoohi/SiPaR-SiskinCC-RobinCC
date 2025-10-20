#ifndef __DIST_CC_C
#define __DIST_CC_C


/*

 Distributed JT, Central CC vals, Non-overlapping
 	
	Uses MPI tags: 23,24
	`exec_info`:
		0: total time in nanoseconds
		1: total load imbalance avg %
		2: 
		3: total net usage
		4: changed ccs avg %
		5: changed ccs std dev %
		6: effective transferred ccs avg %
		7: effective transferred ccs std dev %
*/

void* dist_cc_jt_cent_non_olap_async(
	struct dll_800_graph* g, unsigned long* exec_info, unsigned long* ccs_p, int graph_is_symmetric
)
{
	// Initial checks
		assert(g != NULL);
		
		unsigned long t0 = - get_nano_time();
		unsigned long n0 = - get_net_usage();
		
		fprintf(LL_LOG, "\n\033[3;33mdist_cc_jt_cent_non_olap_async\033[0;37m\n");
		(mw_rank == 0) && printf("\n\033[3;33mdist_cc_jt_cent_non_olap_async\033[0;37m\n");

	// Allocating memory
		unsigned long* cc_ul = NULL;
		unsigned int* cc_ui = NULL;
		if(g->vertices_count < (1UL<<32))
			cc_ui = numa_alloc_interleaved(sizeof(unsigned int) * g->vertices_count);
		else
			cc_ul = numa_alloc_interleaved(sizeof(unsigned long) * g->vertices_count);
		assert(cc_ul != NULL || cc_ui != NULL);
		unsigned long* ttimes = calloc(sizeof(unsigned long), pe->threads_count);
		assert(ttimes != NULL);
		if(cc_ui != NULL)
		{
			(mw_rank == 0) && printf("  Bytes per CC value: %u\n", sizeof(unsigned int));
			fprintf(LL_LOG, "  Bytes per CC value: %u\n", sizeof(unsigned int));
		}
		else
		{
			(mw_rank == 0) && printf("  Bytes per CC value: %u\n", sizeof(unsigned long));
			fprintf(LL_LOG, "  Bytes per CC value: %u\n", sizeof(unsigned long));
		}

	// Dynamic partitioning
		unsigned long partitions_count = g->IP_end_id - g->IP_start_id;
		struct ul_dynamic_partitioning* dp = ul_dynamic_partitioning_initialize(pe, partitions_count);

	// (1) Initializing
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
		PTIP("  (1) Initializing");
	
	// (2) Traversing edges
		mt = - get_nano_time();
		unsigned long edge_values = 0;
		#pragma omp parallel num_threads(pe->threads_count) reduction(+:edge_values)
		{
			unsigned tid = omp_get_thread_num();
			ttimes[tid] = - get_nano_time();
			unsigned long zb_partition = -1UL;		
			
			while(1)
			{
				zb_partition = ul_dynamic_partitioning_get_next_partition(dp, tid, zb_partition);
				if(zb_partition == -1UL)
					break;

				unsigned long partition = g->IP_start_id + zb_partition;

				unsigned long edges_offset = 
					g->offsets_list[g->partitions_start_vertex_ID[partition] - g->offsets_list_first_vertex_ID];

				for(unsigned long v = g->partitions_start_vertex_ID[partition]; v < g->partitions_start_vertex_ID[partition + 1]; v++)
				{
					unsigned long offset_s = g->offsets_list[v - g->offsets_list_first_vertex_ID];
					unsigned long offset_e = g->offsets_list[v + 1 - g->offsets_list_first_vertex_ID];

					for(unsigned long e = offset_s; e < offset_e; e++)
					{
						unsigned long neighbour = g->partition_edges[partition][e - edges_offset];
						edge_values += neighbour;
						
						if(graph_is_symmetric)
							if(neighbour >= v)
								break;

						if(cc_ui != NULL)
						{
							unsigned int x = v;
							unsigned int y = neighbour;

							while(1)
							{
								while(x != cc_ui[x])
									x = cc_ui[x];

								while(y != cc_ui[y])
									y = cc_ui[y];

								if(x == y)
									break;

								if(x < y)
								{
									if(__sync_bool_compare_and_swap(&cc_ui[y], y, x))
										break;
								}
								else
								{
									if(__sync_bool_compare_and_swap(&cc_ui[x], x, y))
										break;
								}

							}
						}
						else
						{
							unsigned long x = v;
							unsigned long y = neighbour;

							while(1)
							{
								while(x != cc_ul[x])
									x = cc_ul[x];

								while(y != cc_ul[y])
									y = cc_ul[y];

								if(x == y)
									break;

								if(x < y)
								{
									if(__sync_bool_compare_and_swap(&cc_ul[y], y, x))
										break;
								}
								else
								{
									if(__sync_bool_compare_and_swap(&cc_ul[x], x, y))
										break;
								}

							}
						}
					}
				}
			}

			ttimes[tid] += get_nano_time();
		}
		mt += get_nano_time();
		ul_dynamic_partitioning_reset(dp);
		PTIP("  (2) Traversing edges:");

		{
			int mret = MPI_Gather(&edge_values, 1, MPI_UNSIGNED_LONG, ul_array, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			if(mw_rank == 0)
			{
				for(int r = 1; r < mw_size; r++)
					edge_values += ul_array[r];
				printf("  Total edge_values: %'lu\n", edge_values);
				fprintf(LL_LOG, "  Total edge_values: %'lu\n", edge_values);
			}
		}

		
	// (3) Pointer jumping & Counting the number of changes
		mt = - get_nano_time();
		unsigned long changed_cc = 0;
		unsigned int* p2t = calloc(sizeof(unsigned int), partitions_count);
		assert(p2t != NULL);
		unsigned long* thread_changes = calloc(sizeof(unsigned long), pe->threads_count);
		assert(thread_changes != NULL);
		
		#pragma omp parallel num_threads(pe->threads_count) reduction(+: changed_cc)
		{
			unsigned tid = omp_get_thread_num();
			ttimes[tid] = - get_nano_time();
			unsigned long partition = -1UL;		

			unsigned long tid_changes = 0;
			
			while(1)
			{
				partition = ul_dynamic_partitioning_get_next_partition(dp, tid, partition);
				if(partition == -1UL)
					break; 

				p2t[partition] = tid;
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

						if(cc_ui[v] != v)
							tid_changes++;
					}
					else
					{
						while(cc_ul[cc_ul[v]] != cc_ul[v])
							cc_ul[v] = cc_ul[cc_ul[v]];

						if(cc_ul[v] != v)
							tid_changes++;
					}
				}
			}

			thread_changes[tid] = tid_changes;
			changed_cc += tid_changes;

			ttimes[tid] += get_nano_time();
		}
		mt += get_nano_time();
		PTIP("  (3) Pointer jumping:");
		fprintf(LL_LOG, "  R%u, changed CCs: %'lu (%.2f%%), @%.2f \n", 
			mw_rank, changed_cc, 100.0 * changed_cc / g->vertices_count, (get_nano_time() + t0)/1e9);
		printf("    R%u, changed CCs: %'lu (%.2f%%), @%.2f\n", 
			mw_rank, changed_cc, 100.0 * changed_cc / g->vertices_count, (get_nano_time() + t0)/1e9);
		ul_dynamic_partitioning_reset(dp);

		{
			double changed_frac = 100.0 * changed_cc / g->vertices_count;
			int mret = MPI_Gather(&changed_frac, 1, MPI_DOUBLE, dl_array, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			assert(mret == MPI_SUCCESS);
			
			if(mw_rank == 0)
			{
				double avg = 0 , std_dev = 0;
				std_dev = get_standard_deviation_and_average(dl_array, mw_size, &avg);
				
				if(exec_info != NULL)
				{
					exec_info[4] = avg * 1e4;
					exec_info[5] = std_dev * 1e4;
				}

				fprintf(LL_LOG, "  Changed CCs Avg.: %.2f%% ± %.2f%%\n", avg, std_dev);
				printf("  Changed CCs Avg.: %.2f%% ± %.2f%%\n", avg, std_dev);
			}
		}

		if(mw_size == 1)
			goto __cenol_async_end;


	// (4) Sending changed CCs to Rank 0
		mt = - get_nano_time();
		
		unsigned long buffer_size = 1024UL * 1024 * 128;
		(mw_rank == 0) && printf("  buffer_size: %'lu\n", buffer_size);
		assert(buffer_size < (1UL<<31));
		assert(buffer_size % sizeof(unsigned long) == 0);

		if(mw_rank != 0)
		{
			unsigned int buffers_count = 0;
			MPI_Request* requests = NULL;
			unsigned int* cc_ui_changes = NULL;
			unsigned long* cc_ul_changes = NULL;
				

			if(changed_cc)
			{
				// Prefix sum of thread_changes
				{
					unsigned long sum = 0;
					for(unsigned int t = 0; t < pe->threads_count; t++)
					{
						unsigned long temp = thread_changes[t];
						thread_changes[t] = sum;
						sum += temp;
					}
					assert(sum == changed_cc);
				}

				assert(changed_cc <= g->vertices_count);
				if(cc_ui != NULL)
					cc_ui_changes = numa_alloc_interleaved(sizeof(unsigned int) * 2 * changed_cc);
				else
					cc_ul_changes = numa_alloc_interleaved(sizeof(unsigned long) * 2 * changed_cc);
				assert(cc_ui_changes != NULL || cc_ul_changes != NULL);

				// Preparing the buffer
				#pragma omp parallel num_threads(pe->threads_count)
				{
					unsigned tid = omp_get_thread_num();
					ttimes[tid] = - get_nano_time();
					
					unsigned long offset = 2 * thread_changes[tid];
					
					for(unsigned long partition = 0; partition < partitions_count; partition++)
					{
						if(p2t[partition] != tid)
							continue;

						unsigned long v_start = (g->vertices_count / partitions_count) * partition;
						unsigned long v_end = (g->vertices_count / partitions_count) * (1 + partition);
						if(partition == partitions_count - 1)
							v_end = g->vertices_count;
						
						if(cc_ui != NULL)
						{
							for(unsigned long v = v_start; v < v_end; v++)
								if(cc_ui[v] != v)
								{
									cc_ui_changes[offset++] = v;
									cc_ui_changes[offset++] = cc_ui[v];
								}
						}
						else
						{
							for(unsigned long v = v_start; v < v_end; v++)
								if(cc_ul[v] != v)
								{
									cc_ul_changes[offset++] = v;
									cc_ul_changes[offset++] = cc_ul[v];
								}
						}
					}

					ttimes[tid] += get_nano_time();
				}

				// Sending the changes
				unsigned long pairs_per_buffer = 0;
				if(cc_ui != NULL)
				{
					buffers_count = ceil((2.0 * changed_cc *  sizeof(unsigned int)) / buffer_size );
					pairs_per_buffer = buffer_size / (2 * sizeof(unsigned int));
				}
				else
				{
					buffers_count = ceil((2.0 * changed_cc * sizeof(unsigned long)) / buffer_size );
					pairs_per_buffer = buffer_size / (2 * sizeof(unsigned long));
				}
				assert(buffers_count > 0);
				// printf("  R%u, buffers_count: %u\n", mw_rank, buffers_count);
				fprintf(LL_LOG, "  R%u, buffers_count: %u\n", mw_rank, buffers_count);

				requests = calloc(sizeof(MPI_Request), buffers_count);
				assert(requests != NULL);

				for(unsigned long b = 0; b < buffers_count; b++)
				{
					char* start = NULL;
					unsigned int len = buffer_size;
					if(cc_ui != NULL)
					{
						start = (char*)(cc_ui_changes + b * 2 * pairs_per_buffer);
						if(b + 1 == buffers_count)
							len = 2 * sizeof(unsigned int) * (changed_cc - b * pairs_per_buffer);
					}
					else
					{
						start = (char*)(cc_ul_changes + b * 2 * pairs_per_buffer);
						if(b + 1 == buffers_count)
							len = 2 * sizeof(unsigned long) * (changed_cc - b * pairs_per_buffer);
					}

					assert(len != 0);
					assert(len < (1U<<31));

					requests[b] = MPI_REQUEST_NULL;
					int mret = MPI_Isend(start, len, MPI_CHAR, 0, 23, MPI_COMM_WORLD, requests + b);
					assert(mret == MPI_SUCCESS);
					assert(requests[b] != MPI_REQUEST_NULL);
				}
			}

			// Sending tag 24
			MPI_Request req_24 = MPI_REQUEST_NULL;
			{
				int mret = MPI_Isend(&buffers_count, sizeof(unsigned int), MPI_CHAR, 0, 24, MPI_COMM_WORLD, &req_24);
				assert(mret == MPI_SUCCESS);
				assert(req_24 != MPI_REQUEST_NULL);
			}

			// printf("    R%u, initiated sending %u buffers @%'.3f\n", mw_rank, buffers_count, (get_nano_time() + t0)/1e9);
			fprintf(LL_LOG, "    R%u, initiated sending %u buffers @%'.3f\n", mw_rank, buffers_count, (get_nano_time() + t0)/1e9);

			// Waiting on sent requests
			while(1)
			{
				int wait_more = 0;
			
				for(int r = 0; r < buffers_count; r++)
					if(requests[r] != MPI_REQUEST_NULL)
					{
						int flag = 0;
						int mret = MPI_Test(requests + r, &flag, MPI_STATUS_IGNORE);
						assert(mret == MPI_SUCCESS);
						if(flag)
							assert(requests[r] == MPI_REQUEST_NULL);
						else
							wait_more = 1;
					}

				if(req_24 != MPI_REQUEST_NULL)
				{
					int flag = 0;
					int mret = MPI_Test(&req_24, &flag, MPI_STATUS_IGNORE);
					assert(mret == MPI_SUCCESS);
					if(flag)
						assert(req_24 == MPI_REQUEST_NULL);
					else
						wait_more = 1;
				}
	
				if(!wait_more)
					break;

				usleep(1000);
			}

			// printf("    R%u, finished sending %u buffers @%'.3f\n", mw_rank, buffers_count, (get_nano_time() + t0)/1e9);
			fprintf(LL_LOG, "    R%u, finished sending %u buffers @%'.3f\n", mw_rank, buffers_count, (get_nano_time() + t0)/1e9);

			// Releasing mem
			if(requests != NULL)
				free(requests);
			requests = NULL;
			if(cc_ui_changes != NULL)
					numa_free(cc_ui_changes, sizeof(unsigned int) * 2 * changed_cc);
			if(cc_ul_changes != NULL)
					numa_free(cc_ul_changes, sizeof(unsigned long) * 2 * changed_cc);
			cc_ui_changes = NULL;
			cc_ul_changes = NULL;
		}
		else
		{
			// Allocating mem
				unsigned int* sent_buffs = calloc(sizeof(unsigned int), mw_size);
				unsigned int* recv_buffs = calloc(sizeof(unsigned int), mw_size);
				assert(sent_buffs != NULL && recv_buffs != NULL);
				for(int r = 0; r < mw_size; r++)
					sent_buffs[r] = -1U;

				unsigned int buffers_per_node = ceil(1.0 * g->vertices_count * sizeof(unsigned long) / buffer_size);
				unsigned int buffers_count = (mw_size - 1) * buffers_per_node ;
				unsigned long fm = get_free_mem();
				while(fm < buffers_count * buffer_size)
					buffers_count -= buffers_per_node;
				assert(buffers_count >= 1);
				printf("  R0, buffers_per_node: %u, buffers_count: %u\n", buffers_per_node, buffers_count);

				char* mem = numa_alloc_interleaved(buffers_count * buffer_size);
				assert(mem != NULL);
				#pragma omp parallel for 
				for(unsigned long c = 0; c < buffers_count * buffer_size; c+= 4096)
					mem[c] = 0;

				MPI_Request* requests = calloc(sizeof(MPI_Request), buffers_count);
				assert(requests != NULL);
				MPI_Request* req_24 = calloc(sizeof(MPI_Request), mw_size);
				assert(req_24 != NULL);

			// Initiating Irecv requests
				for(int b = 0; b < buffers_count; b++)
				{
					requests[b] = MPI_REQUEST_NULL;
					int mret = MPI_Irecv(mem + buffer_size * b, buffer_size, MPI_CHAR, MPI_ANY_SOURCE, 23, MPI_COMM_WORLD, requests + b);
					assert(mret == MPI_SUCCESS);
					assert(requests[b] != MPI_REQUEST_NULL);
				}

				for(int r = 1; r < mw_size; r++)
				{
					req_24[r] = MPI_REQUEST_NULL;
					int mret = MPI_Irecv(sent_buffs + r, 1, MPI_UNSIGNED, r, 24, MPI_COMM_WORLD, req_24 + r);
					assert(mret == MPI_SUCCESS);
					assert(req_24[r] != MPI_REQUEST_NULL);
				}

			// Receiving data and processig it
			unsigned long total_effective_changes = 0;
			unsigned long total_received = 0;
			int needs_another_iteration = 0;
			do
			{	
				// Finished nodes, tag: 24
					for(int r = 1; r < mw_size; r++)
						if(req_24[r] != MPI_REQUEST_NULL)
						{
							int flag;
							MPI_Status status;
							int mret = MPI_Test(req_24 + r, &flag, &status);
							assert(mret == MPI_SUCCESS);

							if(flag)
							{
								assert(sent_buffs[r] != -1U);
								assert(status.MPI_SOURCE == r);
								assert(req_24[r] == MPI_REQUEST_NULL);
								// printf("    R0, received tag 24 for %u buffers sent from R%u @%'.3f\n", sent_buffs[r], r, (get_nano_time() + t0)/1e9);
								fprintf(LL_LOG, "    R0, received tag 24 for %u buffers sent from R%u @%'.3f\n", sent_buffs[r], r, (get_nano_time() + t0)/1e9);
							}
							else
								assert(req_24[r] != MPI_REQUEST_NULL);
						}

				// Data, tag: 23
					for(unsigned int b = 0; b < buffers_count; b++)
					{
						// Receiving data
							assert(requests[b] != MPI_REQUEST_NULL);
							int flag;
							MPI_Status status;
							int mret = MPI_Test(requests + b, &flag, &status);
							assert(mret == MPI_SUCCESS);
							if(!flag)
								continue;
							
							assert(requests[b] == MPI_REQUEST_NULL);
							int count = 0;
							mret = MPI_Get_count(&status, MPI_CHAR, &count);
							assert(mret == MPI_SUCCESS);

							assert(status.MPI_SOURCE < mw_size);
							assert(status.MPI_SOURCE > 0);
							recv_buffs[status.MPI_SOURCE]++;
							// printf("      R0, Received a buffer of %'u B from R%u @%'.3f.\n", count, status.MPI_SOURCE, (get_nano_time() + t0)/1e9);

						// Attaching trees
							unsigned int* buff_ui = NULL;
							unsigned long* buff_ul = NULL;
							if(cc_ui != NULL)
							{
								buff_ui = (unsigned int*)(mem + buffer_size * b);
								count /= sizeof(unsigned int);
							}
							else
							{
								buff_ul = (unsigned long*)(mem + buffer_size * b);
								count /= sizeof(unsigned long);
							}
							assert(count % 2 == 0);
							total_received += count/2;
							unsigned long mt = -get_nano_time();

							unsigned long effective_changes = 0;

							#pragma omp parallel for reduction(+:effective_changes)
							for(unsigned long c = 0; c < count/2; c++)
							{
								if(cc_ui != NULL)
								{
									unsigned int x = buff_ui[ 2 * c];
									unsigned int y = buff_ui[ 2 * c + 1];
									assert(x < g->vertices_count && y < g->vertices_count);
									int round = 0;

									while(1)
									{
										while(x != cc_ui[x])
											x = cc_ui[x];

										while(y != cc_ui[y])
											y = cc_ui[y];

										if(round++ == 0 && x != y)
											effective_changes++;

										if(x == y)
											break;

										if(x < y)
										{
											if(__sync_bool_compare_and_swap(&cc_ui[y], y, x))
												break;
										}
										else
										{
											if(__sync_bool_compare_and_swap(&cc_ui[x], x, y))
												break;
										}
									}
								}
								else
								{
									unsigned long x = buff_ul[ 2 * c];
									unsigned long y = buff_ul[ 2 * c + 1];
									assert(x < g->vertices_count && y < g->vertices_count);

									int round = 0;

									while(1)
									{
										while(x != cc_ul[x])
											x = cc_ul[x];

										while(y != cc_ul[y])
											y = cc_ul[y];

										if(round++ == 0 && x != y)
											effective_changes++;

										if(x == y)
											break;

										if(x < y)
										{
											if(__sync_bool_compare_and_swap(&cc_ul[y], y, x))
												break;
										}
										else
										{
											if(__sync_bool_compare_and_swap(&cc_ul[x], x, y))
												break;
										}
									}
								}
							}
							mt += get_nano_time();
							total_effective_changes += effective_changes;
							// printf("    R0, processed in %'.3f seconds, effective: %.2f %%.\n", mt/1e9, 100.0 * effective_changes/ (count/2));

						// Initiating the Irecv
							mret = MPI_Irecv(mem + buffer_size * b, buffer_size, MPI_CHAR, MPI_ANY_SOURCE, 23, MPI_COMM_WORLD, requests + b);
							assert(mret == MPI_SUCCESS);
							assert(requests[b] != MPI_REQUEST_NULL);
					}

				// The Conclusion Procedure
					needs_another_iteration = 0;
					for(int r = 1; r < mw_size; r++)
						if(sent_buffs[r] == -1U || sent_buffs[r] > recv_buffs[r])
						{
							needs_another_iteration = 1;
							break;
						}
			}
			while(needs_another_iteration);

			// Cancelling Irecv
			for(int b = 0; b < buffers_count; b++)
			{
				assert(requests[b] != MPI_REQUEST_NULL);
				int mret = MPI_Cancel(requests + b);
				assert(mret == MPI_SUCCESS);

				mret = MPI_Request_free(requests + b);
				assert(mret == MPI_SUCCESS);
				assert(requests[b] == MPI_REQUEST_NULL);
			}

			// Some stats
			{
				double avg = 100.0 * total_effective_changes / total_received;
				double avg_per_vertex = 1.0 * total_received / g->vertices_count;
				if(exec_info != NULL)
				{
					exec_info[6] = avg * 1e4;
					exec_info[7] = avg_per_vertex * 1e4;
				}
				fprintf(LL_LOG, "  Effective transferred CCs Avg.: %.2f%%, avg. per vertex: %.2f \n", avg, avg_per_vertex);
				printf("  Effective transferred CCs Avg.: %.2f%%, avg. per vertex: %.2f \n", avg, avg_per_vertex);
			}

			// Free mem
				numa_free(mem, buffer_size * buffers_count);
				mem = NULL;

				free(sent_buffs);
				sent_buffs = NULL;
				free(recv_buffs);
				recv_buffs = NULL;

				free(requests);
				requests = NULL;
				free(req_24);
				req_24 = NULL;

			// Final pointer jumping
				#pragma omp parallel num_threads(pe->threads_count) reduction(+: changed_cc)
				{
					unsigned tid = omp_get_thread_num();
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
						
						for(unsigned long v = v_start; v < v_end; v++)
							if(cc_ui != NULL)
								while(cc_ui[cc_ui[v]] != cc_ui[v])
									cc_ui[v] = cc_ui[cc_ui[v]];
							else
								while(cc_ul[cc_ul[v]] != cc_ul[v])
									cc_ul[v] = cc_ul[cc_ul[v]];
					}
				}

				ul_dynamic_partitioning_reset(dp);
		}

		mt += get_nano_time();
		PTIP("  (4) Transfering CCs");

		free(p2t);
		p2t = NULL;
		free(thread_changes);
		thread_changes = NULL;

	// (5) Counting number of CCs
	__cenol_async_end:
		if(mw_rank == 0 && ccs_p != NULL)
		{
			unsigned long ccs = 0;
			#pragma omp parallel for num_threads(pe->threads_count) reduction(+:ccs)
			for(unsigned long v = 0; v < g->vertices_count; v++)
			{
				if(cc_ui != NULL && cc_ui[v] == v)
					ccs++;

				if(cc_ul != NULL && cc_ul[v] == v)
					ccs++;
			}
			printf("  |CCs|: %'lu\n",ccs);
			*ccs_p = ccs;
		}

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

	// Releasing memory
		
		ul_dynamic_partitioning_release(dp);
		dp = NULL;

		free(ttimes);
		ttimes = NULL;
	
	if(mw_rank == 0)
	{
		printf("\n");
		if(cc_ui != NULL)
			return cc_ui;
		else
			return cc_ul;
	}
	else
	{
		if(cc_ui != NULL)
			numa_free(cc_ui, sizeof(unsigned int) * g->vertices_count);
		else
			numa_free(cc_ul, sizeof(unsigned long) * g->vertices_count);
		
		cc_ui = NULL;
		cc_ul = NULL;

		return NULL;
	}
}

#endif
