#ifndef __PARTITIONING_C
#define __PARTITIONING_C

int dll_800_parallel_edge_partitioning(struct dll_800_graph* g, unsigned long* borders, unsigned long partitions)
{
	assert(partitions > 0 && g->vertices_count > 0);

	borders[0] = 0;
	borders[partitions] = g->vertices_count;

	if(g->vertices_count <= 4 * partitions)
	{
		unsigned long remained_vertices = g->vertices_count;
		unsigned long remained_partitions = partitions;
		for(unsigned long i = 1; i < partitions; i++)
		{
			if(remained_vertices > 0)
			{
				unsigned long q = max(1, remained_vertices/remained_partitions);
				borders[i] = borders[i-1] + q;
				remained_vertices -= q;
			}
			else
				borders[i] = g->vertices_count;

			//printf("%d %d %d\n",i,borders[i], remained_vertices);
			remained_partitions--;
		}
	}
	else
	{
		unsigned long edges_per_thread = (g->edges_count + g->vertices_count) / partitions;
		
		#pragma omp parallel for
		for(unsigned long t = 1; t < partitions; t++)
		{
			unsigned long start = 0;
			unsigned long end = g->vertices_count;
			unsigned long target = t * edges_per_thread;
			unsigned long m = (start + end)/2;
			while(1)
			{
				unsigned long m_val = g->whole_offsets_list[m] + m;
				if( m_val == target )
					break;

				unsigned long b_val = g->whole_offsets_list[m - 1] + m - 1;
				if( b_val < target && m_val > target )
					break;

				if(m_val > target)
					end = m;
				if(m_val < target)
					start = m;

				unsigned long new_m = (start + end)/2;
				if(new_m == m)
					break;
				m = new_m;
			}
			assert( m <= g->vertices_count );
			
			borders[t] = m;
		}	

		unsigned long last_m = 0;
		for(unsigned long t = 1; t < partitions; t++)
		{
			if(borders[t] <= last_m)
				borders[t] = last_m + 1;
			if(borders[t] > g->vertices_count)
				borders[t] = g->vertices_count;
			
			last_m = borders[t];
		}
	}

	//verify_ll_800_serial_edge_partitioning(g, borders, partitions);

	return 0;
}

struct ul_dynamic_partitioning
{
	struct par_env* pe;

	unsigned long partitions_count;
	unsigned long partitions_remained;
	unsigned long threads_count;

	unsigned long* threads_partitions_start;
	unsigned long* threads_partitions_end;
	unsigned long* threads_partitions_remained;

	unsigned long* threads_last_steal_offset;
	
	unsigned char* partitions_status;
};

struct ul_dynamic_partitioning* ul_dynamic_partitioning_initialize(struct par_env* pe, unsigned long partitions_count)
{
	assert(partitions_count > 0 && pe != NULL);
	
	struct ul_dynamic_partitioning* dp = calloc(sizeof(struct ul_dynamic_partitioning), 1);
	assert(dp != NULL);
	dp->partitions_count = partitions_count;
	dp->pe = pe;
	dp->threads_count = pe->threads_count;

	dp->threads_partitions_start = calloc(pe->threads_count ,sizeof(unsigned long));
	dp->threads_partitions_end = calloc(pe->threads_count ,sizeof(unsigned long));
	dp->threads_partitions_remained = calloc(pe->threads_count,sizeof(unsigned long));
	assert(dp->threads_partitions_start != NULL && dp->threads_partitions_end != NULL && dp->threads_partitions_remained != NULL);

	dp->threads_last_steal_offset = calloc(pe->threads_count ,sizeof(unsigned long));
	assert(dp->threads_last_steal_offset != NULL);

	dp->partitions_status = calloc(dp->partitions_count, sizeof(unsigned char));
	assert(dp->partitions_status != NULL);

	dp->partitions_remained = partitions_count;

	// Partitioning for nodes and threads
	unsigned long remained_threads = pe->threads_count;
	unsigned long remained_partitions = dp->partitions_count;
	unsigned long offset = 0;

	for(long t=0; t<pe->threads_count; t++)
	{	
		dp->threads_partitions_start[t] = offset;
		unsigned long quota = remained_partitions/remained_threads;
		offset += quota;
		dp->threads_partitions_end[t] = offset;
		dp->threads_partitions_remained[t] = quota;

		remained_threads--;
		remained_partitions -= quota;

		dp->threads_last_steal_offset[t] = 0;

		// printf("\t\tPartitions for thread %u: %'u - %'u\n",t, dp->threads_partitions_start[t], dp->threads_partitions_end[t]);
	}

	return dp;
}


void ul_dynamic_partitioning_release(struct ul_dynamic_partitioning* dp)
{
	assert(dp != NULL);
	dp->pe = NULL;

	free(dp->threads_partitions_start);
	dp->threads_partitions_start = NULL;

	free(dp->threads_partitions_end);
	dp->threads_partitions_end = NULL;

	free(dp->threads_partitions_remained);
	dp->threads_partitions_remained = NULL;

	free(dp->partitions_status);
	dp->partitions_status = NULL;

	free(dp->threads_last_steal_offset);
	dp->threads_last_steal_offset = NULL;

	free(dp);

	return;
}

void ul_dynamic_partitioning_reset(struct ul_dynamic_partitioning* dp)
{
	assert(dp->partitions_remained == 0);
	dp->partitions_remained =  dp->partitions_count;
	
	for(unsigned long i=0; i<dp->partitions_count; i++)
	{
		assert(dp->partitions_status[i] == 1);
		dp->partitions_status[i] = 0;
	}

	for(unsigned long t=0; t<dp->threads_count; t++)
	{
		assert(dp->threads_partitions_remained[t] == 0);
		dp->threads_partitions_remained[t] =  dp->threads_partitions_end[t] - dp->threads_partitions_start[t];
		dp->threads_last_steal_offset[t] = 0;
	}

	return;
}


unsigned long ul_dynamic_partitioning_get_next_partition(struct ul_dynamic_partitioning* dp, unsigned long thread_id, unsigned long prev_partition)
{
	if(dp->partitions_remained == 0)
		return -1UL;

	if(prev_partition >= dp->partitions_count)
		prev_partition = dp->threads_partitions_start[thread_id];	
	
	for(unsigned long pi = 0; pi < dp->partitions_count; pi++)
	{
		unsigned long p = pi + prev_partition;
		if(p >= dp->partitions_count)
			p -= dp->partitions_count;

		if(dp->partitions_status[p] == 1)
			continue;

		if(!__sync_bool_compare_and_swap(&dp->partitions_status[p], 0, 1))
			continue;

		__sync_add_and_fetch(&dp->partitions_remained, -1UL);

		int s = 0;
		for(unsigned long t = 0; t < dp->threads_count; t++)
			if(p < dp->threads_partitions_end[t])
			{
				s=1;
				__sync_add_and_fetch(&dp->threads_partitions_remained[t], -1UL);
				break;
			}
		assert(s == 1);
		return p;
	}

	return -1UL;
}
#endif 