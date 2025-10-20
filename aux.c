#ifndef __AUX_C
#define __AUX_C

#ifndef _DEFAULT_SOURCE // _BSD_SOURCE
	#define _DEFAULT_SOURCE // _BSD_SOURCE
#endif
#define _GNU_SOURCE

#include <unistd.h>
#include <numa.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <asm/unistd.h>
#include <fcntl.h>
#include <sys/syscall.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <linux/hw_breakpoint.h>
#include <errno.h>
#include <stdint.h>
#include <linux/limits.h>
#include <inttypes.h>
#include <linux/futex.h>
#include <emmintrin.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <libgen.h>
#include <pthread.h>
#include <cpuid.h>
#include <sched.h>

#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

long run_command(char* in_cmd, char* output, unsigned int output_size);
unsigned int get_available_cpus_count();

char* LL_INPUT_GRAPH_PATH = "data/test_csr.txt";
char* LL_INPUT_GRAPH_TYPE = "text";
char LL_INPUT_GRAPH_BASE_NAME [256] = {0};
unsigned int LL_INPUT_GRAPH_BATCH_ORDER = 0;
unsigned int LL_INPUT_GRAPH_IS_SYMMETRIC = 0;
unsigned int LL_STORE_INPUT_GRAPH_IN_SHM = 0;
char* LL_OUTPUT_REPORT_PATH = NULL;
char* LL_GRAPH_RA_BIN_FOLDER = "data/RA_bin_arrays";
unsigned int LL_NUM_THREADS = 0;
unsigned int LL_NO_HT = 0;
FILE* LL_LOG = NULL;

struct par_env;
struct par_env* pe = NULL;

int mw_rank; // MPI World rank
int mw_size; // MPI World size
unsigned int mt_upperbound = 0; // MPI tag upperbound
unsigned int* nodes_cores = NULL;
unsigned int total_cores = 0;
unsigned long total_mem = 0;
unsigned int lower_ranks_cores = 0;
unsigned long* ul_array = NULL;
double* dl_array = NULL;

unsigned long get_nano_time();

void read_env_vars()
{
	char* temp = getenv("LL_INPUT_GRAPH_PATH");
	if(temp != NULL)
		LL_INPUT_GRAPH_PATH = strdup(temp);

	temp = getenv("LL_INPUT_GRAPH_TYPE");
	if(temp != NULL)
		LL_INPUT_GRAPH_TYPE = strdup(temp);
	
	{
		int len = strlen(LL_INPUT_GRAPH_PATH);
		temp = calloc(len + 64, 1);
		assert(temp != NULL);

		if(strcmp(LL_INPUT_GRAPH_TYPE, "binary") == 0)
		{
			char* bn = strrchr(LL_INPUT_GRAPH_PATH, '/');
			if(bn == NULL)
				bn = LL_INPUT_GRAPH_PATH;
			else
				bn++;
			int bn_len = min(256, strlen(bn));
			sprintf(LL_INPUT_GRAPH_BASE_NAME, "%.*s", bn_len, bn);
		}
		else if(strcmp(LL_INPUT_GRAPH_TYPE, "text") == 0)
		{
			sprintf(temp, "%s", LL_INPUT_GRAPH_PATH);
			char* bn = basename(temp);
			int bn_len = strlen(bn);
			bn_len = min(bn_len - strlen(".txt"), 256);
			sprintf(LL_INPUT_GRAPH_BASE_NAME, "%.*s", bn_len, bn);
		}
		else if( strcmp(LL_INPUT_GRAPH_TYPE, "PARAGRAPHER_CSX_WG_400_AP") == 0 ||
			strcmp(LL_INPUT_GRAPH_TYPE, "PARAGRAPHER_CSX_WG_404_AP") == 0 ||
		 	strcmp(LL_INPUT_GRAPH_TYPE, "PARAGRAPHER_CSX_WG_800_AP") == 0
		)
		{
			sprintf(temp, "%s.properties", LL_INPUT_GRAPH_PATH);
			char* bn = basename(temp);
			int bn_len = strlen(bn);
			bn_len = min(bn_len - strlen(".properties"), 256);
			sprintf(LL_INPUT_GRAPH_BASE_NAME, "%.*s", bn_len, bn);
		}

		free(temp);
		temp = NULL;
	}

	temp = getenv("LL_INPUT_GRAPH_BATCH_ORDER");
	if(temp != NULL)
		LL_INPUT_GRAPH_BATCH_ORDER = atoi(temp);

	temp = getenv("LL_INPUT_GRAPH_IS_SYMMETRIC");
	if(temp != NULL && strcmp(temp, "1") == 0)
		LL_INPUT_GRAPH_IS_SYMMETRIC = 1;

	temp = getenv("LL_STORE_INPUT_GRAPH_IN_SHM");
	if(temp != NULL && strcmp(temp, "1") == 0)
		LL_STORE_INPUT_GRAPH_IN_SHM = 1;

	temp = getenv("LL_OUTPUT_REPORT_PATH");
	if(temp != NULL && strlen(temp) > 0)
		LL_OUTPUT_REPORT_PATH = strdup(temp);

	temp = getenv("LL_GRAPH_RA_BIN_FOLDER");
	if(temp != NULL)
		LL_GRAPH_RA_BIN_FOLDER = strdup(temp);

	temp = getenv("LL_NUM_THREADS");
	if(temp != NULL)
		LL_NUM_THREADS = atoi(temp);

	temp = getenv("LL_NO_HT");
	if(temp != NULL && strcmp(temp, "1") == 0)
		LL_NO_HT = 1;

	{
		char hostname[128];
		gethostname(hostname, 127);
		for(int i = 0; i < min(strlen(hostname), 128); i++)
			if(hostname[i] == '.')
			{
				hostname[i] = 0;
				break;
			}

		char output_file[PATH_MAX];
		sprintf(output_file, "/dev/null");

		if(LL_OUTPUT_REPORT_PATH != NULL && access(LL_OUTPUT_REPORT_PATH, F_OK) == 0)
		{	
			char* dn = dirname(strdup(LL_OUTPUT_REPORT_PATH));

			sprintf(output_file, "%s/%u-%s-%u-%s-%u.log", 
				dn, 
				LL_INPUT_GRAPH_BATCH_ORDER, 
				LL_INPUT_GRAPH_BASE_NAME,
				(unsigned int)(get_nano_time()/60e9),
				hostname,
				syscall(SYS_gettid)
			);
		}

		LL_LOG = fopen(output_file, "w");
		assert(LL_LOG != NULL);

		printf("Started on %s, cpus: %u, log file: %s\n", hostname, get_available_cpus_count(), output_file);
		fprintf(LL_LOG, "Started on %s, cpus: %u, log file: %s\n", hostname, numa_num_configured_cpus(), output_file);
	}

	fprintf(LL_LOG, "\n\033[2;32mEnv Vars:\033[0;37m\n");
	fprintf(LL_LOG, "    LL_INPUT_GRAPH_PATH:         %s\n", LL_INPUT_GRAPH_PATH);
	fprintf(LL_LOG, "    LL_INPUT_GRAPH_BASE_NAME:    %s\n", LL_INPUT_GRAPH_BASE_NAME);
	fprintf(LL_LOG, "    LL_INPUT_GRAPH_TYPE:         %s\n", LL_INPUT_GRAPH_TYPE);
	fprintf(LL_LOG, "    LL_INPUT_GRAPH_BATCH_ORDER:  %u\n", LL_INPUT_GRAPH_BATCH_ORDER);
	fprintf(LL_LOG, "    LL_INPUT_GRAPH_IS_SYMMETRIC: %u\n", LL_INPUT_GRAPH_IS_SYMMETRIC);
	fprintf(LL_LOG, "    LL_STORE_INPUT_GRAPH_IN_SHM: %u\n", LL_STORE_INPUT_GRAPH_IN_SHM);
	fprintf(LL_LOG, "    LL_OUTPUT_REPORT_PATH:       %s\n", LL_OUTPUT_REPORT_PATH);
	fprintf(LL_LOG, "    LL_GRAPH_RA_BIN_FOLDER:      %s\n", LL_GRAPH_RA_BIN_FOLDER);
	fprintf(LL_LOG, "    LL_NUM_THREADS:              %u\n", LL_NUM_THREADS);
	fprintf(LL_LOG, "    LL_NO_HT:                    %u\n", LL_NO_HT);
	fflush(LL_LOG);
	
	return;
}

double get_standard_deviation(double* data, int count)
{
	double sum = 0;
	for (int i = 0; i < count; i++) 
		sum += data[i];
	double avg = sum / count;

	sum = 0;
	for (int i = 0; i < count; i++) 
		sum += pow(data[i] - avg, 2);

	return sqrt(sum / count);
}

double get_standard_deviation_and_average(double* data, int count, double* avg)
{
	assert(avg != NULL);
	*avg = 0;
	double sum = 0;
	for (int i = 0; i < count; i++) 
		sum += data[i];
	*avg = sum / count;

	sum = 0;
	for (int i = 0; i < count; i++) 
		sum += pow(data[i] - *avg, 2);

	return sqrt(sum / count);
}

unsigned long get_nano_time()
{
	struct timespec ts;
	timespec_get(&ts,TIME_UTC);
	return ts.tv_sec*1e9+ts.tv_nsec;
}

void quick_sort_uint(unsigned int* bucket, unsigned long start, unsigned long end)
{
	if(end - start < 16)
	{
		do{
			unsigned long new_end = start;

			for(unsigned long i=start+1; i<=end; i++)	
				if(bucket[i-1] > bucket[i])
				{
					unsigned int t=bucket[i];
					bucket[i]=bucket[i-1];
					bucket[i-1]=t;

					new_end = i;
				}

			end = new_end;
		}
		while(end > start);

		return;
	}

	// moving the middle index to the end to suffle for semi-sorted arrays
	unsigned int t=bucket[(end+start)/2];
	bucket[(end+start)/2] = bucket[end];
	bucket[end] = t;

	unsigned long pivot_index=end;
	unsigned long front_index=start;

	unsigned int bp;
	unsigned int p=bucket[pivot_index];
	unsigned int f;

	while(pivot_index > front_index)
	{
		//printf("%2d %2d ** \t",pivot_index, front_index);
		if(bucket[front_index] >= p)
		{
			bp=bucket[pivot_index-1];
			f=bucket[front_index];

			bucket[pivot_index]=f;
			bucket[front_index]=bp;

			pivot_index--;
			
		}
		else
			front_index++;
	}

	bucket[pivot_index]=p;

	if(pivot_index > start + 1)
		quick_sort_uint(bucket, start, pivot_index-1);

	if(pivot_index < end - 1)
		quick_sort_uint(bucket, pivot_index+1, end);

	return;
}

void quick_sort_ushort(unsigned short* bucket, unsigned long start, unsigned long end)
{
	if(end - start < 16)
	{
		do{
			unsigned long new_end = start;

			for(unsigned long i=start+1; i<=end; i++)	
				if(bucket[i-1] > bucket[i])
				{
					unsigned short t=bucket[i];
					bucket[i]=bucket[i-1];
					bucket[i-1]=t;

					new_end = i;
				}

			end = new_end;
		}
		while(end > start);

		return;
	}

	// moving the middle index to the end to suffle for semi-sorted arrays
	unsigned short t=bucket[(end+start)/2];
	bucket[(end+start)/2] = bucket[end];
	bucket[end] = t;

	unsigned long pivot_index=end;
	unsigned long front_index=start;

	unsigned short bp;
	unsigned short p=bucket[pivot_index];
	unsigned short f;

	while(pivot_index > front_index)
	{
		//printf("%2d %2d ** \t",pivot_index, front_index);
		if(bucket[front_index] >= p)
		{
			bp=bucket[pivot_index-1];
			f=bucket[front_index];

			bucket[pivot_index]=f;
			bucket[front_index]=bp;

			pivot_index--;
			
		}
		else
			front_index++;
	}

	bucket[pivot_index]=p;

	if(pivot_index > start + 1)
		quick_sort_ushort(bucket, start, pivot_index-1);

	if(pivot_index < end - 1)
		quick_sort_ushort(bucket, pivot_index+1, end);

	return;
}

char* get_date_time(char* in)
{
	time_t t0 = time(NULL);
	struct tm * t = localtime(&t0);
	sprintf(in, "%04d/%02d/%02d %02d:%02d:%02d", 1900 + t->tm_year, 1 + t->tm_mon, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
	return in;
}

unsigned long get_file_size(char* input_file_name)
{
	int fd = open(input_file_name,O_RDONLY); 
	assert(fd > 0);

	struct stat st;
	assert(0 == fstat(fd, &st));
	close(fd);
	fd = -1;
	
	return st.st_size;
}

int get_file_contents(char* file_name, char* buff, int buff_size)
{
	int fd=open(file_name, O_RDONLY);
	if(fd<0)
	{
		fprintf(LL_LOG, "Can't open the %s : %d - %s\n",file_name, errno, strerror(errno));
		return -1;
	}

	int count = read(fd, buff, buff_size);
	if(count<buff_size)
		buff[count]=0;

	close(fd);

	return count;
}

long run_command(char* in_cmd, char* output, unsigned int output_size)
{
	assert(in_cmd != NULL);

	int in_cmd_size = strlen(in_cmd);
	char* cmd = malloc(in_cmd_size + 2 * 1024);
	assert(cmd != NULL);
	
	char hostname[128];
	gethostname(hostname, 127);

	char* res_file = cmd + in_cmd_size + 1024;
	sprintf(res_file, "_temp_res_%s_%u_%lu_%lu.txt", hostname, mw_rank, syscall(SYS_gettid),get_nano_time());	
	sprintf(cmd, "%s 1>%s 2>&1", in_cmd, res_file);

	int ret = system(cmd);
	long ret2 = 0;
	if(output != NULL && output_size != 0)
		ret2 = get_file_contents(res_file, output, output_size);
	ret2 = (ret2 << 32) + ret;

	if(access(res_file, F_OK) == 0)
	{
		ret = unlink(res_file);
		assert(ret == 0);
	}
	
	memset(cmd, 0, in_cmd_size + 2 * 64);
	free(cmd);
	cmd = NULL;

	return ret2;
}

// Binary Search for `val` in [`start`, `end`) of `vals` as a sorted array  
// returns -1UL if it cannot find
// otherwise, returns index of the element that is equal to val
unsigned long uint_binary_search(unsigned int* vals, unsigned long in_start, unsigned long in_end, unsigned int val)
{
	unsigned long end = in_end;
	unsigned long start = in_start;
	assert(start <= end);
	if(start == end)
		return -1UL;

	if(vals[start] > val) 
		return -1UL;

	if(vals[end - 1] < val)
		return -1UL; 

	unsigned long m = (start + end)/2;
	unsigned int r = 0;
	unsigned int r_max = 1 + log2(end + 1 - start);
	while(1)
	{
		unsigned int m_val = vals[m];
		if( m_val == val )
			return m;

		if(m + 1 < in_end)
		{
			if(vals[m+1] == val)
				return m+1;

			if(m_val < val && vals[m + 1] > val)
				return -1UL;
		}

		if(m_val > val)
			end = m;
		if(m_val < val)
			start = m;
		m = (start + end)/2;

		assert(r++ <= r_max);
	}

	assert("Don't reach here.");
	return -1UL;
}

void flush_os_cache()
{
	char* ts = calloc(1024 + get_nprocs()/4, 1);
	assert(ts != NULL);
	sprintf(ts, "taskset 0x");
	for(int i = 0; i < get_nprocs() / 4 + 1; i++)
		sprintf(ts + strlen(ts), "F");
	
	sprintf(ts + strlen(ts), " ./flushcache.sh");
	run_command(ts, NULL, 0);

	free(ts);
	ts = NULL;

	return;
}

void* create_shm(char* shm_file_name, unsigned long length)
{
	assert(shm_file_name != NULL);
	assert(length > 0);

 	int shm_fd = shm_open(shm_file_name, O_RDWR|O_CREAT, 0644);
 	if(shm_fd == -1)
	{
		fprintf(LL_LOG, "create_shm(), error in shm_open() %d, %s .\033[0;37m \n", errno, strerror(errno));
		assert(shm_fd != -1);
	}	

	int ret = ftruncate(shm_fd, length);
	if(ret != 0)
	{
		fprintf(LL_LOG, "create_shm(), error in ftruncate() %d, %s .\033[0;37m \n", errno, strerror(errno));
		ret = shm_unlink(shm_file_name);
		assert(ret == 0);
		return NULL;
	}

	void* mem = mmap(NULL, length, PROT_READ | PROT_WRITE, MAP_SHARED , shm_fd, 0);
	if(mem == MAP_FAILED)
	{
		fprintf(LL_LOG, "create_shm(), error in mmap() %d, %s .\033[0;37m \n", errno, strerror(errno));
		ret = shm_unlink(shm_file_name);
		assert(ret == 0);
		return NULL;
	}

	close(shm_fd);
	shm_fd = -1;

	return mem;
}

char* ul2s(unsigned long in, char* out)
{
	assert(out != NULL);

	if(in > 1e21)
		sprintf(out, "%'.1fZ", in/1e21);
	else if(in > 1e18)
		sprintf(out, "%.1fE", in/1e18);
	else if(in > 1e15)
		sprintf(out, "%.1fP", in/1e15);
	else if(in > 1e12)
		sprintf(out, "%.1fT", in/1e12);
	else if(in > 1e9)
		sprintf(out, "%.1fG", in/1e9);
	else if(in > 1e6)
		sprintf(out, "%.1fM", in/1e6);
	else if(in > 1e3)
		sprintf(out, "%.1fk", in/1e3);
	else
		sprintf(out, "%lu", in);

	return out;
}

char* ul2s_bin(unsigned long in, char* out)
{
	assert(out != NULL);

	if(in >= (1UL<<60))
		sprintf(out, "%luE", in/(1UL<<60));
	else if(in >= (1UL<<50))
		sprintf(out, "%luP", in/(1UL<<50));
	else if(in >= (1UL<<40))
		sprintf(out, "%luT", in/(1UL<<40));
	else if(in >= (1UL<<30))
		sprintf(out, "%luG", in/(1UL<<30));
	else if(in >= (1UL<<20))
		sprintf(out, "%luM", in/(1UL<<20));
	else if(in >= (1UL<<10))
		sprintf(out, "%luk", in/(1UL<<10));
	else
		sprintf(out, "%lu", in);

	return out;
}


// SPLITMIXX64 (https://prng.di.unimi.it/splitmix64.c) to fill the seeds using the seed x
void rand_initialize_splitmix64(unsigned long s[4], unsigned long x)
{	
	for(int si = 0; si < 4; si++)
	{
		unsigned long z = (x += 0x9e3779b97f4a7c15);
		z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
		z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
		s[si] = z ^ (z >> 31);
	}

	return;
}

/* 
	XOSHIRO, 
	
	https://prng.di.unimi.it/
	https://prng.di.unimi.it/xoshiro256plusplus.c
	https://prng.di.unimi.it/xoroshiro64star.c

	@article
	{
		10.1145/3460772,
		author = {Blackman, David and Vigna, Sebastiano},
		title = {Scrambled Linear Pseudorandom Number Generators},
		year = {2021},
		issue_date = {December 2021},
		publisher = {Association for Computing Machinery},
		address = {New York, NY, USA},
		volume = {47},
		number = {4},
		issn = {0098-3500},
		url = {https://doi.org/10.1145/3460772},
		doi = {10.1145/3460772},
		journal = {ACM Trans. Math. Softw.},
		month = {sep},
		articleno = {36},
		numpages = {32}
	}
*/
unsigned long rand_xoshiro256(unsigned long s[4]) 
{
	unsigned long rand_val = s[0] + (((s[0] + s[3]) << 23) | ((s[0] + s[3]) >> (64 - 23)));

	unsigned long t = s[1] << 17;
	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];
	s[2] ^= t;
	s[3] = (s[3] << 45) | (s[3] >> (64 - 45));

	return rand_val;
}

unsigned int  rand_xoroshiro64star(unsigned int s[2]) {
	const unsigned int  s0 = s[0];
	unsigned int  s1 = s[1];
	const unsigned int  result = s0 * 0x9E3779BB;

	s1 ^= s0;
	
	// s[0] = rotl(s0, 26) ^ s1 ^ (s1 << 9); // a, b
	s[0] = ((s0 << 26) | (s0 >> (32 - 26))) ^ s1 ^ (s1 << 9); // a, b
	
	// s[1] = rotl(s1, 13); // c
	s[1] = ((s1 << 13) | (s1 >> (32 - 13))); // c

	return result;
}

unsigned int get_available_cpus_count()
{
	cpu_set_t set;
	unsigned int count = 0;
	
	int ret = pthread_getaffinity_np(pthread_self(), sizeof (set), &set);
	assert(ret == 0);

	for (int i = 0; i < CPU_SETSIZE; i++)
		if (CPU_ISSET (i, &set))
			count++;

	return count;
}

unsigned long get_net_usage()
{
    char input[1024] = {0};
    
    // https://enterprise-support.nvidia.com/s/article/understanding-mlx5-linux-counters-and-status-parameters
    int ret = run_command("cat /sys/class/infiniband/*/ports/1/counters/port_*_data | xargs | sed -r 's/ /+/g'|bc", input, 1024);
    unsigned long sum = 4 * atol(input);
    ret = run_command("cat /sys/class/net/*/statistics/*_bytes | xargs | sed -r 's/ /+/g' | bc", input, 1024);
    sum += atol(input);

    return sum;    
}

#endif