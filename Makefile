MPIEXEC := mpiexec
MPICC := mpicc
GCC := gcc
GXX := g++
LIB := paragrapher/lib64:$(shell echo "$(LD_LIBRARY_PATH)"|sed 's/:*$$//'):
SHELL := /bin/bash

OBJ := obj
INCLUDE_LIBS := $(addprefix -L , $(subst :, ,$(LIB))) 
INCLUDE_HEADER := $(addprefix -I , $(subst :,/../include ,$(LIB)))

OMP_WAIT_POLICY := passive
ifeq "$(wait_active)" "1" 
	OMP_WAIT_POLICY := active
endif

HT_DISABLED := 1
ifeq "$(no_ht)" "0" 
	HT_DISABLED := 0
endif

ifeq "$(energy)" "1"
	ENERGY_FLAG := -D_ENERGY_MEASUREMENT
endif

ifeq "$(papi)" "1"
	PAPI_FLAG := -D_USE_PAPI -lpapi
endif

OMP_VARS := OMP_DYNAMIC=false OMP_WAIT_POLICY=$(OMP_WAIT_POLICY)

COMPILE_TYPE := -O3 # -DNDEBUG
ifeq "$(debug)" "1"	
	COMPILE_TYPE := -g
endif	

FLAGS := $(ENERGY_FLAG) $(PAPI_FLAG) -Wfatal-errors -lm -lpthread -fopenmp -lnuma -lparagrapher -lrt # -lpfm
COMMON_C_FILES := $(shell find *.c ! -name 'alg*')
__timestamp := $(shell date +"%Y%m%d-%H%M%S")

$(OBJ)/alg%.obj: alg%.c $(COMMON_C_FILES) Makefile paragrapher
	mkdir -p $(OBJ)
	@if [ `$(GCC) -dumpversion | cut -f1 -d.` -le 8 ]; then\
		$(GCC) -dumpversion; \
		echo -e "\033[0;33mError:\033[0;37m Version 9 or newer is required for gcc.\n\n";\
		exit -1;\
	fi

	@echo -e "\n\nCompiling" $@
	$(MPICC) $(INCLUDE_HEADER) $(FLAGS) $< -std=gnu11  $(COMPILE_TYPE) -c -o $@ 
	@echo 

alg%: $(OBJ)/alg%.obj Makefile 
	@echo -e "\n\nBuilding" $@-$(__timestamp).o
	$(MPICC) $(INCLUDE_LIBS) $(OBJ)/$@.obj $(FLAGS) -o $(OBJ)/alg$*-$(__timestamp).o
	@echo ""
	@echo -e "\n\nExecuting $@-$(__timestamp).o"
	@echo -e "args: "$(args)
	PARAGRAPHER_LIB_FOLDER=paragrapher/lib64 LD_LIBRARY_PATH=$(LIB) LL_NUM_THREADS=$(threads) LL_NO_HT=$(HT_DISABLED) $(OMP_VARS) $(MPIEXEC) $(mpiexec-args) $(OBJ)/alg$*-$(__timestamp).o $(args) 
	rm $(OBJ)/alg$*-$(__timestamp).o
	
all: paragrapher Makefile
	
paragrapher: FORCE
	make -C paragrapher all

download_sample:
	@if [ ! -f data/cnr-2000.graph ]; then \
		mkdir -p data ;\
		echo -e "--------------------\n\033[1;34mDownloading cnr-2000\033[0;37m"; \
		wget -P data "http://data.law.di.unimi.it/webdata/cnr-2000/cnr-2000.graph"; \
		wget -P data "http://data.law.di.unimi.it/webdata/cnr-2000/cnr-2000.properties"; \
		echo -e "--------------------\n";\
	fi

siskincc_sample: download_sample
	mkdir -p logs
	LL_INPUT_GRAPH_PATH=data/cnr-2000 LL_INPUT_GRAPH_TYPE=PARAGRAPHER_CSX_WG_800_AP LL_OUTPUT_REPORT_PATH=logs/siskincc.txt make alg1_siskincc mpiexec-args="--report-bindings -n ${SLURM_NNODES} -npernode 1 -oversubscribe"

robincc_sample: download_sample
	mkdir -p logs
	LL_INPUT_GRAPH_PATH=data/cnr-2000 LL_INPUT_GRAPH_TYPE=PARAGRAPHER_CSX_WG_800_AP LL_OUTPUT_REPORT_PATH=logs/robincc.txt make alg2_robincc mpiexec-args="--report-bindings -n ${SLURM_NNODES} -npernode 1 -oversubscribe"


clean: unmount
	rm -f $(OBJ)/*.obj $(OBJ)/*.o ||  echo
	make -C paragrapher clean  || echo
	
unmount:
	make -C paragrapher unmount || echo
	killall pg_fuse || echo
	for f in `ls /tmp | grep pg_fuse`; do rm -rf /tmp/$$f; done || echo
	for f in `ls /tmp | grep paragrapher`; do rm -rf /tmp/$$f; done || echo
	/opt/service/bin/dropcache both

clean_shm_graphs:
	rm -f /dev/shm/ll_graph*
	
touch:
	touch *.c 

FORCE: ;

# .SUFFIXES:

.keep_exec_files: $(addprefix $(OBJ)/,$(subst .c,.obj, $(shell ls *.c)))
