
# ==========================
# circScan Makefile
# (c) 2016 Jian-Hua Yang
# yangjh7@mail.sysu.edu.cn
# ==========================
# define our object and binary directories
BIN_DIR	 = bin
CXX		 = g++
CC       = gcc
CXXFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC $(INCLUDE)
LIBS	 = -lz -lm
BAM_DIR  = thirdUtils/BamTools
UTILS_DIR = bioUtils
PROGRAM_DIR = src

all:
	cd $(BAM_DIR); make api; make
	cd $(UTILS_DIR); make
	cd $(PROGRAM_DIR)/circScan; make

clean:
	cd $(BAM_DIR); make clean_api
	cd $(UTILS_DIR); make clean
	cd $(PROGRAM_DIR)/circScan; make clean
	cd $(BIN_DIR); rm circScan
		
test:
	@cd test_data; bash run_test.sh
