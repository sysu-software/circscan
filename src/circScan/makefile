CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type
UTILITIES_DIR = ../../thirdUtils
BIN_DIR = ../../bin
BIO_DIR = ../../bioUtils
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/include/api \
           -I$(BIO_DIR)
BIO_LIBS   = -L$(UTILITIES_DIR)/BamTools/lib/ -lbamtools \
             -L$(BIO_DIR)/ -lbiotools \
             
circScan: circScan.o circScanMain.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o ${BIN_DIR}/circScan circScanMain.o circScan.o \
	$(BIO_LIBS) $(LIBS) 

circScan.o: circScan.cpp circScan.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c circScan.cpp
	
circScanMain.o: circScanMain.cpp circScan.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c circScanMain.cpp
	
clean:
	rm -f *.o
