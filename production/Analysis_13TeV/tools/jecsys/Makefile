#ROOT lib
ROOT_LIB = /usr/local/Cellar/root/5.34.32/lib/root/

#Where we're going to plave our shared libs
SHARED_LIB = $(shell pwd)/lib/


#PYTHON stuff
PYTHON_VERSION = 2.7
#PYTHON_INCLUDE = /usr/local/Cellar/python/2.7.10_2/Frameworks/Python.framework/Versions/2.7/include/python2.7/
#PYTHON_LIB = /usr/local/Cellar/python/2.7.10_2/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config/
PYTHON_INCLUDE = /home/gunter/anaconda2/include/python2.7/
PYTHON_LIB = /home/gunter/anaconda2/lib/python2.7/config/

#Location of pybind header files
#PYBIND_INC = /Users/thothgunter/C++_wrapping/pydbind11/pybind11/include/ 
PYBIND_INC = /home/gunter/pybind/pybind11-master/include/ 

#Adding python and pybind include libs
CFLAGS += -I$(PYTHON_INCLUDE) -I$(PYBIND_INC) -Iinterface -I/home/gunter/build/include/

#pre-compilation and compliation object defs and shared lib flags
SHARED_FLAGS = -lpthread -ldl -lutil `root-config --libs` -lMinuit -lSpectrum -lEve -lGeom

INCLUDE_DIR = /home/gunter/WW_analysis/production/Analysis_13TeV/tools/jecsys 

OBJS  		= $(patsubst CondFormats/JetMETObjects/src/%.cc,lib/%.o,$(wildcard CondFormats/JetMETObjects/src/*.cc))
EXECS 		= $(patsubst bin/%.cc,%,$(wildcard bin/*.cc))
EXEOBJS  	= $(patsubst bin/%.cc,lib/%.o,$(wildcard bin/*.cc))



# compile mesh classes
TARGET = JECUncertainty
all: $(OBJS) $(TARGET).o ErrorTypes.o $(TARGET).so
 

 
lib/%.o : CondFormats/JetMETObjects/src/%.cc
	g++ -std=c++11 -fPIC -Wall $(CFLAGS) $(INCLUDE) -I$(INCLUDE_DIR) -c $< -o $@

ErrorTypes.o : ErrorTypes.cpp
	g++ -std=c++11 -fPIC -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@

$(TARGET).o: $(TARGET).cpp 
	g++ -std=c++11 -fPIC -I$(PYTHON_INCLUDE) $(CFLAGS) -I$(PYBIND_INC)  -c $(TARGET).cpp

$(TARGET).so: $(TARGET).o $(OBJS)
	g++ -shared -fPIC $(TARGET).o $(OBJS) ErrorTypes.o -L$(PYTHON_LIB) -lpython$(PYTHON_VERSION) $(SHARED_FLAGS)  -o $(TARGET).so
#
#lib/%.o : bin/%.cc
#	g++ -Wall $(CFLAGS) $(SHARED) $(INCLUDE) -c $< -o $@
#
#
#% : $(OBJS) lib/%.o
#  $(LD) $(LIBS) $(OBJS) $(LINKDEFOBJ) lib/$@.o -o $@














