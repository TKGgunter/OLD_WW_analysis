#ROOT lib
ROOT_LIB = /home/gunter/build/include/ 

#Where we're going to place our shared libs
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



# compile mesh classes
TARGET = BTagCalibrationStandalone 
 
 
 


#$(TARGET).o:
#	g++ -std=c++11 -I$(PYTHON_INCLUDE) $(CFLAGS) -I$(PYBIND_INC) -fPIC -c BTagCalibrationStandalone.cpp 

$(TARGET).so: BTagCalibrationStandalone.o 
	g++ -shared   BTagCalibrationStandalone.o -L$(ROOT_LIB) -L$(PYTHON_LIB) -lpython$(PYTHON_VERSION) $(SHARED_FLAGS)  -o BTagCalibrationStandalone.so 















