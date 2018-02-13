All: test 

#testUpdate
#bitwiseQuery OthelloQuery
#testSimpleInsert gentest bitwiseQuery
#testUpdate OthelloQuery OthelloKey testSimpleInsert

CXXFLAGS0= -std=c++11 -D USEBYTE -D NEWHASH -m64 -Ofast -flto -march=native -funroll-loops -D PREFETCH  -O3 -D NO_VERSION_CHECK
CXXFLAGS= -std=c++11 -D USEBYTE -D NEWHASH -D PREFETCH -g 

BUILD_NUMBER_FILE=version.txt

LDFLAGS= -lpthread
SOURCES= main.cpp othello.cc lfsrgen.cc


test: *.cc *.h 
	echo $$(($$(cat $(BUILD_NUMBER_FILE)) + 1)) > $(BUILD_NUMBER_FILE)
	VERSION= $(cat $(BUILD_NUMBER_FILE))
	g++ $(CXXFLAGS) -D VERSIONID=`cat $(BUILD_NUMBER_FILE)` $^ -o $@ 
	g++ $(CXXFLAGS0) -D VERSIONID=`cat $(BUILD_NUMBER_FILE)`	$^ -o $@O

VERSIONID: *.cc *.h


debug: CFLAGS += -D TEST -g
debug: All	
HEADERS= $(SOURCES:.cc=.h)
OBJECTS=$(SOURCES:.cc=.o)

