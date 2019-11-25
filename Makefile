include mach.inc

EXE = cloudFV
SOURCES = driver.cpp FileIO.cpp params.cpp YAKL.cpp Exchange.cpp Initializer.cpp Parser.cpp Tendencies.cpp TimeIntegrator.cpp\
           TendenciesThetaConsADER.cpp TendenciesThetaConsSD.cpp TendenciesThetaPrimADER.cpp TendenciesThetaPrimSD.cpp
OBJECTS = $(SOURCES:.cpp=.o)

default: main

main: $(EXE)

$(EXE): $(OBJECTS)
	$(LINK) $(OBJECTS) $(LIB) -o $(EXE) $(LDFLAGS)

%.o:%.cpp *.h
	$(CXX) $(CXXFLAGS) -I./cub -c $< -o $(notdir $@)

clean:
	rm -f *.o $(EXE)

realclean: clean

distclean: clean


