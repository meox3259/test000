PROGRAM = main
CXXFLAGS = -O3 -Wall -std=c++20 -lz
OBJECTS = sequence_handler.o sam.o handle_file.o main.o utils.o

main : $(OBJECTS)
	g++ -o $(PROGRAM) $(OBJECTS) $(CXXFLAGS) 

sequence_handler.o : sequence_handler.cc sequence_handler.h
sam.o : sam.cc sam.h
handle_file.o : handle_file.cc handle_file.h
main.o : main.cc
utils.o : utils.cc utils.h

clean : 
	rm main sequence_handler.o sam.o handle_file.o main.o utils.o