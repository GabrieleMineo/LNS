CFLAGS = -Wall -O3 --std=c++11

main: main.o classi.o
	g++ main.o classi.o -o main

main.o: main.cpp
	g++ -c main.cpp -o main.o $(CFLAGS)

classi.o: classi.cpp classi.h
	g++ -c classi.cpp -o classi.o $(CFLAGS)

clean:
	rm -f *.o
