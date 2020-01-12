CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.cpp TSM.o random.o
	mpic++ main.cpp TSM.o random.o -o main.exe
random.o : random.h random.cpp
	$(CC) -c random.cpp -o random.o $(CFLAGS)	
TSM.o : TSM.cpp TSM.h 
	$(CC) -c TSM.cpp -o TSM.o $(CFLAGS)
clean :
	rm *.o es10.exe
