CC = g++
CFLAGS = -Wall -O3 --std=c++17

mcsim.exe : mcsim.o
	$(CC) mcsim.o -o mcsim.exe
mcsim.o : mcsim.cpp
	$(CC) -c mcsim.cpp -o mcsim.o $(CFLAGS)
clear :
	rm *.o mcsim.exe
esegui :
	./mcsim.exe
clear_results :
	rm *.out
clear_all : clear clear_results
clear_frames :
	rm -rf xyz_frames/*.xyz
	rm -rf mgl_frames/*.mgl
	