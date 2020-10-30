CC=g++
CFLAGS= -Wall -g 

gustavson: gustavson.cpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm *.o gustavson 
