CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall -lm $(DEBUG)

tomo:
	$(CC) $(LFLAGS) reconstruct.cpp -o tomo

clean:
	rm -f *.o *~ *.out tomo
