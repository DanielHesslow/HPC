OBJS   = newton.o
CC     = gcc
CFLAGS  = -pthread -O3 -Wall -lm -march=native
.PHONY : clean

newton : $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS)

newton.o : newton.c

clean:
	rm -f $(OBJS) newton
