OBJS   = main.o
CC     = gcc
CFLAGS  = -O2 -Wall -lm -march=native
.PHONY : clean

main : $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS)

main.o : main.c

clean:
	rm -f $(OBJS) main
