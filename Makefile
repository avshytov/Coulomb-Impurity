INC = -I/usr/include/python2.7
LIB = -lm -lpython2.7
CC = gcc
LD = gcc -shared
LDFLAGS = $(LIB)
CFLAGS = -fPIC -std=gnu99 $(INC) -Wall
TARGET=_rpam.so
SRC = _rpam.c
OBJ = $(SRC:.c=.o)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<
	
_rpam.so: $(OBJ)
	$(LD)  -o $@ $(OBJ) $(LDFLAGS)
	
