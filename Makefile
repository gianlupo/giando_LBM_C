CC = gcc

CFLAGS = -O3 -Wall

LDFLAGS = -lm

SRC = \
    src/main.c \
    src/init.c \
    src/force.c \
    src/lbm.c \
    src/bound.c \
    src/output.c \
    src/param.c

TARGET = run

all:
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
