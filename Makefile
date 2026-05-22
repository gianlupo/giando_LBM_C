CC = gcc

CFLAGS = -Wall
RELEASE_FLAGS = -O3
DEBUG_FLAGS = -g -O0 -DDEBUG

LDFLAGS = -lm

SRC = \
    src/main.c \
    src/obstacle.c \
    src/init.c \
    src/force.c \
    src/lbm.c \
    src/bound.c \
    src/output.c \
    src/param.c

TARGET = run

all:
	$(CC) $(CFLAGS) $(RELEASE_FLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

debug:
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
