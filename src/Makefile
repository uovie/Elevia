CC      := g++
LIBS    := ../lib/

TARGETS := Elevia
OBJ     := main.o read.o core.o print.o

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(OBJ)

$(OBJ): %.o : %.cpp
	$(CC) -I $(LIBS) -c $< -o $@

$(TARGETS): $(OBJ)
	$(CC) $^ -o $@