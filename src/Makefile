CC      := g++
LIBS    := ../lib/

TARGETS := Elevia
OBJ     := main.o read.o core.o print.o
DEPS    := Global.h misc.h umath.h Mol_Int.h SCF.h

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -I $(LIBS) -c $< -o $@

$(TARGETS): $(OBJ)
	$(CC) $^ -o $@