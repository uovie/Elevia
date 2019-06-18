CC      := g++
INCL    := /usr/local/libint/2.5.0/include
LIBS    := /usr/local/libint/2.5.0/lib

TARGETS := Elevia
OBJ     := main.o read.o intro_basis.o core.o print.o
DEPS    := Global.h umath.h Mol_Int.h

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -I $(INCL) -L $(LIBS) -c $< -o $@

$(TARGETS): $(OBJ)
	$(CC) $^ -o $@