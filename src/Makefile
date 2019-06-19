CC      := g++
INCS    := -I /usr/local/libint/2.5.0/include \
           -I /usr/local/libint/2.5.0/include/libint2
LIBS    := -L /usr/local/libint/2.5.0/lib

TARGETS := Elevia
OBJ     := main.o open.o intro_basis.o core.o close.o
DEPS    := Global.h umath.h Mol_Int.h

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -c $< -o $@ $(INCS) -lint2 $(LIBS)

$(TARGETS): $(OBJ)
	$(CC)  $^ -o $@ $(INCS) -lint2 $(LIBS)