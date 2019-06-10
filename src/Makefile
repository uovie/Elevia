CC      := g++
CCFLAGS := 
LDFLAGS :=

TARGETS := Elevia
MAINS   := $(addsuffix .o, $(TARGETS) )
OBJ     := read.o core.o print.o $(MAINS)
DEPS    := main.h sub.h

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJ)

$(OBJ): %.o : %.c $(DEPS)
	$(CC) -c -o $@ $< $(CCFLAGS)

$(TARGETS): % : $(filter-out $(MAINS), $(OBJ)) %.o
	$(CC) -o $@ $(LIBS) $^ $(CCFLAGS) $(LDFLAGS)