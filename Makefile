.RECIPEPREFIX = >
CC = gcc
SRC = src
OBJDIR = objdir
CFLAGS = -g -std=gnu11
LINK_LIB = -lgsl -lgslcblas -lm

#optimization level
O = -O2

VPATH = $(SRC)

_OBJS = $(patsubst $(SRC)/%.c, %.o, $(wildcard $(SRC)/*.c))
OBJS = $(patsubst %, $(OBJDIR)/%, $(_OBJS))

mcmc: $(OBJS)
> $(CC) -o $@ $^ $(CFLAGS) $(LINK_LIB)

$(OBJDIR)/main.o: main.c
> $(CC) -c -o $@ $< $(CFLAGS) $(O) $(LINK_LIB)

$(OBJDIR)/func.o: func.c
> $(CC) -c -o $@ $< $(CFLAGS) $(O) $(LINK_LIB)

$(OBJDIR)/likelihood.o: likelihood.c
> $(CC) -c -o $@ $< $(CFLAGS) $(O) $(LINK_LIB)

$(OBJDIR)/mh_mcmc.o: mh_mcmc.c
> $(CC) -c -o $@ $< $(CFLAGS) $(O) $(LINK_LIB)

$(OBJDIR)/read_val.o: read_val.c
> $(CC) -c -o $@ $< $(CFLAGS) $(O) $(LINK_LIB)

$(OBJS): | $(OBJDIR)
$(OBJDIR):
> mkdir $(OBJDIR)

.PHONY: clean
clean:
> rm -r mcmc $(OBJDIR)
