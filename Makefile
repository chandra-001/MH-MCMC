.RECIPEPREFIX = >
CC = gcc
SRC = src
OBJDIR = objdir
CFLAGS = -g -std=gnu11 -Wall -Wextra
LINK_LIB = -lopenblas -lm

#optimization level
O = -O2

VPATH = $(SRC)

_OBJS = $(patsubst $(SRC)/%.c, %.o, $(wildcard $(SRC)/*.c))
OBJS = $(patsubst %, $(OBJDIR)/%, $(_OBJS))

mcmc: $(OBJS)
> $(CC) -o $@ $^ $(CFLAGS) $(LINK_LIB)

$(OBJDIR)/%.o: %.c $(SRC)/mcmc.h
> $(CC) -c -o $@ $< $(CFLAGS) $(O)

$(OBJS): | $(OBJDIR)
$(OBJDIR):
> mkdir $(OBJDIR)

.PHONY: clean
clean:
> rm -rf mcmc $(OBJDIR)
