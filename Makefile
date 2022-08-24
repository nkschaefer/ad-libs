CC=gcc -std=c99 -W -Wall -Wextra
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CELLAR=$(shell brew info argp-standalone | grep Cellar | cut -d' ' -f1)
	FLAGS += -I$(CELLAR)/include/
	LDFLAGS += -L$(CELLAR)/lib/ -largp
endif

EXECUTABLES=bin/adlibs_score bin/calc_fhat bin/pi_pops

all: $(EXECUTABLES)

bin/adlibs_score: src/adlibs_score.c
	$(CC) src/adlibs_score.c $(FLAGS) $(LDFLAGS) -lm -lz  -o bin/adlibs_score

bin/calc_fhat: src/calc_fhat.c
	$(CC) src/calc_fhat.c $(FLAGS) $(LDFLAGS) -lm -lz -o bin/calc_fhat

bin/pi_pops: src/pi_pops.c
	$(CC) src/pi_pops.c $(FLAGS) $(LDFLAGS) -lm -lz -o bin/pi_pops

clean:
	rm -vf $(EXECUTABLES)
