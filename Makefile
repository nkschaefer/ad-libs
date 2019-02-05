CC=gcc -std=c99 -W -Wall -Wextra

EXECUTABLES=bin/adlibs_score bin/calc_fhat bin/pi_pops

all: $(EXECUTABLES)

bin/adlibs_score: src/adlibs_score.c
	$(CC) src/adlibs_score.c -lm -lz  -o bin/adlibs_score

bin/calc_fhat: src/calc_fhat.c
	$(CC) src/calc_fhat.c -lm -lz -o bin/calc_fhat

bin/pi_pops: src/pi_pops.c
	$(CC) src/pi_pops.c -lm -lz -o bin/pi_pops

clean:
rm -vf $(EXECUTABLES)