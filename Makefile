CC=gcc

all: adlibs_flip_pops adlibs_score calc_fhat pi_pops

adlibs_flip_pops:
	$(CC) src/adlibs_flip_pops.c -lm -lz -o bin/adlibs_flip_pops

adlibs_score:
	$(CC) src/adlibs_score.c -lm -lz  -o bin/adlibs_score

calc_fhat:
	$(CC) src/calc_fhat.c -lm -lz -o bin/calc_fhat

pi_pops:
	$(CC) src/pi_pops.c -lm -lz -o bin/pi_pops