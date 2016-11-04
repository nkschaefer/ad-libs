CC=gcc

all: bin/adlibs_flip_pops bin/adlibs_score bin/calc_fhat bin/pi_pops

bin/adlibs_flip_pops: src/adlibs_flip_pops.c
	$(CC) src/adlibs_flip_pops.c -lm -lz -o bin/adlibs_flip_pops

bin/adlibs_score: src/adlibs_score.c
	$(CC) src/adlibs_score.c -lm -lz  -o bin/adlibs_score

bin/calc_fhat: src/calc_fhat.c
	$(CC) src/calc_fhat.c -lm -lz -o bin/calc_fhat

bin/pi_pops: src/pi_pops.c
	$(CC) src/pi_pops.c -lm -lz -o bin/pi_pops
