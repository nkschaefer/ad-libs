/*
Given output from adlibs_hmm.py (via stdin), flips ancestral states, so that
AA becomes BB and BB becomes AA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    // Buffer to store lines
    char line[100];
    while(fgets(line, sizeof(line), stdin)){
        // State name should be the last 2 letters on the line.
        int linelen = strlen(line);
        if (line[linelen-3] == 'A' && line[linelen-2] == 'A'){
            line[linelen-3] = 'B';
            line[linelen-2] = 'B';
        }
        else if (line[linelen-3] == 'B' && line[linelen-2] == 'B'){
            line[linelen-3] = 'A';
            line[linelen-2] = 'A';
        }
        fprintf(stdout, "%s", line);
    }

    return 0;
}
