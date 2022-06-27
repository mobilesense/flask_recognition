#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include <ivfgeo/ivfgeo.h>

#define MAX_INVERTED_FILES (1 << 15)  /* 32K */

typedef struct {
  int   format;
  char *output;
  char *input[MAX_INVERTED_FILES];
  int   nb_input;
} parameters_t;

int main(int argc, char **argv)
{
  parameters_t params;

  ivfgeo_t *ivf_final = NULL;
  ivfgeo_t *ivf_tmp   = NULL;

  int character;
  int error = 0;
  int count = 0;

  /* Set defaults */
  params.format   = 0;
  params.output   = "out.ivfgeo";
  params.nb_input = 0;

  /* Parse command line */
  opterr = 0;

  while ((character = getopt(argc, argv, "f:o:i:")) >= 0) switch(character) {
  case 'f':
    params.format = atoi(optarg);
    break;
  case 'o':
    params.output = optarg;
    break;
  case 'i':
    if (params.nb_input >= MAX_INVERTED_FILES) {
      fprintf(stderr, "maximum number of input files reached, dropping %s\n", optarg);
    }
    params.input[params.nb_input++] = optarg;
    break;
  case '?':
    fprintf(stderr, "unknown option -%c\n", (char) optopt);
    error++;
  }

  if (error > 0) {
    fprintf(stderr, "%d error(s) found\n", error);
    fprintf(stdout, "Usage: %s [-f format] -o outfile [-i infile]+\n", argv[0]);
    fprintf(stdout, "    -f: ivfgeo format\n");
    fprintf(stdout, "    -o: output file name\n");
    fprintf(stdout, "    -i: input file name (may be used several times)\n");

    return error;
  }

  /* Do your job */
  if (params.nb_input > 0) {
    printf("reading %s\n", params.input[0]);
    ivf_final = ivfgeo_read(params.input[0], params.format);
  }
  for (count = 1; count < params.nb_input; count++) {
    printf("reading %s\n", params.input[count]);    
    ivf_tmp = ivfgeo_read(params.input[count], params.format);
    ivfgeo_merge(ivf_final, ivf_tmp, 0);
  }
  if (ivf_final != NULL) {
    printf("writing %s\n", params.output);
    ivfgeo_write(params.output, ivf_final, params.format);
  }

  return 0;
}
