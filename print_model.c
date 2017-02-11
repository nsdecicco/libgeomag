
/**
 * @file
 * @author Nicholas DeCicco <nsd.cicco@gmail.com>
 *                          <deciccon0@students.rowan.edu>
 */

#include <stdio.h>
#include <math.h>
#include "geomag.h"

void print_model(FILE *const os, BFieldModel const*const model);

int main()
{
	BFieldModel model;
	FILE *fp;

	if (!read_model(&model, "IGRF12.COF")) {
		fprintf(stderr, "Fatal: failed to open coefficient file\n");
		return 0;
	}

	if (!(fp = fopen("igrf.c", "w"))) {
		fprintf(stderr, "Fatal: failed to open file for writing\n");
		return 0;
	}

	print_model(fp, &model);

	fclose(fp);

	return 0;
}

/**
 * Generates source code from the B field model so that the model may be
 * built statically into the program.
 *
 * @param os Output file stream.
 * @param model
 */
void print_model(FILE *const os, BFieldModel const*const model)
{
	int i, j;

	fprintf(os, "#include <math.h>\n"
	            "#include \"geomag.h\"\n"
	            "#define nan NAN\n"
	            "BFieldModel model = {\n");

#define PRINT_VAR(TYPE, VAR, LEN) \
	fprintf(os, "{ "); \
	for (j = 0; j < LEN; j++) { \
		fprintf(os, TYPE "%s", VAR[j], (j == LEN-1) ? " },\n" : ", "); \
	}

#define PRINT_DOUBLE(VAR, LEN) PRINT_VAR("%lf", VAR, LEN)
#define PRINT_INT(VAR, LEN) PRINT_VAR("%d", VAR, LEN)
#define PRINT_LONG(VAR, LEN) PRINT_VAR("%ld", VAR, LEN)
#define PRINT_STR(VAR, LEN) PRINT_VAR("\"%.9s\"", VAR, LEN)

	PRINT_STR(model->name, MAXMOD)

	PRINT_DOUBLE(model->epoch, MAXMOD)
	PRINT_DOUBLE(model->yrmin, MAXMOD)
	PRINT_DOUBLE(model->yrmax, MAXMOD)

	fprintf(os, "%lf, %lf,\n", model->minyr, model->maxyr);

	PRINT_DOUBLE(model->altmin, MAXMOD)
	PRINT_DOUBLE(model->altmax, MAXMOD)

	fprintf(os, "%d,\n", model->nmodel);

	PRINT_INT(model->max1, MAXMOD);
	PRINT_INT(model->max2, MAXMOD);
	PRINT_INT(model->max3, MAXMOD);

	PRINT_LONG(model->irec_pos, MAXMOD);

#if 1
	/* Print gh1 */
	fprintf(os, "{\n");
	for (i = 0; i < MAXMOD; i++) {
		PRINT_DOUBLE(model->gh1[i], MAXCOEFF)
	}
	fprintf(os, "},\n");

	/* Print gh2 */
	fprintf(os, "{\n");
	for (i = 0; i < MAXMOD; i++) {
		PRINT_DOUBLE(model->gh2[i], MAXCOEFF)
	}
	fprintf(os, "},\n");
#endif

	fprintf(os, "};\n");
}
