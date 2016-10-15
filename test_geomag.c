
/**
 * @file
 * @author Nicholas DeCicco <nsd.cicco@gmail.com>
 *                          <deciccon0@students.rowan.edu>
 */

#include <stdio.h>
#include <math.h>
#include <netcdf.h>
#include "geomag.h"

#define N_LAT 40
#define N_LON 90

#define LIST_VARIABLES(VAR) \
	VAR(d) \
	VAR(i) \
	VAR(h) \
	VAR(f) \
	VAR(x) \
	VAR(y) \
	VAR(z) \
	VAR(ddot) \
	VAR(idot) \
	VAR(hdot) \
	VAR(fdot) \
	VAR(xdot) \
	VAR(ydot) \
	VAR(zdot)

#define NUM_VARIABLES 14

#define DEF_VAR(x) double x[N_LAT][N_LON];
#define DEFINE_VARIABLES LIST_VARIABLES(DEF_VAR)

#define DEF_NC_VAR(x) int x ## _var;
#define DEFINE_NC_VARS LIST_VARIABLES(DEF_NC_VAR)

#define LAT_DIM 0
#define LON_DIM 1

#define TRYNC(x) if ((status = (x)) != NC_NOERR) goto ncErr;

int main()
{
	BFieldModel model;
	BField bfield[N_LAT][N_LON];

	// Define arrays for each variable
	DEFINE_VARIABLES

	// Define
	DEFINE_NC_VARS

	int ncid;

	const double alt = 500.0;
	const double date = julday(10, 10, 2016);

	double epsilon = 10; /* We can't get the H field at the poles, so only go
	                        within 'epsilon' (degrees) of them. */
	int ii, ij;
	int dims[2];
	int status;
	double lat[N_LAT], lon[N_LON];
	int lat_var, lon_var;

	if (!read_model(&model, "IGRF12.COF")) {
		fprintf(stderr, "Fatal: failed to open coefficient file\n");
		return 0;
	}

#define COPY_VAR(x) \
	x[ii][ij] = bfield[ii][ij].x;

	for (ii = 0; ii < N_LAT; ii++) {
		for (ij = 0; ij < N_LON; ij++) {
			lat[ii] = -90.0f + epsilon + ((double) ii)/(180.0f - epsilon);
			lon[ij] = -180.0f + ((double) ij)/360.0f;
			get_field_components(bfield[ii]+ij, &model, alt, kUnitsKilometers,
			                     kCoordSysGeodetic, lat[ii], lon[ij], date,
			                     "IGRF12.COF");
			LIST_VARIABLES(COPY_VAR);
		}
	}

	if ((status = nc_create("test.nc", NC_NETCDF4, &ncid)) != NC_NOERR)
		goto ncErr_dontclose;

	TRYNC(nc_def_dim(ncid, "lat", (long) N_LAT, dims+LAT_DIM));
	TRYNC(nc_def_dim(ncid, "lon", (long) N_LON, dims+LON_DIM));

	TRYNC(nc_def_var(ncid, "lat", NC_DOUBLE, 1, dims+LAT_DIM, &lat_var));
	TRYNC(nc_def_var(ncid, "lon", NC_DOUBLE, 1, dims+LON_DIM, &lon_var));

	const char degrees_north[] = "degrees_north";
	const char degrees_east[] = "degrees_east";
	TRYNC(nc_put_att_text(ncid, lat_var, "units", sizeof(degrees_north)/sizeof(char), degrees_north));
	TRYNC(nc_put_att_text(ncid, lon_var, "units", sizeof(degrees_east)/sizeof(char), degrees_east));

	TRYNC(nc_put_var_double(ncid, lat_var, lat));
	TRYNC(nc_put_var_double(ncid, lon_var, lon));

#define CREATE_NC_VAR(x) \
	TRYNC(nc_def_var(ncid, #x, NC_DOUBLE, 2, dims, & x ## _var));

	LIST_VARIABLES(CREATE_NC_VAR);

#define FILL_VAR(x) \
	TRYNC(nc_put_var_double(ncid, x ## _var, (double*) x));

	LIST_VARIABLES(FILL_VAR);

	nc_close(ncid);

	return 0;

ncErr:
	nc_close(ncid);

ncErr_dontclose:
	fprintf(stderr, "NetCDF error occured: \"%s\"", nc_strerror(status));
	return 1;
}
