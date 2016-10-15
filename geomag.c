/****************************************************************************/
/*                                                                          */
/*     NGDC's Geomagnetic Field Modeling software for the IGRF and WMM      */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Disclaimer: This program has undergone limited testing. It is        */
/*     being distributed unofficially. The National Geophysical Data        */
/*     Center does not guarantee it's correctness.                          */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Version 7.0:                                                         */
/*     - input file format changed to                                       */
/*            -- accept new DGRF2005 coeffs with 0.01 nT precision          */
/*            -- make sure all values are separated by blanks               */
/*            -- swapped n and m: first is degree, second is order          */
/*     - new my_isnan function improves portablility                        */
/*     - corrected feet to km conversion factor                             */
/*     - fixed date conversion errors for yyyy,mm,dd format                 */
/*     - fixed lon/lat conversion errors for deg,min,sec format             */
/*     - simplified leap year identification                                */
/*     - changed comment: units of ddot and idot are arc-min/yr             */
/*     - added note that this program computes the secular variation as     */
/*            the 1-year difference, rather than the instantaneous change,  */
/*            which can be slightly different                               */
/*     - clarified that height is above ellipsoid, not above mean sea level */
/*            although the difference is negligible for magnetics           */
/*     - changed main(argv,argc) to usual definition main(argc,argv)        */
/*     - corrected rounding of angles close to 60 minutes                   */
/*     Thanks to all who provided bug reports and suggested fixes           */
/*                                                                          */
/*                                          Stefan Maus Jan-25-2010         */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Version 6.1:                                                         */
/*     Included option to read coordinates from a file and output the       */
/*     results to a new file, repeating the input and adding columns        */
/*     for the output                                                       */
/*                                          Stefan Maus Jan-31-2008         */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Version 6.0:                                                         */
/*     Bug fixes for the interpolation between models. Also added warnings  */
/*     for declination at low H and corrected behaviour at geogr. poles.    */
/*     Placed print-out commands into separate routines to facilitate       */
/*     fine-tuning of the tables                                            */
/*                                          Stefan Maus Aug-24-2004         */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      This program calculates the geomagnetic field values from           */
/*      a spherical harmonic model.  Inputs required by the user are:       */
/*      a spherical harmonic model data file, coordinate preference,        */
/*      altitude, date/range-step, latitude, and longitude.                 */
/*                                                                          */
/*         Spherical Harmonic                                               */
/*         Model Data File       :  Name of the data file containing the    */
/*                                  spherical harmonic coefficients of      */
/*                                  the chosen model.  The model and path   */
/*                                  must be less than PATH chars.           */
/*                                                                          */
/*         Coordinate Preference :  Geodetic (WGS84 latitude and altitude   */
/*                                  above ellipsoid (WGS84),                */
/*                                  or geocentric (spherical, altitude      */
/*                                  measured from the center of the Earth). */
/*                                                                          */
/*         Altitude              :  Altitude above ellipsoid (WGS84). The   */
/*                                  program asks for altitude above mean    */
/*                                  sea level, because the altitude above   */
/*                                  ellipsoid is not known to most users.   */
/*                                  The resulting error is very small and   */
/*                                  negligible for most practical purposes. */
/*                                  If geocentric coordinate preference is  */
/*                                  used, then the altitude must be in the  */
/*                                  range of 6370.20 km - 6971.20 km as     */
/*                                  measured from the center of the earth.  */
/*                                  Enter altitude in kilometers, meters,   */
/*                                  or feet                                 */
/*                                                                          */
/*         Date                  :  Date, in decimal years, for which to    */
/*                                  calculate the values of the magnetic    */
/*                                  field.  The date must be within the     */
/*                                  limits of the model chosen.             */
/*                                                                          */
/*         Latitude              :  Entered in decimal degrees in the       */
/*                                  form xxx.xxx.  Positive for northern    */
/*                                  hemisphere, negative for the southern   */
/*                                  hemisphere.                             */
/*                                                                          */
/*         Longitude             :  Entered in decimal degrees in the       */
/*                                  form xxx.xxx.  Positive for eastern     */
/*                                  hemisphere, negative for the western    */
/*                                  hemisphere.                             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      Subroutines called :  julday,getshc,interpsh,                       */
/*                            extrapsh,shval3,dihf                          */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>            
#include <string.h>
#include <ctype.h>
#include <math.h> 

int my_isnan(double d)
{
	return (d != d);              /* IEEE: only NaN is not equal to itself */
}

#define NaN log(-1.0)
#define FT2KM (1.0/0.0003048)
#define PI 3.141592654
#define RAD2DEG (180.0/PI)

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#define IEXT 0
#define FALSE 0
#define TRUE 1                  /* constants */
#define RECL 81

#define MAXINBUFF RECL+14

/** Max size of in buffer **/

#define MAXREAD MAXINBUFF-2
/** Max to read 2 less than total size (just to be safe) **/

#define MAXMOD 30
/** Max number of models in a file **/

#define PATH MAXREAD
/** Max path and filename length **/

#define EXT_COEFF1 (double)0
#define EXT_COEFF2 (double)0
#define EXT_COEFF3 (double)0

#define MAXDEG 13
#define MAXCOEFF (MAXDEG*(MAXDEG+2)+1) /* index starts with 1!, (from old Fortran?) */
double d=0,f=0,h=0,i=0;
double dtemp,ftemp,htemp,itemp;
double x=0,y=0,z=0;
double xtemp,ytemp,ztemp;
double gh1[MAXCOEFF]; /** Schmidt quasi-normal internal spherical harmonic coeff. */
double gh2[MAXCOEFF]; /** Schmidt quasi-normal internal spherical harmonic coeff. */
double gha[MAXCOEFF]; /** Coefficients of resulting model. */
double ghb[MAXCOEFF]; /** Coefficients of rate of change model. */

/*
 *
 *                             Program Geomag
 *
 ***************************************************************************
 *
 *      This program, originally written in FORTRAN, was developed using
 *      subroutines written by
 *      A. Zunde
 *      USGS, MS 964, Box 25046 Federal Center, Denver, Co.  80225
 *      and
 *      S.R.C. Malin & D.R. Barraclough
 *      Institute of Geological Sciences, United Kingdom.
 *
 *      Translated
 *      into C by    : Craig H. Shaffer
 *                     29 July, 1988
 *
 *      Rewritten by : David Owens
 *                     For Susan McLean
 *
 *      Maintained by: Adam Woods
 *      Contact      : geomag.models@noaa.gov
 *                     National Geophysical Data Center
 *                     World Data Center-A for Solid Earth Geophysics
 *                     NOAA, E/GC1, 325 Broadway,
 *                     Boulder, CO  80303
 *
 *
 ***************************************************************************
 *
 *      Some variables used in this program
 *
 *    Name         Type                    Usage
 * ------------------------------------------------------------------------
 *
 *   alt        Scalar Double          altitude above WGS84 Ellipsoid
 *
 *   epoch      Double array of MAXMOD epoch of model.
 *
 *   ext        Scalar Double          Three 1st-degree external coeff.
 *
 *   latitude   Scalar Double          Latitude.
 *
 *   longitude  Scalar Double          Longitude.
 *
 *   gh1        Double array           Schmidt quasi-normal internal
 *                                     spherical harmonic coeff.
 *
 *   gh2        Double array           Schmidt quasi-normal internal
 *                                     spherical harmonic coeff.
 *
 *   gha        Double array           Coefficients of resulting model.
 *
 *   ghb        Double array           Coefficients of rate of change model.
 *
 *   irec_pos   Integer array of MAXMOD Record counter for header
 *
 */

double julday();
int interpsh();
int extrapsh();
int shval3();
int dihf();
int getshc();

typedef enum {
	kUnitsKilometers = 1,
	kUnitsMeters     = 2,
	kUnitsFeet       = 3
} Units;

typedef enum {
	kCoordSysGeodetic, /* WGS-84 */
	kCoordSysGeocentric
} CoordinateSystem;

typedef struct {
	double d;    /** Declination of the field from the
	                 geographic north (deg). */
	double i;    /** Inclination of the field (deg). */
	double h;    /** */
	double x;    /** */
	double y;    /** */
	double z;    /** */
	double ddot; /** Annual rate of change of declination. (arc-min/yr) */
	double fdot; /** */
	double hdot; /** */
	double idot; /** Annual rate of change of inclination. (arc-min/yr). */
	double xdot; /** */
	double ydot; /** */
	double zdot; /** */
} BField;

/**
 * @param alt Altitude, in units specified by altUnits.
 * @param altUnits
 * @param latitude North latitude, in degrees.
 * @param longitude East (?) longitude, in degrees.
 * @param sdate Start date.
 */
int get_field_components(const double alt,
                         const Units altUnits,
                         const CoordinateSystem coordSys,
                         const double latitude,
                         const double longitude,
                         const double sdate,
                         )
{
	int warn_H, warn_H_strong, warn_P;
	
	int modelI;       /* Which model (Index) */
	int nmodel;       /* Number of models in file */
	int max1[MAXMOD]; /* Main field coefficient. */
	int max2[MAXMOD]; /* Secular variation coefficient. */
	int max3[MAXMOD]; /* Acceleration coefficient. */
	int nmax;
	long  irec_pos[MAXMOD];
	

	char mdfile[PATH];
	char inbuff[MAXINBUFF];
	char model[MAXMOD][9];

	double epoch[MAXMOD];
	double yrmin[MAXMOD]; /* Min year of model. */
	double yrmax[MAXMOD]; /* Max year of model. */
	double minyr;         /* Min year of all models. */
	double maxyr;         /* Max year of all models. */
	double altmin[MAXMOD]; /* Minimum height of each model. */
	double altmax[MAXMOD]; /* Maximum height of each model. */
	double minAlt;         /* Minimum height of selected model. */
	double maxAlt;         /* Maximum height of selected model. */
	double alt=-999999;
	double warn_H_val, warn_H_strong_val;

	/* Initializations. */

	inbuff[MAXREAD+1]   = '\0';  /* Just to protect mem. */
	inbuff[MAXINBUFF-1] = '\0';  /* Just to protect mem. */

	/*  Obtain the desired model file and read the data  */

	warn_H = 0;
	warn_H_val = 99999.0;
	warn_H_strong = 0;
	warn_H_strong_val = 99999.0;
	warn_P = 0;

	if (sdate < minyr || sdate > maxyr+1) {
		return 0;
	} 

	/* Pick model */
	for (modelI = 0; modelI < nmodel; modelI++) {
		if (sdate < yrmax[modelI]) break;
	}
	/* if beyond end of last model use last model */
	if (modelI == nmodel) modelI--;

	/* Get altitude min and max for selected model. */
	minalt = altmin[modelI];
	maxalt = altmax[modelI];

	if (coordSys == kCoordSysGeocentric) {
		/* Add Earth radius to ranges. */
		minalt += 6371.2;
		maxalt += 6371.2;
	}

	if (coordSys == kCoordSysGeocentric && altUnits != kUnitsKilometers) {
		fprintf(stderr, "get_field_components: altitude must be specified in "
		                "kilometers with geocentric coordinate system\n");
		return 0;
	}

	/* Do unit conversions if neccessary */
	switch (altUnits) {
		case kUnitsMeters:
			minAlt *= 1000.0;
			maxAlt *= 1000.0;
			break;

		case kUnitsFeet:
			minAlt *= FT2KM;
			maxAlt *= FT2KM;
			break;
	}

	if (alt < minAlt || alt > maxAlt) {
		return 0;
	}

	/* Convert altitude to km */
	switch (altUnits) {
		case kUnitsMeters:
			alt *= 0.001;
			break;

		case kUnitsFeet:
			alt /= FT2KM;
			break;
	}

	if (latitude < -90 || latitude > 90 ||
	    longitude < -180 || longitude > 180)
	{
		return 0;
	}

	/* This will compute everything needed for 1 point in time. */

	if (max2[modelI] == 0) {
		getshc(mdfile, 1, irec_pos[modelI], max1[modelI], 1);
		getshc(mdfile, 1, irec_pos[modelI+1], max1[modelI+1], 2);
		nmax = interpsh(sdate, yrmin[modelI], max1[modelI],
		                       yrmin[modelI+1], max1[modelI+1], 3);
		nmax = interpsh(sdate+1, yrmin[modelI] , max1[modelI],
		                         yrmin[modelI+1], max1[modelI+1],4);
	} else {
		getshc(mdfile, 1, irec_pos[modelI], max1[modelI], 1);
		getshc(mdfile, 0, irec_pos[modelI], max2[modelI], 2);
		nmax = extrapsh(sdate, epoch[modelI], max1[modelI], max2[modelI], 3);
		nmax = extrapsh(sdate+1, epoch[modelI], max1[modelI], max2[modelI], 4);
	}

	/* Do the first calculations */
	shval3(coordSys, latitude, longitude, alt, nmax, 3,
	       IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3);
	dihf(3);
	shval3(coordSys, latitude, longitude, alt, nmax, 4,
	       IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3);
	dihf(4);

	bfield->ddot = ((dtemp - d)*RAD2DEG);
	if (bfield->ddot > 180.0) ddot -= 360.0;
	if (bfield->ddot <= -180.0) ddot += 360.0;
	bfield->ddot *= 60.0;

	bfield->idot = ((itemp - i)*RAD2DEG)*60;
	d = d*(RAD2DEG);   i = i*(RAD2DEG);
	bfield->hdot = htemp - h;
	bfield->xdot = xtemp - x;
	bfield->ydot = ytemp - y;
	bfield->zdot = ztemp - z;
	bfield->fdot = ftemp - f;

	/* Deal with geographic and magnetic poles */
	if (h < 100.0) /* at magnetic poles */
	{
		d = NaN;
		ddot = NaN;
		/* while rest is ok */
	}

	if (h < 1000.0) {
		warn_H = 0;
		warn_H_strong = 1;
		if (h < warn_H_strong_val) {
			warn_H_strong_val = h;
		}
	} else if (h < 5000.0 && !warn_H_strong)  {
		warn_H = 1;
		if (h<warn_H_val) {
			warn_H_val = h;
		}
	}

	/* at geographic poles */
	if (90.0 - fabs(latitude) <= 0.001) {
		bfield->x = NaN;
		bfield->y = NaN;
		bfield->d = NaN;
		bfield->xdot = NaN;
		bfield->ydot = NaN;
		bfield->ddot = NaN;
		warn_P = 1;
		warn_H = 0;
		warn_H_strong = 0;
		/* while rest is ok */
	}

	return 1;
}

/**
 *
 *
 * @param mdfile Model file name.                   
 */
int read_model(const char mdfile[])
{
	int fileline = 0; /* First line will be 1 */
	FILE *stream = fopen(mdfile, "rt");

	rewind(stream);

	modelI = -1;                             /* First model will be 0 */
	while (fgets(inbuff, MAXREAD, stream))
	{
		fileline++;                           /* On new line */

		if (strlen(inbuff) != RECL)       /* IF incorrect record size */
		{
			printf("Corrupt record in file %s on line %d.\n", mdfile, fileline);
			fclose(stream);
			return 0;
		}

#if 0
		/* old statement Dec 1999 */
		if (!strncmp(inbuff, "    ", 4)){ /* If 1st 4 chars are spaces */
#else
		/* New statement Dec 1999 changed by wmd  required by year 2000 models */
		if (!strncmp(inbuff, "   ", 3)) /* If 1st 3 chars are spaces */
#endif
		{
			/* New model. */
			modelI++;

			/* If too many headers */
			if (modelI > MAXMOD)
			{
				printf("Too many models in file %s on line %d.", mdfile, fileline);
				fclose(stream);
				return 0;
			}

			irec_pos[modelI]=ftell(stream);
			/* Get fields from buffer into individual vars.  */
			sscanf(inbuff, "%s%lg%d%d%d%lg%lg%lg%lg", model[modelI],
			       &epoch[modelI], &max1[modelI], &max2[modelI], &max3[modelI],
			       &yrmin[modelI], &yrmax[modelI], &altmin[modelI],
			       &altmax[modelI]);

			/* Compute date range for all models */
			if (modelI == 0) { /* If first model */
				minyr = yrmin[0];
				maxyr = yrmax[0];
			} else {
				if (yrmin[modelI] < minyr) {
					minyr=yrmin[modelI];
				}
				if (yrmax[modelI] > maxyr){
					maxyr = yrmax[modelI];
				}
			}
		}
	}

	nmodel = modelI + 1;
	fclose(stream);

	/* if date specified in command line then warn if past end of validity */

	if ((sdate > maxyr) && (sdate < maxyr + 1)) {
		printf("\nWarning: The date %4.2f is out of range,\n"
		       "         but still within one year of model expiration date.\n"
		       "         An updated model file is available before 1.1.%4.0f\n",
		       sdate, maxyr);
	}

	return 1;
}

/**
 * Computes the decimal day of year from month, day, year.
 * @author Daniel Bergstrom
 *
 * References:
 *
 * 1. Nachum Dershowitz and Edward M. Reingold, Calendrical Calculations,
 *    Cambridge University Press, 3rd edition, ISBN 978-0-521-88540-9.
 *
 * 2. Claus Tøndering, Frequently Asked Questions about Calendars,
 *    Version 2.9, http://www.tondering.dk/claus/calendar.html
 *
 * @param month
 * @param day
 * @param year
 */
double julday(const int month, const int day, const int year)
{
	int days[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };

	int leap_year = (((year % 4) == 0) &&
	                (((year % 100) != 0) || ((year % 400) == 0)));

	double day_in_year = (days[month - 1] + day + (month > 2 ? leap_year : 0));

	return ((double)year + (day_in_year / (365.0 + leap_year)));
}

/**
 * Reads spherical harmonic coefficients from the specified model into an
 * array.
 *
 * @param stream     Logical unit number
 * @param iflag      Flag for SV equal to ) or not equal to 0
 *                   for designated read statements
 * @param strec      Starting record number to read from model
 * @param nmax_of_gh Maximum degree and order of model
 *
 * @param gh1 or 2   Schmidt quasi-normal internal spherical
 *                   harmonic coefficients
 *
 *     FORTRAN
 *           Bill Flanagan
 *           NOAA CORPS, DESDIS, NGDC, 325 Broadway, Boulder CO.  80301
 *
 *     C
 *           C. H. Shaffer
 *           Lockheed Missiles and Space Company, Sunnyvale CA
 *           August 15, 1988
 *
 */
int getshc(char file[PATH], int iflag, long int strec, int nmax_of_gh, int gh)
{
	char  inbuff[MAXINBUFF];
	char irat[9];
	int ii,m,n,mm,nn;
	int ios;
	int line_num;
	double g,hh;
	double trash;

	if (!(stream = fopen(file, "rt"))) {
		fprintf(stderr, "\nError on opening file %s", file);
		return ios;
	}

	ii = 0;
	ios = 0;
	fseek(stream,strec,SEEK_SET);
	for (nn = 1; nn <= nmax_of_gh; nn++)
	{
		for (mm = 0; mm <= nn; mm++)
		{
			if (iflag == 1)
			{
				fgets(inbuff, MAXREAD, stream);
				sscanf(inbuff, "%d%d%lg%lg%lg%lg%s%d",
				       &n, &m, &g, &hh, &trash, &trash, irat, &line_num);
			}
			else
			{
				fgets(inbuff, MAXREAD, stream);
				sscanf(inbuff, "%d%d%lg%lg%lg%lg%s%d",
				       &n, &m, &trash, &trash, &g, &hh, irat, &line_num);
			}
			if ((nn != n) || (mm != m))
			{
				ios = -2;
				fclose(stream);
				return ios;
			}
			ii = ii + 1;
			switch(gh)
			{
				case 1:
					gh1[ii] = g;
					break;
				case 2:
					gh2[ii] = g;
					break;
				default:
					fprintf(stderr, "\nError in subroutine getshc");
					break;
			}
			if (m != 0)
			{
				ii = ii+ 1;
				switch (gh)
				{
					case 1:
						gh1[ii] = hh;
						break;
					case 2:
						gh2[ii] = hh;
						break;
					default:
						fprintf(stderr, "\nError in subroutine getshc");
						break;
				}
			}
		}
	}

	fclose(stream);
	return ios;
}


/**
 * Extrapolates linearly a spherical harmonic model with a rate-of-change
 * model.
 *
 * @param date   date of resulting model (in decimal year)
 * @param dte1   date of base model
 * @param nmax1  maximum degree and order of base model
 * @param gh1    Schmidt quasi-normal internal spherical harmonic coefficients
 *               of base model
 * @param nmax2  maximum degree and order of rate-of-change model
 * @param gh2    Schmidt quasi-normal internal spherical harmonic coefficients
 *               of rate-of-change model
 *
 * @param gha    Schmidt quasi-normal internal spherical harmonic coefficients
 * @param ghb    Schmidt quasi-normal internal spherical harmonic coefficients
 * @param nmax   maximum degree and order of resulting model
 *
 *     FORTRAN
 *           A. Zunde
 *           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225
 *
 *     C
 *           C. H. Shaffer
 *           Lockheed Missiles and Space Company, Sunnyvale CA
 *           August 16, 1988
 *
 */
int extrapsh(double date, double dte1, int nmax1, int nmax2, int gh)
{
	int    nmax;
	int    k, l;
	int    ii;
	double factor;

	factor = date - dte1;
	if (nmax1 == nmax2)
	{
		k =  nmax1 * (nmax1 + 2);
		nmax = nmax1;
	}
	else
	{
		if (nmax1 > nmax2)
		{
			k = nmax2 * (nmax2 + 2);
			l = nmax1 * (nmax1 + 2);
			switch(gh)
			{
				case 3:
					for (ii = k + 1; ii <= l; ii++)
					{
						gha[ii] = gh1[ii];
					}
					break;
				case 4:
					for ( ii = k + 1; ii <= l; ++ii)
					{
						ghb[ii] = gh1[ii];
					}
					break;
				default:
					printf("\nError in subroutine extrapsh");
					break;
			}
			nmax = nmax1;
		}
		else
		{
			k = nmax1 * (nmax1 + 2);
			l = nmax2 * (nmax2 + 2);
			switch(gh)
			{
				case 3:
					for (ii = k + 1; ii <= l; ii++)
					{
						gha[ii] = factor * gh2[ii];
					}
					break;
				case 4:
					for (ii = k + 1; ii <= l; ii++)
					{
						ghb[ii] = factor * gh2[ii];
					}
					break;
				default:
					printf("\nError in subroutine extrapsh");
					break;
			}
			nmax = nmax2;
		}
	}
	switch(gh)
	{
		case 3:
			for (ii = 1; ii <= k; ii++)
			{
				gha[ii] = gh1[ii] + factor * gh2[ii];
			}
			break;
		case 4:
			for (ii = 1; ii <= k; ++ii)
			{
				ghb[ii] = gh1[ii] + factor * gh2[ii];
			}
			break;
		default:
			printf("\nError in subroutine extrapsh");
			break;
	}
	return(nmax);
}

/**
 * Interpolates linearly, in time, between two spherical harmonic models.
 *
 * @param date  date of resulting model (in decimal year)
 * @param dte1  date of earlier model
 * @param nmax1 maximum degree and order of earlier model
 * @param gh1   Schmidt quasi-normal internal spherical harmonic coefficients
 *              of earlier model
 * @param dte2  date of later model
 * @param nmax2 maximum degree and order of later model
 * @param gh2   Schmidt quasi-normal internal spherical harmonic coefficients
 *              of internal model
 *
 * @param gha Coefficients of resulting model
 * @param ghb Coefficients of resulting model
 * @param nmax Maximum degree and order of resulting model
 *
 * @author Original Fortran code by
 *         A. Zunde
 *         USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225
 *
 * @author Conversion to C by
 *         C. H. Shaffer
 *         Lockheed Missiles and Space Company, Sunnyvale CA
 *         August 17, 1988
 */
int interpsh(double date, double dte1, int nmax1, double dte2, int nmax2,
             int gh)
{
	int   nmax;
	int   k, l;
	int   ii;
	double factor;

	factor = (date - dte1) / (dte2 - dte1);
	if (nmax1 == nmax2)
	{
		k =  nmax1 * (nmax1 + 2);
		nmax = nmax1;
	}
	else
	{
		if (nmax1 > nmax2)
		{
			k = nmax2 * (nmax2 + 2);
			l = nmax1 * (nmax1 + 2);
			switch(gh)
			{
				case 3:
					for (ii = k + 1; ii <= l; ii++)
					{
						gha[ii] = gh1[ii] + factor * (-gh1[ii]);
					}
					break;
				case 4:
					for (ii = k + 1; ii <= l; ii++)
					{
						ghb[ii] = gh1[ii] + factor * (-gh1[ii]);
					}
					break;
				default:
					printf("\nError in subroutine extrapsh");
					break;
			}
			nmax = nmax1;
		}
		else
		{
			k = nmax1 * (nmax1 + 2);
			l = nmax2 * (nmax2 + 2);
			switch (gh)
			{
				case 3:
					for (ii = k + 1; ii <= l; ii++)
					{
						gha[ii] = factor * gh2[ii];
					}
					break;
				case 4:
					for (ii = k + 1; ii <= l; ii++)
					{
						ghb[ii] = factor * gh2[ii];
					}
					break;
				default:
					fprintf(stderr, "\nError in subroutine extrapsh");
					break;
			}
			nmax = nmax2;
		}
	}
	switch(gh)
	{
		case 3:
			for (ii = 1; ii <= k; ii++)
			{
				gha[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
			}
			break;
	    case 4:
			for (ii = 1; ii <= k; ii++)
			{
				ghb[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
			}
			break;
		default:
			printf("\nError in subroutine extrapsh");
			break;
	}
	return(nmax);
}

/**
 * Calculates field components from spherical harmonic (sh) models.
 *
 * based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin,
 * report no. 71/1, institute of geological sciences, U.K.
 *
 *     FORTRAN
 *           Norman W. Peddie
 *           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225
 *
 *     C
 *           C. H. Shaffer
 *           Lockheed Missiles and Space Company, Sunnyvale CA
 *           August 17, 1988
 *
 * @param coordSys  indicates coordinate system used
 * @param latitude  north latitude, in degrees
 * @param longitude east longitude, in degrees
 * @param elev      WGS84 altitude above ellipsoid (coordSys is
 *                  kCoordSysGeodetic), or radial distance from earth's center
 *                  (coordSys is kCoordSysGeocentric)
 * @param a2,b2     squares of semi-major and semi-minor axes of the reference
 *                  spheroid used for transforming between geodetic and
 *                  geocentric coordinates or components
 * @param nmax      maximum degree and order of coefficients
 * @param iext      external coefficients flag (=0 if none)
 * @param ext1,2,3  the three 1st-degree external coefficients
 *                  (not used if iext = 0)
 *
 * @param x northward component
 * @param y eastward component
 * @param z vertically-downward component
 */
int shval3(const CoordinateSystem coordSys, double flat, double flon,
           double elev, int lnmax, int gh,
           int iext, double ext1, double ext2, double ext3)
{
	double earths_radius = 6371.2;
	double dtr = 0.01745329;
	double slat;
	double clat;
	double ratio;
	double aa, bb, cc, dd;
	double sd;
	double cd;
	double r;
	double a2;
	double b2;
	double rr;
	double fm,fn;
	double sl[14];
	double cl[14];
	double p[119];
	double q[119];
	int ii,j,k,l,m,n;
	int npq;
	int ios;
	double argument;
	double power;
	a2 = 40680631.59; /* Square of the semi-major axes of the WGS84 reference
	                     sphereoid used for transforming between geodetic and
	                     geocentric coordinates. */
	b2 = 40408299.98; /* Square of the semi-minor axes of the WGS84 reference
	                     sphereoid used for transforming between geodetic and
	                     geocentric coordinates. */
	ios = 0;
	r = elev;
	argument = flat * dtr;
	slat = sin( argument );
	if ((90.0 - flat) < 0.001)
	{
		aa = 89.999;            /*  300 ft. from North pole  */
	}
	else
	{
		if ((90.0 + flat) < 0.001)
		{
			aa = -89.999;        /*  300 ft. from South pole  */
		}
		else
		{
			aa = flat;
		}
	}
	argument = aa * dtr;
	clat = cos(argument);
	argument = flon * dtr;
	sl[1] = sin(argument);
	cl[1] = cos(argument);
	switch (gh)
	{
		case 3:
			x = 0;
			y = 0;
			z = 0;
			break;
		case 4:
			xtemp = 0;
			ytemp = 0;
			ztemp = 0;
			break;
		default:
			fprintf(stderr, "\nError in subroutine shval3");
			break;
	}

	sd = 0.0;
	cd = 1.0;
	l = 1;
	n = 0;
	m = 1;
	npq = (nmax * (nmax + 3)) / 2;
	if (coordSys == kCoordSysGeodetic)
	{
		aa = a2 * clat * clat;
		bb = b2 * slat * slat;
		cc = aa + bb;
		argument = cc;
		dd = sqrt( argument );
		argument = elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc;
		r = sqrt( argument );
		cd = (elev + dd) / r;
		sd = (a2 - b2) / dd * slat * clat / r;
		aa = slat;
		slat = slat * cd - clat * sd;
		clat = clat * cd + aa * sd;
	}
	ratio = earths_radius / r;
	argument = 3.0;
	aa = sqrt( argument );
	p[1] = 2.0 * slat;
	p[2] = 2.0 * clat;
	p[3] = 4.5 * slat * slat - 1.5;
	p[4] = 3.0 * aa * clat * slat;
	q[1] = -clat;
	q[2] = slat;
	q[3] = -3.0 * clat * slat;
	q[4] = aa * (slat * slat - clat * clat);
	for ( k = 1; k <= npq; ++k)
	{
		if (n < m)
		{
			m = 0;
			n = n + 1;
			argument = ratio;
			power =  n + 2;
			rr = pow(argument,power);
			fn = n;
		}
		fm = m;
		if (k >= 5)
		{
			if (m == n)
			{
				argument = (1.0 - 0.5/fm);
				aa = sqrt( argument );
				j = k - n - 1;
				p[k] = (1.0 + 1.0/fm) * aa * clat * p[j];
				q[k] = aa * (clat * q[j] + slat/fm * p[j]);
				sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
				cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
			}
			else
			{
				argument = fn*fn - fm*fm;
				aa = sqrt( argument );
				argument = ((fn - 1.0)*(fn-1.0)) - (fm * fm);
				bb = sqrt( argument )/aa;
				cc = (2.0 * fn - 1.0)/aa;
				ii = k - n;
				j = k - 2 * n + 1;
				p[k] = (fn + 1.0) * (cc * slat/fn * p[ii] - bb/(fn - 1.0) * p[j]);
				q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j];
			}
		}
		switch (gh)
		{
			case 3:
				aa = rr * gha[l];
				break;
			case 4:
				aa = rr * ghb[l];
				break;
			default:
				printf(stderr, "\nError in subroutine shval3");
			break;
		}
		if (m == 0)
		{
			switch(gh)
			{
				case 3:
					x = x + aa * q[k];
					z = z - aa * p[k];
					break;
				case 4:
					xtemp = xtemp + aa * q[k];
					ztemp = ztemp - aa * p[k];
					break;
				default:
					printf("\nError in subroutine shval3");
					break;
			}
			l = l + 1;
		}
		else
		{
			switch(gh)
			{
				case 3:
					bb = rr * gha[l+1];
					cc = aa * cl[m] + bb * sl[m];
					x = x + cc * q[k];
					z = z - cc * p[k];
					if (clat > 0)
					{
						y = y + (aa * sl[m] - bb * cl[m]) *
						fm * p[k]/((fn + 1.0) * clat);
					}
					else
					{
						y = y + (aa * sl[m] - bb * cl[m]) * q[k] * slat;
					}
					l = l + 2;
					break;
				case 4:
					bb = rr * ghb[l+1];
					cc = aa * cl[m] + bb * sl[m];
					xtemp = xtemp + cc * q[k];
					ztemp = ztemp - cc * p[k];
					if (clat > 0)
					{
						ytemp += (aa * sl[m] - bb * cl[m]) *
						         fm * p[k]/((fn + 1.0) * clat);
					}
					else
					{
						ytemp += (aa * sl[m] - bb * cl[m]) * q[k] * slat;
					}
					l = l + 2;
					break;
				default:
					fprintf(stderr, "\nError in subroutine shval3");
					break;
			}
		}
		m = m + 1;
	}
	if (iext != 0)
	{
		aa = ext2 * cl[1] + ext3 * sl[1];
		switch (gh)
		{
			case 3:
				x = x - ext1 * clat + aa * slat;
				y = y + ext2 * sl[1] - ext3 * cl[1];
				z = z + ext1 * slat + aa * clat;
				break;
			case 4:
				xtemp = xtemp - ext1 * clat + aa * slat;
				ytemp = ytemp + ext2 * sl[1] - ext3 * cl[1];
				ztemp = ztemp + ext1 * slat + aa * clat;
				break;
			default:
				printf("\nError in subroutine shval3");
				break;
		}
	}
	switch(gh)
	{
		case 3:
			aa = x;
			x = x * cd + z * sd;
			z = z * cd - aa * sd;
			break;
		case 4:
			aa = xtemp;
			xtemp = xtemp * cd + ztemp * sd;
			ztemp = ztemp * cd - aa * sd;
			break;
		default:
			fprintf(stderr, "\nError in subroutine shval3");
			break;
	}
	return ios;
}

/**
 * Computes the geomagnetic d, i, h, and f from x, y, and z.
 *
 * @param gh
 * @param x northward component
 * @param y eastward component
 * @param z vertically-downward component
 * @param d declination
 * @param i inclination
 * @param h horizontal intensity
 * @param f total intensity
 *
 * @author Fortran written by A. Zunde
 *         USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225
 * @author C implementation written by
 *         C. H. Shaffer
 *         Lockheed Missiles and Space Company, Sunnyvale CA
 *         August 22, 1988
 * @author Modernization by
 *         Nicholas DeCicco
 *         Rowan University
 *         <deciccon0@students.rowan.edu>
 *         <nsd.cicco@gmail.com>
 */
int dihf (int gh)
{
	int ios;
	int j;
	double sn;
	double h2;
	double hpx;
	double argument, argument2;

	ios = gh;
	sn = 0.0001;

	switch(gh)
	{
		case 3:
			for (j = 1; j <= 1; j++)
			{
				h2 = x*x + y*y;
				argument = h2;
				h = sqrt(argument);       /* calculate horizontal intensity */
				argument = h2 + z*z;
				f = sqrt(argument);      /* calculate total intensity */
				if (f < sn)
				{
					d = NaN;        /* If d and i cannot be determined, */
					i = NaN;        /*       set equal to NaN         */
				}
				else
				{
					argument = z;
					argument2 = h;
					i = atan2(argument,argument2);
					if (h < sn)
					{
						d = NaN;
					}
					else
					{
						hpx = h + x;
						if (hpx < sn)
						{
							d = PI;
						}
						else
						{
							argument = y;
							argument2 = hpx;
							d = 2.0 * atan2(argument,argument2);
						}
					}
				}
			}
			break;
		case 4:
			for (j = 1; j <= 1; j++)
			{
				h2 = xtemp*xtemp + ytemp*ytemp;
				argument = h2;
				htemp = sqrt(argument);
				argument = h2 + ztemp*ztemp;
				ftemp = sqrt(argument);
				if (ftemp < sn)
				{
					dtemp = NaN;    /* If d and i cannot be determined, */
					itemp = NaN;    /*       set equal to 999.0         */
				}
				else
				{
					argument = ztemp;
					argument2 = htemp;
					itemp = atan2(argument,argument2);
					if (htemp < sn)
					{
						dtemp = NaN;
					}
					else
					{
						hpx = htemp + xtemp;
						if (hpx < sn)
						{
							dtemp = PI;
						}
						else
						{
							argument = ytemp;
							argument2 = hpx;
							dtemp = 2.0 * atan2(argument,argument2);
						}
					}
				}
			}
			break;
		default:
			fprintf(stderr, "\nError in subroutine dihf");
			break;
	}
	return ios;
}
