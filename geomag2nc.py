
# Author: Nicholas DeCicco <nsd.cicco@gmail.com>
#                          <deciccon0@students.rowan.edu>

from subprocess import call
import io
from numpy import *
from netCDF4 import Dataset

dataset = "IGRF12.COF"
date = "2016,10,10"
coordSys = 'D'
height = 'K500'

n_lat = 40
n_lon = 90

coordsFileName = 'coords.txt'
outFileName = 'out.txt'

coordsFile = io.open(coordsFileName, 'w')

# We can't get the H field at the poles, so only go within 'epsilon' (degrees)
# of them.
epsilon = 10

latValues = linspace(-90+epsilon, 90-epsilon, n_lat)
lonValues = linspace(-180, 180, n_lon)
latGrid,lonGrid = meshgrid(latValues, lonValues)

for i in range(0,n_lon):
    for j in range(0,n_lat):
        line = '{} {} {} {} {}\n'.format(date, coordSys, height,
                                         latValues[j], lonValues[i])
        coordsFile.write(line)

coordsFile.close()

call(["geomag70.exe", dataset, 'f', coordsFileName, outFileName])

n = r"(-?\d(?:\.?\d*)+)" # number
d = r"(-?\d+)d" # minute
m = r"(-?\d+)m" # second
regex = r"\s+".join([r"(\d+),(\d+),(\d+)", # Date
                     r"[CD]",              # Coordinate system
                     r"[KMF]\d+",          # Altitude
                     n,                    # Latitude
                     n,                    # Longitude
                     d, m,                 # Declination (D)
                     d, m,                 # Inclination (I)
                     n,                    # Horizontal field strength (H)
                     n,                    # North component (X)
                     n,                    # East component (Y)
                     n,                    # Down component (Z)
                     n,                    # Total field strength (F)
                     n,                    # dD/dt
                     n,                    # dI/dt
                     n,                    # dH/dt
                     n,                    # dX/dt
                     n,                    # dY/dt
                     n,                    # dZ/dt
                     n                     # dF/dt
                     ]) + '.*'
variables = ['h', 'x', 'y', 'z', 'f',
             'ddot', 'idot', 'hdot', 'xdot', 'ydot', 'zdot', 'fdot']
data = fromregex('out.txt', regex, dtype=[('year', int), ('month', int), ('day', int),
                                          ('lat', float), ('lon', float),
                                          ('D_deg', float), ('D_min', float),
                                          ('I_deg', float), ('I_min', float)] +
                                         [(var, float) for var in variables])

try:
    out = Dataset('out.nc', 'w', format='NETCDF4')

    latDim = out.createDimension('lat', n_lat)
    lonDim = out.createDimension('lon', n_lon)
    lat = out.createVariable('lat', 'f8', ('lat'))
    lon = out.createVariable('lon', 'f8', ('lon'))

    lat.units = 'degrees_north'
    lon.units = 'degrees_east'

    lat[:] = latValues
    lon[:] = lonValues

    D = out.createVariable('d', 'f8', ('lon', 'lat'))
    I = out.createVariable('i', 'f8', ('lon', 'lat'))

    D[:] = data['D_deg'] + data['D_min']/60.0
    I[:] = data['I_deg'] + data['I_min']/60.0

    for varName in variables:
        var = out.createVariable(varName, 'f8', ('lon', 'lat'))
        var[:] = data[varName]

finally:
    out.close()
