import numpy as np

lons = np.arange(-180, 190, 40)
lats = np.arange(0, 190, 40)

num_stations = len(lons)*len(lats)

header0 = str(num_stations)
with open('STATIONS', 'w') as f:

    f.write(header0)
    f.write('\n')

    for i in lats:

        for j in lons:

            lat = str(i)

            lon = str(j)

            f.write('{} {}' .format(lat, lon + '\n'))
