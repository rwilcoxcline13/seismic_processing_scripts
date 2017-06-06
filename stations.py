import numpy as np

lons = np.arange(-180, 190, 40)
lats = np.arange(-30, 65, 5)
colats = 90-lats


station_names = []

with open('station_names.txt', 'r') as f1:

    lines = f1.readlines()

    for x in lines:

        station_names.append(x.split('\n')[0])





k = 0

with open('STATIONS', 'w') as f:

    for i in colats:

        for j in lons:

            colat = str(i)

            lon = str(j)

            k = k+1

            network = 'II'

            station = station_names[k]

            f.write('{}      {}       {}       {}      {}       {}' .format(str(station), network, colat, lon, str(0), str(0) + '\n'))
