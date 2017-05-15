import numpy as np
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel


files = np.sort([ 'gradient_depth358'])
turning_points = np.zeros((len(files), len(distance)))

for d in range(len(distance)):
    for f in range(len(files)):

        print(files[f])
        model = TauPyModel(model=str(files[f]))
        arrivals = model.get_pierce_points(source_depth_in_km= 250, distance_in_degree= distance[d], phase_list= ["S"])
        s_arr = arrivals[0]
        turning_point = np.max(s_arr.pierce['depth'])
        turning_points[f, d] = turning_point

        print(turning_point)

u = np.arange(1, 7, 1)

for j in range(len(distance)):

    plt.plot(u, turning_points[:, j], 'o--', label = '{}' .format(distance[j]))


##Need to update to define actual figure for publication quality

plt.xlabel('Perturbation Distance Above CMB (100km)')
plt.ylabel('Turning Point Depth (km)')
plt.legend(loc = 0)
plt.xlim([0, 7])
plt.show()
