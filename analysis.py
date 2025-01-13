import numpy as np
import matplotlib.pyplot as plt
import os

directory = 'asymm_mesh_2'

if not os.path.exists('images/' + directory):
    os.makedirs('images/' + directory)
    print(f"Directory 'images/{directory}' created.")


filename1 = os.path.join('results', directory, 'residual_first_iteration.csv')
filename2 = os.path.join('results', directory, 'relative_distances.csv')

residuals = np.loadtxt(filename1, delimiter=',')
relative_distances = np.loadtxt(filename2, delimiter=',')

initial_index = 2 # 110 or 2, depending on the fact that continuation is used
residuals = residuals[initial_index:]
relative_distances = relative_distances[initial_index:]

plt.plot(residuals, label='Residuals')
plt.yscale('log')
plt.xlabel('Time step')
plt.ylabel('Value')
plt.title('Residuals at first iteration')
plt.legend()
plt.savefig('images/' + directory + '/residuals.png')
plt.show()
plt.close()
plt.clf()

plt.plot(relative_distances, label='Relative distance')
plt.yscale('log')
plt.xlabel('Time step')
plt.ylabel('Value')
plt.title('Relative distances')
plt.legend()
plt.savefig('images/' + directory + '/relative_distances.png')
plt.show()
plt.close()
plt.clf()
