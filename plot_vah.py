import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('VAH_data.txt', comments='#')
E = data[:, 0]
v_y = data[:, 1]
j = data[:, 2]

plt.figure(figsize=(10, 6))
plt.plot(E, np.abs(v_y), 'bo-', linewidth=2, markersize=8)
plt.xlabel('Electric Field E (V/m)', fontsize=12)
plt.ylabel('Average Velocity v_y (m/s)', fontsize=12)
plt.title('Velocity-Field Characteristic for Superlattice', fontsize=14)
plt.grid(True, alpha=0.3)
plt.savefig('velocity_field.png', dpi=300)
plt.show()