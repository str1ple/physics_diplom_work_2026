import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('VAH_data.txt', comments='#')
E = data[:, 0]
v_y = data[:, 1]
v_x = data[:, 2]

plt.figure(figsize=(10, 6))
plt.plot(E, np.abs(v_y), 'bo-', linewidth=2, markersize=8, label='v_y (продольная)')
plt.plot(E, np.abs(v_x), 'rs-', linewidth=2, markersize=8, label='v_x (поперечная)')

plt.xlabel('Electric Field E (V/m)', fontsize=12)
plt.ylabel('Average Velocity (m/s)', fontsize=12)
plt.title('Velocity-Field Characteristic for Superlattice', fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)

plt.savefig('velocity_field.png', dpi=300)
plt.show()