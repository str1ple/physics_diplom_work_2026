import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-v0_8-whitegrid')

# пост поле
dc_data = np.loadtxt('dc_field_dependence.txt', comments='#')
E_dc = dc_data[:, 0]
v_x_dc, v_y_dc = dc_data[:, 2], dc_data[:, 3]

plt.figure(figsize=(9, 6))
plt.plot(E_dc, np.abs(v_x_dc), 'bo-', lw=2, ms=6, label='v_x')
plt.plot(E_dc, np.abs(v_y_dc), 'rs-', lw=2, ms=6, label='v_y')
plt.xlabel('DC Electric Field E (V/m)')
plt.ylabel('Average Velocity (m/s)')
plt.title('Velocity vs DC Field (with AC field)')
plt.legend()
plt.savefig('dc_field_plot.png', dpi=300)
plt.show()

# перем поле
ac_data = np.loadtxt('ac_field_dependence.txt', comments='#')
E0_ac = ac_data[:, 0]
v_x_ac, v_y_ac = ac_data[:, 2], ac_data[:, 3]

plt.figure(figsize=(9, 6))
plt.plot(E0_ac, np.abs(v_x_ac), 'bo-', lw=2, ms=6, label='v_x')
plt.plot(E0_ac, np.abs(v_y_ac), 'rs-', lw=2, ms=6, label='v_y')
plt.xlabel('AC Field Amplitude E₀ (V/m)')
plt.ylabel('Average Velocity (m/s)')
plt.title('Velocity vs AC Field Amplitude')
plt.legend()
plt.savefig('ac_field_plot.png', dpi=300)
plt.show()

# фаза
phase_data = np.loadtxt('phase_dependence.txt', comments='#')
phi = phase_data[:, 0]
v_x_p, v_y_p = phase_data[:, 2], phase_data[:, 3]

plt.figure(figsize=(12, 5))
plt.plot(phi, v_x_p, 'bo-', lw=2, ms=6, label='v_x')
plt.plot(phi, v_y_p, 'rs-', lw=2, ms=6, label='v_y')
plt.xlabel('Phase φ (rad)')
plt.ylabel('Average Velocity (m/s)')
plt.title('Velocity vs Phase')
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi], ['0', 'π/2', 'π', '3π/2', '2π'])
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('phase_plot.png', dpi=300)

# Polar
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='polar')
ax.plot(phi, v_x_p, 'b-', lw=2, label='v_x')
ax.plot(phi, v_y_p, 'r-', lw=2, label='v_y')
ax.set_title('Polar Plot: Velocity vs Phase', pad=20)
ax.legend(loc='upper right')
plt.savefig('phase_polar_plot.png', dpi=300)
plt.show()