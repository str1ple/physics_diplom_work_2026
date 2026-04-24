import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'

# для легенд
TEMP = 300          # K
FREQ = 1.0          # ТГц
DC_FIXED = 500      # В/м
AC_FIXED = 1000     # В/м
PHASE_FIXED = 45    # град

# пост поле
data1 = np.loadtxt('dc_only_transverse.txt', comments='#')
E_dc = data1[:, 0]
v_perp1 = data1[:, 1]

plt.figure(figsize=(8, 6))
plt.plot(E_dc, v_perp1, 'o-', color='darkblue', lw=2, ms=6,
         label=f'Поперечная скорость\n(без переменного поля)')
plt.xlabel('Постоянное электрическое поле $E_{dc}$ (В/м)')
plt.ylabel('Средняя поперечная скорость $v_\\perp$ (м/с)')
plt.title(f'Зависимость поперечной скорости от постоянного поля\n'
          f'$E_x = E_y$, $T = {TEMP}$ К')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('dc_only_transverse.png', dpi=300)
plt.show()

# перем поле
data2 = np.loadtxt('ac_amplitude_transverse.txt', comments='#')
E0 = data2[:, 0]
v_perp2 = data2[:, 1]

plt.figure(figsize=(8, 6))
plt.plot(E0, v_perp2, 's-', color='crimson', lw=2, ms=6,
         label=f'С переменным полем\n'
               f'$E_{{dc}} = {DC_FIXED}$ В/м, $\\phi = {PHASE_FIXED}°$')
plt.xlabel('Амплитуда переменного поля $E_0$ (В/м)')
plt.ylabel('Средняя поперечная скорость $v_\\perp$ (м/с)')
plt.title(f'Влияние амплитуды переменного поля на поперечную скорость\n'
          f'$T = {TEMP}$ К, $f = {FREQ}$ ТГц')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('ac_amplitude_transverse.png', dpi=300)
plt.show()

# сравнение
plt.figure(figsize=(9, 6))
# только пост поле
plt.plot(E_dc, v_perp1, 'o-', color='darkblue', lw=2, ms=6,
         label='Только постоянное поле\n(изменяется $E_{dc}$)')
# оба поля
plt.plot(E0, v_perp2, 's-', color='crimson', lw=2, ms=6,
         label=f'Постоянное + переменное поле\n'
               f'$E_{{dc}} = {DC_FIXED}$ В/м, $\\phi = {PHASE_FIXED}°$')
plt.xlabel('Напряжённость поля (В/м)\n'
           '(для кривой 1 – $E_{dc}$, для кривой 2 – амплитуда $E_0$)')
plt.ylabel('Средняя поперечная скорость $v_\\perp$ (м/с)')
plt.title(f'Сравнение поперечной скорости с переменным полем и без него\n'
          f'$T = {TEMP}$ К, $f = {FREQ}$ ТГц')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('comparison_transverse.png', dpi=300)
plt.show()

# фаза
data3 = np.loadtxt('phase_transverse.txt', comments='#')
phi_rad = data3[:, 0]
phi_deg = data3[:, 1]
v_perp3 = data3[:, 2]

plt.figure(figsize=(10, 5))
plt.plot(phi_deg, v_perp3, 'o-', color='green', lw=2, ms=6,
         label=f'$E_{{dc}} = {DC_FIXED}$ В/м, $E_0 = {AC_FIXED}$ В/м')
plt.xlabel('Фаза $\\phi$ (градусы)')
plt.ylabel('Средняя поперечная скорость $v_\\perp$ (м/с)')
plt.title(f'Зависимость поперечной скорости от фазы переменного поля\n'
          f'$T = {TEMP}$ К, $f = {FREQ}$ ТГц')
plt.xticks(np.arange(0, 361, 45))
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('phase_transverse_linear.png', dpi=300)
plt.show()

# поляр диаг
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='polar')
ax.plot(phi_rad, v_perp3, 'o-', color='green', lw=2, ms=6,
        label=f'$E_{{dc}} = {DC_FIXED}$ В/м\n$E_0 = {AC_FIXED}$ В/м')
ax.set_title('Полярная диаграмма: поперечная скорость vs фаза', pad=20)
ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
plt.tight_layout()
plt.savefig('phase_transverse_polar.png', dpi=300)
plt.show()