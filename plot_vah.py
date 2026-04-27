import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# Настройки стиля
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11

def load_comparison_data(pattern):
    files = sorted(glob.glob(pattern))
    data = {}
    for f in files:
        name = os.path.basename(f)
        # Извлекаем значение параметра из имени файла (например, E0_500 -> 500)
        parts = name.replace('.txt', '').split('_')
        val_str = parts[-1]
        try:
            val = float(val_str)
        except:
            val = val_str
        data[val] = np.loadtxt(f, comments='#')
    return data

def plot_comparison_series(fig, ax, data, title, xlabel, param_name, param_unit, is_freq=False):
    # Сортируем ключи для порядка в легенде
    keys = sorted(data.keys(), key=lambda x: float(str(x).replace('_', '.')))

    # Цвета и маркеры
    colors = plt.cm.tab10(np.linspace(0, 1, len(data)))
    markers_dc = ['o', 's', '^', 'D', 'p', 'h'] # Для чистого DC
    markers_ac = ['o', 's', '^', 'D', 'p', 'h'] # Для DC+AC

    # Сначала рисуем ВСЕ линии "Только DC", чтобы они были на заднем плане (или одинаковые)
    # Но так как DC only одинаков для всех, нарисуем его один раз жирным серым/синим
    first_key = keys[0]
    d_ref = data[first_key]
    ax.plot(d_ref[:,0], d_ref[:,1], color='navy', linewidth=2.5, linestyle='-',
            marker='o', markersize=6, label='Только постоянное поле (DC only)', zorder=10)

    # Теперь рисуем линии DC+AC для каждого параметра
    for i, val in enumerate(keys):
        d = data[val]
        color = colors[i]
        label_val = val if not is_freq else f"{val}"

        ax.plot(d[:,0], d[:,2], color=color, linewidth=2, linestyle='--',
                marker=markers_ac[i%len(markers_ac)], markersize=6,
                label=f'{param_name} = {label_val} {param_unit}')

    ax.set_xlabel('$E_{dc}$ (В/м)\n(напряжённость постоянного поля)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Средняя поперечная скорость $v_perp$ (м/с)', fontweight='bold', fontsize=12)
    ax.set_title(title, fontweight='bold', fontsize=14, pad=15)

    # Легенда
    ax.legend(loc='best', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

# ================= СЕРИЯ 1: Амплитуда E0 =================
print("Plotting Series 1: Comparison by AC Amplitude...")
data_s1 = load_comparison_data("comp_s1_E0_*.txt")
fig1, ax1 = plt.subplots(figsize=(10, 7))
plot_comparison_series(fig1, ax1, data_s1,
                       'Влияние амплитуды переменного поля $E_0$\nна поперечную скорость',
                       '$E_{dc}$ (В/м)', '$E_0$', 'В/м')
fig1.savefig('series1_compare_amplitude.png', dpi=300, bbox_inches='tight')
plt.show()

# ================= СЕРИЯ 2: Частота f =================
print("Plotting Series 2: Comparison by AC Frequency...")
data_s2 = load_comparison_data("comp_s2_f_*.txt")
fig2, ax2 = plt.subplots(figsize=(10, 7))
plot_comparison_series(fig2, ax2, data_s2,
                       'Влияние частоты переменного поля $f$\nна поперечную скорость',
                       '$E_{dc}$ (В/м)', '$10*f$', 'ТГц', is_freq=True)
fig2.savefig('series2_compare_frequency.png', dpi=300, bbox_inches='tight')
plt.show()

# ================= СЕРИЯ 3: Фаза phi =================
print("Plotting Series 3: Comparison by AC Phase...")
data_s3 = load_comparison_data("comp_s3_phi_*.txt")
fig3, ax3 = plt.subplots(figsize=(10, 7))
plot_comparison_series(fig3, ax3, data_s3,
                       'Влияние фазы переменного поля $phi$\nна поперечную скорость',
                       '$E_{dc}$ (В/м)', '$phi$', '°')
fig3.savefig('series3_compare_phase.png', dpi=300, bbox_inches='tight')
plt.show()

print("Done. Figures saved as series*_compare_*.png")