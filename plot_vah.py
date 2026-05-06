import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from scipy.interpolate import CubicSpline
from matplotlib.ticker import MaxNLocator

# Настройки стиля
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 11

# Константы моделирования для легенды
COMMON_PARAMS = {
    'T': 300,          # К
    't_mod': 200,      # пс
    'N': 2000          # число частиц
}

def load_comparison_data(pattern):
    files = sorted(glob.glob(pattern))
    data = {}
    for f in files:
        name = os.path.basename(f)
        # Извлекаем значение параметра из имени файла
        parts = name.replace('.txt', '').split('_')
        val_str = parts[-1]
        try:
            val = float(val_str)
        except:
            val = val_str
        data[val] = np.loadtxt(f, comments='#')
    return data

def add_gost_axes(ax, x_max, y_min, y_max, xlabel_text, ylabel_text):
    # Скрываем стандартную рамку
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # ось X со стрелкой
    # ИЗМЕНЕНИЕ: Начинаем линию строго от x=0, чтобы не перечеркивать ноль на оси Y
    ax.annotate('', xy=(x_max * 1.03, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))

    # ось Y со стрелкой (начинается от y_min и идет выше y_max)
    y_arrow_tip = y_max * 1.05
    ax.annotate('', xy=(0, y_arrow_tip), xytext=(0, y_min),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))

    # Подпись оси X (справа, у стрелки)
    ax.text(x_max * 1.05, -y_max * 0.04, xlabel_text,
            ha='left', va='center', fontsize=13, fontweight='regular')

    # Подпись оси Y (ГОРИЗОНТАЛЬНО, сверху оси)
    ax.text(-x_max * 0.05, y_arrow_tip * 1.05, ylabel_text,
            ha='center', va='bottom', rotation=0, fontsize=13, fontweight='regular')

    # Настройка тиков: убираем черточки (length=0), оставляем только текст
    ax.tick_params(axis='both', which='major', labelsize=11, direction='out', length=0)

def plot_comparison_series_gost(fig, ax, data, title, param_name, param_unit, is_freq=False, series_type='amp'):
    keys = sorted(data.keys(), key=lambda x: float(str(x).replace('_', '.')))
    colors = plt.cm.tab10(np.linspace(0, 1, len(data)))
    markers = ['o', 's', '^', 'D', 'p', 'h']

    # --- 1. Подготовка данных для определения масштаба осей ---
    all_x = []
    all_y = []
    for val in keys:
        d = data[val]
        all_x.extend(d[:, 0])
        all_y.extend(d[:, 1]) # DC only
        all_y.extend(d[:, 2]) # DC + AC

    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)

    # Отступы
    x_range = x_max - x_min if x_max != x_min else 1
    y_range = y_max - y_min if y_max != y_min else 1

    # Увеличиваем диапазоны для осей, чтобы влезли стрелки и подписи
    x_plot_max = x_max + 0.15 * x_range
    y_plot_max = y_max + 0.15 * y_range
    y_plot_min = y_min - 0.1 * y_range

    # --- 2. Построение кривых (Сплайн + Маркеры) ---

    legend_handles = []
    legend_labels = []

    # кривая "Только DC"
    first_key = keys[0]
    d_ref = data[first_key]
    x_ref = d_ref[:, 0]
    y_ref = d_ref[:, 1]

    # Сплайн
    x_smooth_ref = np.linspace(x_ref.min(), x_ref.max(), 300)
    cs_ref = CubicSpline(x_ref, y_ref)
    y_smooth_ref = cs_ref(x_smooth_ref)

    # Линия
    line_ref, = ax.plot(x_smooth_ref, y_smooth_ref, color='navy', linewidth=2.5, linestyle='-', zorder=5)
    # Точки
    ax.plot(x_ref, y_ref, color='navy', marker='o', linestyle='', markersize=7, zorder=6)

    # ДОБАВЛЯЕМ В ЛЕГЕНДУ
    legend_handles.append(line_ref)
    legend_labels.append('Только постоянное поле $E_{dc}$')

    # Кривые "DC + AC"
    for i, val in enumerate(keys):
        d = data[val]
        x = d[:, 0]
        y = d[:, 2]

        color = colors[i]
        label_val = val if not is_freq else f"{val}"

        # Сплайн
        x_smooth = np.linspace(x.min(), x.max(), 300)
        cs = CubicSpline(x, y)
        y_smooth = cs(x_smooth)

        # Линия
        line, = ax.plot(x_smooth, y_smooth, color=color, linewidth=2, linestyle='--', zorder=3)
        # Точки
        ax.plot(x, y, color=color, marker=markers[i%len(markers)], linestyle='', markersize=7, zorder=4)

        legend_handles.append(line)
        legend_labels.append(f'{param_name} = {label_val} {param_unit}')

    # --- 3. Оформление осей ---
    add_gost_axes(ax, x_plot_max, y_plot_min, y_plot_max,
                  xlabel_text='$E_{dc}$, В/м',
                  ylabel_text='$v_\perp$, м/с')

    ax.set_title(title, fontsize=15, pad=25, fontweight='regular')

    # Сетка
    ax.grid(True, alpha=0.5, color='gray', linestyle='--')

    # Устанавливаем лимиты, оставляя место слева для подписей
    ax.set_xlim(-x_plot_max * 0.05, x_plot_max * 1.15)
    ax.set_ylim(y_plot_min, y_plot_max * 1.2)

    # ИЗМЕНЕНИЕ: Делаем деления (ticks) более частыми, особенно по оси Y
    ax.xaxis.set_major_locator(MaxNLocator(nbins=8))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=12)) # Увеличено число делений

    # Применяем локаторы к холсту, чтобы сформировать значения
    fig.canvas.draw()

    # Фильтруем значения, убирая вылазящие за концы и отрицательные по оси X
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()

    ax.set_xticks([t for t in xticks if 0 <= t < x_plot_max])
    ax.set_yticks([t for t in yticks if t < y_plot_max])

    # --- 4. Формирование окна легенды ---

    # дополнительные параметры для текущей серии
    extra_params = ""
    if series_type == 'amp':
        extra_params = f"$f = 1.0$ ТГц, $\phi = 90^\circ$"
    elif series_type == 'freq':
        extra_params = f"$E_0 = 3000$ В/м, $\phi = 90^\circ$"
    elif series_type == 'phase':
        extra_params = f"$E_0 = 3000$ В/м, $f = 1.0$ ТГц"

    common_text = f"$T = {COMMON_PARAMS['T']}$ К, $t_{{mod}} = {COMMON_PARAMS['t_mod']}$ пс, $N = {COMMON_PARAMS['N']}$"
    params_text = f"Параметры:\n{common_text}\n{extra_params}"

    from matplotlib.patches import Patch
    patch = Patch(color='none', edgecolor='none')
    legend_handles.append(patch)
    legend_labels.append(params_text)

    # Размещаем легенду справа от графика
    leg = ax.legend(handles=legend_handles, labels=legend_labels,
                    loc='center left', bbox_to_anchor=(1.05, 0.5),
                    frameon=True, shadow=True, borderpad=0.8,
                    handlelength=2.5, handletextpad=0.8)

    text_params = leg.get_texts()[-1]
    text_params.set_fontweight('regular')
    text_params.set_fontsize(12)
    text_params.set_color('darkred')

    # Явно задаем отступы, чтобы легенда точно влезла и график не обрезался
    fig.subplots_adjust(left=0.08, right=0.65, top=0.85, bottom=0.1)

# ================= СЕРИЯ 1: Амплитуда E0 =================
print("Plotting Series 1: Comparison by AC Amplitude (GOST)...")
data_s1 = load_comparison_data("comp_s1_E0_*.txt")
if data_s1:
    fig1, ax1 = plt.subplots(figsize=(15, 8))
    plot_comparison_series_gost(fig1, ax1, data_s1,
                                'Влияние амплитуды переменного поля $E_0$\nна поперечную скорость',
                                '$E_0$', 'В/м', series_type='amp')
    fig1.savefig('series1_compare_amplitude_gost.png', dpi=300, bbox_inches='tight')
    plt.show()

# ================= СЕРИЯ 2: Частота f =================
print("Plotting Series 2: Comparison by AC Frequency (GOST)...")
data_s2 = load_comparison_data("comp_s2_f_*.txt")
if data_s2:
    fig2, ax2 = plt.subplots(figsize=(15, 8))
    plot_comparison_series_gost(fig2, ax2, data_s2,
                                'Влияние частоты переменного поля $f$\nна поперечную скорость',
                                '$f$', 'ТГц', is_freq=True, series_type='freq')
    fig2.savefig('series2_compare_frequency_gost.png', dpi=300, bbox_inches='tight')
    plt.show()

# ================= СЕРИЯ 3: Фаза phi =================
print("Plotting Series 3: Comparison by AC Phase (GOST)...")
data_s3 = load_comparison_data("comp_s3_phi_*.txt")
if data_s3:
    fig3, ax3 = plt.subplots(figsize=(15, 8))
    plot_comparison_series_gost(fig3, ax3, data_s3,
                                'Влияние фазы переменного поля $\phi$\nна поперечную скорость',
                                '$\phi$', '°', series_type='phase')
    fig3.savefig('series3_compare_phase_gost.png', dpi=300, bbox_inches='tight')
    plt.show()

print("Готово!")