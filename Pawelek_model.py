import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import os

from Parameters import model_params

output_dir = "results_pony"
os.makedirs(output_dir, exist_ok=True)

def pawelek_model(y, t, p):
    T, I, R, V, F = y
    d = p['dI']
    if t > p['mu']:
        d = p['dI'] * np.exp(p['s'] * (t - p['mu']))
    dTdt = -p['b'] * V * T - p['w'] * F * T + p['r'] * R
    dIdt = p['b'] * V * T - d * I - p['k'] * I * F
    dRdt = p['w'] * F * T - p['r'] * R
    dVdt = p['p'] * I - p['c'] * V
    dFdt = p['q'] * I - p['d'] * F
    return [dTdt, dIdt, dRdt, dVdt, dFdt]

def saenz_model(y, t, p):
    T, I, R, V, F = y
    d = p['dI']
    if t > p['mu']:
        d = p['dI'] * np.exp(p['s'] * (t - p['mu']))
    dTdt = -p['b'] * V * T - p['w'] * F * T + p['r'] * R
    dIdt = p['b'] * V * T - d * I  # No k * I * F
    dRdt = p['w'] * F * T - p['r'] * R
    dVdt = p['p'] * I - p['c'] * V
    dFdt = p['q'] * I - p['d'] * F
    return [dTdt, dIdt, dRdt, dVdt, dFdt]

y0 = [model_params['T0'], model_params['I0'], model_params['R0'], model_params['V0'], model_params['F0']]

t = np.linspace(0, 10, 1000)

solution_pawelek = odeint(pawelek_model, y0, t, args=(model_params,))
T_p, I_p, R_p, V_p, F_p = solution_pawelek.T

solution_saenz = odeint(saenz_model, y0, t, args=(model_params,))
T_s, I_s, R_s, V_s, F_s = solution_saenz.T

V_p = np.maximum(V_p, 1e-10)  # Avoid log10(0)
V_s = np.maximum(V_s, 1e-10)

days = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
V_data_log10_approx = np.array([3, 6, 3, 4, 3, 3, 2, 0, 0, 0])
V_data = 10 ** V_data_log10_approx

plt.figure(figsize=(10, 6))
plt.plot(t, np.log10(V_p), 'r-', linewidth=2, label='Ур. (1)')
plt.plot(t, np.log10(V_s), 'g--', linewidth=2, label='Модель Саенса')
plt.plot(days, V_data_log10_approx, 'ro', markersize=8, label='Данные по вирусной нагрузке')
plt.axhline(np.log10(100), ls='--', color='b', label='Предел обнаружения (100 копий РНК/мл)')
plt.xlabel('Время (дни)')
plt.ylabel(r'$\log_{10}$ копий NS мл$^{-1}$')
plt.xlim(0, 10)
plt.ylim(0, 8)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_position(('axes', 0))
ax.spines['bottom'].set_position(('axes', 0))
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.title('Сравнение моделей с данными по вирусной нагрузке для Пони 1')
plt.legend()
plt.grid(True, ls="--")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'viral_load_comparison_pony.png'), dpi=300)
plt.close()
print(f"Graph 1 saved: {os.path.join(output_dir, 'viral_load_comparison_pony.png')}")

scale = 1e11
plt.figure(figsize=(10, 6))
plt.plot(t, T_p / scale, 'b-', linewidth=2, label='Чувствительные клетки')
plt.plot(t, I_p / scale, 'g-', linewidth=2, label='Инфицированные клетки')
plt.plot(t, R_p / scale, 'r--', linewidth=2, label='Рефрактерные клетки')
plt.plot(t, (T_p + I_p + R_p) / scale, 'k:', linewidth=2, label='Общее количество клеток')
plt.xlabel('Время (дни)')
plt.ylabel(r'Клетки ($\times 10^{11}$)')
plt.xlim(0, 10)
plt.ylim(0, 3.5)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_position(('axes', 0))
ax.spines['bottom'].set_position(('axes', 0))
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.title('Изменения популяций клеток, предсказанные Ур. (1) для Пони 1')
plt.legend()
plt.grid(True, ls="--")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'cell_populations_pony.png'), dpi=300)
plt.close()
print(f"Graph 2 saved: {os.path.join(output_dir, 'cell_populations_pony.png')}")