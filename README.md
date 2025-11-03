# Математическое моделирование врождённого иммунного ответа при вирусной инфекции с учётом интерферона (Модель Павелек)

Данный проект реализует численное моделирование динамики взаимодействия вируса и клеток организма с учётом действия интерферона. В основе лежит [модель Павелек (Pawelek et al.)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3386161/), описывающая изменения во времени количества здоровых, инфицированных, рефрактерных клеток, вирусных частиц и молекул интерферона.

---

## Структура проекта

```
├── Parameters.py          # Файл с параметрами модели
├── Pawelek_model.py       # Основной скрипт моделирования и построения графиков
├── results_pony/          # Папка, в которую сохраняются все результаты моделирования
│   ├── viral_load_comparison_pony.png
│   └── cell_populations_pony.png
```

---

##  Файлы проекта

### `Parameters.py`
Содержит все параметры модели в виде словаря `model_params`:

```python
model_params = {
    'b': 1.9e-4,   # коэффициент инфицирования
    'w': 2.5e-2,   # скорость перехода клеток в рефрактерное состояние под действием интерферона
    'r': 2.0e-1,   # скорость восстановления клеток
    'dI': 2,       # скорость гибели инфицированных клеток
    'k': 2.4,      # эффект интерферона на гибель инфицированных клеток
    'p': 8.1e-6,   # скорость продукции вируса
    'c': 5.8,      # скорость очистки вируса
    'q': 2.1e-9,   # скорость продукции интерферона
    'd': 2.4,      # скорость распада интерферона
    's': 2.2,      # коэффициент экспоненциального роста гибели клеток после порога mu
    'mu': 6,       # порог времени, после которого dI начинает расти экспоненциально
    'T0': 3.5e11,  # начальное число чувствительных клеток
    'I0': 0,       # начальное число инфицированных клеток
    'R0': 0,       # начальное число рефрактерных клеток
    'V0': 10,      # начальное количество вируса
    'F0': 1        # начальный уровень интерферона
}
```

---

### `Pawelek_model.py`
Основной модуль моделирования.  
Реализует два варианта модели:

1. **Pawelek model (ур. 1)** — учитывает влияние интерферона на гибель инфицированных клеток.  
2. **Saenz model** — упрощённый вариант без учёта этого влияния.

Обе модели решаются численно с помощью `scipy.integrate.odeint`.

Скрипт выполняет следующие шаги:
1. Загружает параметры из `Parameters.py`.
2. Определяет систему дифференциальных уравнений для обеих моделей.
3. Интегрирует их на временном интервале `0–10` дней.
4. Сравнивает результаты моделей с экспериментальными данными по вирусной нагрузке.
5. Строит и сохраняет два графика:
   - `viral_load_comparison_pony.png` — сравнение моделирования с данными эксперимента;
   - `cell_populations_pony.png` — изменение популяций клеток во времени.

---

## Результаты

После запуска программы автоматически создаётся папка `results_pony`, где сохраняются:

- **`viral_load_comparison_pony.png`**  
  График сравнения модели Павелек и модели Саенса с экспериментальными данными по вирусной нагрузке.

- **`cell_populations_pony.png`**  
  Динамика изменения числа чувствительных, инфицированных и рефрактерных клеток.

---

## Как запустить проект

### 1. Установите зависимости
Убедитесь, что у вас установлены необходимые библиотеки Python,
все необходимые библиотеки указаны в файле `requirements.txt`.  
```bash
pip install -r requirements.txt
```

### 2. Запустите основной скрипт
```bash
python Pawelek_model.py
```

После выполнения в консоли появятся сообщения:
```
Graph 1 saved: results_pony/viral_load_comparison_pony.png
Graph 2 saved: results_pony/cell_populations_pony.png
```

## Источники

1.	[Modeling Within-Host Dynamics of Influenza Virus Infection Including Immune Responses Kasia A. Pawelek,Giao T. Huynh,Michelle Quinlivan,Ann Cullinane,Libin Rong, Alan S. Perelson, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3386161/)
2.	[Interferon Induction and/or Production and Its Suppression by Influenza A Viruses, Philip I Marcus, Jillian M Rojek, Margaret J Sekellick, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC548469/)
3.	[3.	Dynamics of influenza virus infection and pathology Roberto A Saenz, Michelle Quinlivan, Debra Elton, Shona Macrae, Anthony S Blunden, Jennifer A Mumford, Janet M Daly, Paul Digard, Ann Cullinane, Bryan T Grenfell, John W McCauley, James L N Wood, Julia R Gog, 2010](https://pubmed.ncbi.nlm.nih.gov/20130053/)
