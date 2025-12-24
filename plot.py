import matplotlib
# Используем бэкенд, не требующий Qt/Wayland
matplotlib.use('Agg')  # <-- Добавьте эту строку ДО импорта pyplot
import matplotlib.pyplot as plt

with open("stats.tsv") as f:
    ths = []
    times = []
    
    # Пропускаем заголовок, если он есть
    header = f.readline()
    
    for line in f:
        if line.strip():  # пропускаем пустые строки
            try:
                th, time = line.strip().split("\t")
                ths.append(int(th))
                times.append(float(time))
            except ValueError as e:
                print(f"Ошибка в строке: {line.strip()}")
                print(f"Сообщение: {e}")

# Проверяем, что данные есть
if not ths:
    print("Нет данных для построения графика!")
    exit(1)

fig = plt.figure(figsize=(10, 6))
plt.plot(ths, times, 'o-', linewidth=2, markersize=8)
plt.xlabel("Количество потоков (th)")
plt.ylabel("Время выполнения (сек)")
plt.title("Зависимость времени выполнения от количества потоков")
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Сохраняем с лучшим качеством
fig.savefig("th_x_time.jpg", dpi=150)
print("График сохранен как th_x_time.jpg")