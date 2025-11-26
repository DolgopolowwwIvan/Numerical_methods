import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_splines_full():
    # Чтение данных
    try:
        data = pd.read_csv('output/spline_comparison.csv')
    except Exception as e:
        print(f"Ошибка загрузки: {e}")
        return
    
    # Автоматически определяем колонки сплайнов
    all_columns = list(data.columns)
    spline_columns = [col for col in all_columns if col not in ['Index', 'Original']]
    
    if not spline_columns:
        print("Не найдены данные сплайнов")
        return
    
    print(f"Найдены сплайны: {spline_columns}")
    
    # Создаем график
    plt.figure(figsize=(15, 8))
    
    # Для 1670 точек лучше использовать тонкие линии
    # Исходные данные (каждую 10-ю точку для лучшей визуализации)
    plt.scatter(data['Index'][::10], data['Original'][::10], 
                color='black', label='Исходные данные', 
                alpha=0.5, s=10, zorder=5)
    
    # Цвета для сплайнов
    colors = ['gray', 'green', 'orange', 'red']  
    line_widths = [1.0, 1.5, 2.0, 3.0]  
    
    # Рисуем сплайны (все точки)
    for i, col in enumerate(spline_columns):
        color = colors[i % len(colors)]
        
        # Определяем стиль линии
        if 'Interpolation' in col:
            linestyle = '--'
            linewidth = 1.0
            label = 'Интерполяционный сплайн'
        else:
            linestyle = '-'
            linewidth = 1.2
            label = col
        
        plt.plot(data['Index'], data[col], 
                color=color, linewidth=linewidth, 
                linestyle=linestyle, label=label,
                alpha=0.8)
    
    plt.xlabel('Номер случайного числа', fontsize=12)
    plt.ylabel('Сгенерированные случайные числа', fontsize=12)
    plt.title('Сравнение сплайнов (1670 точек)\nВариант 7: N=1670, M=1.04, σ=3.74', fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Сохраняем
    plt.savefig('output/spline_comparison_full.png', dpi=300, bbox_inches='tight')
    print("График сохранен: output/spline_comparison_full.png")
    plt.show()
    
    # Анализ отклонений
    print("\n=== Анализ отклонений (1670 точек) ===")
    for col in spline_columns:
        deviation = np.mean(np.abs(data[col] - data['Original']))
        print(f"{col}: {deviation:.6f}")

if __name__ == "__main__":
    plot_splines_full()
