import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import seaborn as sns

plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class SmoothingSplineAnalysis:
    
    def __init__(self, N=1670, M=1.04, sigma=3.74, seed=42):
       
        self.N = N
        self.M = M
        self.sigma = sigma
        self.seed = seed
        self.x = None
        self.y = None
        self.generate_data()
    
    def generate_data(self):
    
        np.random.seed(self.seed)
        self.x = np.arange(self.N)  # Номера наблюдений
        self.y = np.random.normal(self.M, self.sigma, self.N)
        print(f"Сгенерировано {self.N} наблюдений")
        print(f"Статистики: mean={np.mean(self.y):.3f}, std={np.std(self.y):.3f}")
    
    def create_splines(self, p_values=[0, 0.4, 0.8, 0.9]):
        self.p_values = p_values
        self.splines = {}
        
        for p in p_values:
            if p == 0:
                spline = UnivariateSpline(self.x, self.y, s=0)
            else:
                s = p * self.N * np.var(self.y)
                spline = UnivariateSpline(self.x, self.y, s=s)
            
            self.splines[p] = spline
            print(f"Создан сплайн с p={p}, s={spline._data[-1]:.2f}")
    
    def plot_splines_comparison(self):
        plt.figure(figsize=(14, 10))
             
        sample_size = min(200, self.N)
        indices = np.linspace(0, self.N-1, sample_size, dtype=int)
        x_vis = self.x[indices]
        y_vis = self.y[indices]
        
        sort_idx = np.argsort(x_vis)
        x_vis = x_vis[sort_idx]
        y_vis = y_vis[sort_idx]
           
        x_dense = np.linspace(self.x.min(), self.x.max(), 1000)
        

        plt.scatter(x_vis, y_vis, alpha=0.6, label='Исходные данные', 
                   color='gray', s=20)
        
        # Сплайны
        colors = ['red', 'blue', 'green', 'orange']
        for (p, spline), color in zip(self.splines.items(), colors):
            y_spline = spline(x_dense)
            plt.plot(x_dense, y_spline, color=color, linewidth=2, 
                    label=f'p = {p}')
        
        plt.xlabel('Номер наблюдения')
        plt.ylabel('Случайная величина')
        plt.title('Сравнение интерполяционного и сглаживающих сплайнов\n(Вариант 7: N=1670, M=1.04, σ=3.74)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def analyze_p0_behavior(self):
        """Анализ поведения при p=0"""
        if 0 not in self.splines:
            print("Сначала создайте сплайн с p=0")
            return
        
        spline_p0 = self.splines[0]
        
        # Проверка интерполяционных свойств
        y_pred_p0 = spline_p0(self.x)
        interpolation_error = np.max(np.abs(self.y - y_pred_p0))
        
        print("\n" + "="*50)
        print("АНАЛИЗ ПОВЕДЕНИЯ ПРИ p = 0")
        print("="*50)
        print(f"Максимальная ошибка интерполяции: {interpolation_error:.2e}")
        print(f"Сплайн точно проходит через точки: {interpolation_error < 1e-10}")
        
        # Анализ производных
        x_test = np.linspace(self.x.min(), self.x.max(), 10)
        derivatives_p0 = spline_p0.derivative()(x_test)
        second_derivatives_p0 = spline_p0.derivative(n=2)(x_test)
        
        print(f"Среднее значение первой производной: {np.mean(np.abs(derivatives_p0)):.4f}")
        print(f"Среднее значение второй производной: {np.mean(np.abs(second_derivatives_p0)):.4f}")
        
        # Визуализация локального поведения
        self.plot_local_behavior()
    
    def plot_local_behavior(self, region_size=100):
        """Визуализация локального поведения сплайнов"""
        # Выбираем случайную область для детального рассмотрения
        start_idx = np.random.randint(0, self.N - region_size)
        x_local = self.x[start_idx:start_idx + region_size]
        y_local = self.y[start_idx:start_idx + region_size]
        
        x_dense_local = np.linspace(x_local.min(), x_local.max(), 500)
        
        plt.figure(figsize=(15, 5))
        
        plt.subplot(1, 2, 1)
        plt.scatter(x_local, y_local, alpha=0.7, label='Данные', color='black', s=30)
        for p, spline in self.splines.items():
            y_spline_local = spline(x_dense_local)
            plt.plot(x_dense_local, y_spline_local, linewidth=2, label=f'p = {p}')
        plt.title('Локальное поведение сплайнов')
        plt.xlabel('Номер наблюдения')
        plt.ylabel('Случайная величина')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        # Сравнение производных
        for p, spline in self.splines.items():
            derivative = spline.derivative()(x_dense_local)
            plt.plot(x_dense_local, derivative, linewidth=2, label=f'p = {p}')
        plt.title('Первые производные сплайнов')
        plt.xlabel('Номер наблюдения')
        plt.ylabel('Производная')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    def calculate_smoothing_metrics(self):
        """Вычисление метрик сглаживания"""
        print("\n" + "="*50)
        print("МЕТРИКИ СГЛАЖИВАНИЯ")
        print("="*50)
        
        x_dense = np.linspace(self.x.min(), self.x.max(), 1000)
        
        metrics = []
        for p, spline in self.splines.items():
            y_spline = spline(x_dense)
            derivative = spline.derivative()(x_dense)
            
            # Метрики
            roughness = np.mean(derivative**2)  # "Шероховатость"
            deviation = np.mean((y_spline - np.mean(y_spline))**2)  # Отклонение от среднего
            
            metrics.append({
                'p': p,
                'roughness': roughness,
                'deviation': deviation,
                'smoothing_param': spline._data[-1]
            })
            
            print(f"p = {p}:")
            print(f"  Параметр сглаживания s = {spline._data[-1]:.2f}")
            print(f"  'Шероховатость' (средний квадрат производной) = {roughness:.4f}")
            print(f"  Отклонение от среднего = {deviation:.4f}")
            print()
        
        return metrics

def main():
    
    print("ПРАКТИЧЕСКОЕ ЗАДАНИЕ №3 - ВАРИАНТ 7")
    print("="*60)
    
    
    analyzer = SmoothingSplineAnalysis(N=1670, M=1.04, sigma=3.74)
    
    p_values = [0, 0.4, 0.8, 10000.0]
    analyzer.create_splines(p_values)
    
 
    analyzer.plot_splines_comparison()
    
    analyzer.analyze_p0_behavior()
    
    metrics = analyzer.calculate_smoothing_metrics()
 
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 2, 1)
    roughness_vals = [m['roughness'] for m in metrics]
    plt.plot(p_values, roughness_vals, 'o-', linewidth=2, markersize=8)
    plt.xlabel('Параметр сглаживания p')
    plt.ylabel('"Шероховатость"')
    plt.title('Зависимость шероховатости от p')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    deviation_vals = [m['deviation'] for m in metrics]
    plt.plot(p_values, deviation_vals, 'o-', linewidth=2, markersize=8, color='orange')
    plt.xlabel('Параметр сглаживания p')
    plt.ylabel('Отклонение от среднего')
    plt.title('Зависимость отклонения от p')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()