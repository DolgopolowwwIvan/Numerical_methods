#include "complex_operations.h"
#include <cmath>
#include <algorithm>

using namespace std;

vector<Complex> dft(const vector<Complex>& input) {
    int N = input.size();
    vector<Complex> output(N, 0);

    for (int m = 0; m < N; m++) {
        for (int n = 0; n < N; n++) {
            double angle = -2 * PI * m * n / N;
            output[m] += input[n] * exp(Complex(0, angle));
        }
    }
    return output;
}

vector<Complex> idft(const vector<Complex>& input) {
    int N = input.size();
    vector<Complex> output(N, 0);

    for (int n = 0; n < N; n++) {
        for (int m = 0; m < N; m++) {
            double angle = 2 * PI * m * n / N;
            output[n] += input[m] * exp(Complex(0, angle));
        }
        output[n] /= N;
    }
    return output;
}

vector<Complex> fft(const vector<Complex>& input) {
    int N = input.size();

    // 1. Копируем и делаем бит-реверс
    vector<Complex> result = input;
    for (int i = 1, j = 0; i < N; i++) {
        int bit = N >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j) swap(result[i], result[j]);
    }

    // 2. Основной цикл
    for (int len = 2; len <= N; len <<= 1) {
        int M = len / 2;

        vector<Complex> twiddles(M);
        for (int m = 0; m < M; m++) {
            double angle = -2 * PI * m / len;
            twiddles[m] = exp(Complex(0, angle));
        }

        // Обрабатываем блоки длины len
        for (int i = 0; i < N; i += len) {
            for (int m = 0; m < M; m++) {
                // Четные в текущем блоке
                Complex u = result[i + m];     
                // Нечетные в текущем блоке   
                Complex v = result[i + m + M];    
                Complex twiddle = twiddles[m];

                result[i + m] = u + twiddle * v;
                result[i + m + M] = u - twiddle * v;
            }
        }
    }

    return result;
}

vector<Complex> ifft(const vector<Complex>& input) {
    int N = input.size();

    // z(-j) = z(N-j) - используем свойство периодичности
    vector<Complex> conjugated_input(N);
    for (int j = 0; j < N; j++) {
        conjugated_input[j] = conj(input[j]);
    }

    // Вычисляем FFT от сопряженного
    vector<Complex> temp = fft(conjugated_input);

    // Сопрягаем результат и делим на N
    vector<Complex> output(N);
    for (int j = 0; j < N; j++) {
        output[j] = conj(temp[j]) / double(N);
    }

    return output;
}