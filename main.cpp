#include "integrator.h"
#include "custom.h"
#include <fstream>
#include <string>
#include <windows.h>

static void SaveCSV(const TMatrix& result_matrix, const std::string& file_path)
{
    std::ofstream out(file_path);
    out.setf(std::ios::scientific);
    out.precision(20);

    for (int row = 0; row < result_matrix.rowCount(); ++row) {
        for (int col = 0; col < result_matrix.colCount(); ++col) {
            out << (double)result_matrix(row, col);
            if (col + 1 != result_matrix.colCount()) out << ",";
        }
        out << "\n";
    }
}

int main()
{
    SetConsoleOutputCP(CP_UTF8);
    // Создание объекта интегратора
    TIntegrator* integrator = new TDormandPrinceIntegrator();
    // Установка максимально допустимой ошибки интегрирования на каждом шаге
    integrator->setRelTol(1e-16L);

    // Интервал выдачи результатов при этом установить равным 0.01 с
    const long double output_step_fine = 0.01L;

    // Создание объектов моделей
    TModel* model_1 = new TArenstorfModel(0.0L, 80.0L, output_step_fine);
    TModel* model_2 = new TArenstorfModel2(0.0L, 55.0L, output_step_fine);

    // Запуск интегрирования системы ДУ модели №1
    integrator->Run(model_1);
    // Копирование результатов из объекта модели №1 в отдельную матрицу
    TMatrix result = model_1->getResult();
    // Запись результатов в файл или построение графика
    SaveCSV(result, "arenstorf_1_dt001.csv");

    // Запуск интегрирования системы ДУ модели №2
    integrator->Run(model_2);
    // Копирование результатов из объекта модели №2 в отдельную матрицу
    result = model_2->getResult();
    SaveCSV(result, "arenstorf_2_dt001.csv");

    // Для проверки корректности реализации алгоритма плотной выдачи:
    // дополнительно построить графики с интервалом выдачи 1 с

    const long double output_step_coarse = 1.0L;

    model_1 = new TArenstorfModel(0.0L, 80.0L, output_step_coarse);
    model_2 = new TArenstorfModel2(0.0L, 55.0L, output_step_coarse);

    integrator->Run(model_1);
    SaveCSV(model_1->getResult(), "arenstorf_1_dt1.csv");

    integrator->Run(model_2);
    SaveCSV(model_2->getResult(), "arenstorf_2_dt1.csv");

    return 0;
}
