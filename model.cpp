//---------------------------------------------------------------------------

#include "model.h"

//---------------------------------------------------------------------------

void TModel::addResult( const TVector& X, long double t )
{
    // Проверим, выходит ли счётчик строк в матрице результатов за пределы последней строки
    // Если да, то увеличим количество строк на 1
    if (N == Result.rowCount())
        Result.resize(N + 1, getOrder() + 1);

    // Поместим результаты в последнюю строку матрицы Result
    // Момент времени помещается в 0-ой столбец, вектор состояния - в остальные столбцы
    Result(N, 0) = t;

    const int state_size = (int)X.size();
    for (int state_i = 0; state_i < state_size; ++state_i)
        Result(N, state_i + 1) = X[(size_t)state_i];

    // Увеличим N
    ++N;
}

void TModel::clearResult()
{
    // Очистим матрицу результатов и сбросим счётчик строк
    Result.resize(0, 0);
    N = 0;
}

void TModel::prepareResult()
{
    // Зададим матрице результатов такой размер, чтобы поместились все значения вектора состояния
    // и соответствующих им моментов времени на интервале [t0; t1] с шагом step
    const long double time_span = t1_ - t0_;
    const int rows_hint = (int)(time_span / step_) + 2;

    Result.resize(rows_hint, getOrder() + 1);

    // Сбросим счётчик строк в матрице результатов
    N = 0;
}
