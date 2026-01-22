//---------------------------------------------------------------------------

#include <cmath>
#include <algorithm>
#include "integrator.h"
#include "model.h"

using std::max;
using std::min;

//===========================================================================
// class TDormandPrinceIntegrator

// Коэффициенты таблицы Дормана-Принса 5(4) :contentReference[oaicite:2]{index=2}
const long double TDormandPrinceIntegrator::c[7] =
{ 0.L, 1.L/5, 3.L/10, 4.L/5, 8.L/9, 1.L, 1.L };

const long double TDormandPrinceIntegrator::a[7][6] = {
    { 0.L },
    { 1.L/5 },
    { 3.L/40, 9.L/40 },
    { 44.L/45, -56.L/15, 32.L/9 },
    { 19372.L/6561, -25360.L/2187, 64448.L/6561, -212.L/729 },
    { 9017.L/3168, -355.L/33, 46732.L/5247, 49.L/176, -5103.L/18656 },
    { 35.L/384, 0.L, 500.L/1113, 125.L/192, -2187.L/6784, 11.L/84 }
};

const long double TDormandPrinceIntegrator::b4[7] =
{ 35.L/384, 0.L, 500.L/1113, 125.L/192, -2187.L/6784, 11.L/84, 0.L };

const long double TDormandPrinceIntegrator::b5[7] =
{ 5179.L/57600, 0.L, 7571.L/16695, 393.L/640,
 -92097.L/339200, 187.L/2100, 1.L/40 };

//---------------------------------------------------------------------------

TDormandPrinceIntegrator::TDormandPrinceIntegrator()
    : TIntegrator()
{
    // Определение ошибки округления (машинный эпсилон) :contentReference[oaicite:3]{index=3}
    long double v = 1.L;
    do {
        u = v;
        v = v / 2.L;
    } while (1.L + v > 1.L);
}

//---------------------------------------------------------------------------
// Плотная выдача: коэффициенты b_j(theta)

static inline long double b1_theta(long double theta) {
    return theta * (1.L + theta * (-1337.L/480 + theta * (1039.L/360 + theta * (-1163.L/1152))));
}
static inline long double b2_theta(long double /*theta*/) {
    return 0.L;
}
static inline long double b3_theta(long double theta) {
    return (100.L * theta * theta *
           (1054.L/9275 + theta * (-4682.L/27825 + theta * (379.L/5565))))
           / 3.L;
}
static inline long double b4_theta(long double theta) {
    return (-5.L * theta * theta *
           (27.L/40 + theta * (-9.L/5 + theta * (83.L/96))))
           / 2.L;
}
static inline long double b5_theta(long double theta) {
    return (18225.L * theta * theta *
           (-3.L/250 + theta * (22.L/375 + theta * (-37.L/600))))
           / 848.L;
}
static inline long double b6_theta(long double theta) {
    return (-22.L * theta * theta *
           (-3.L/10 + theta * (29.L/30 + theta * (-17.L/24))))
           / 7.L;
}

//---------------------------------------------------------------------------

void TDormandPrinceIntegrator::Run(TModel* Model)
{
    // Переключатель: считать ? с множителем h (как в методичке) или без h (как ты просил)
    // В методичке на стр. 3 формула ? содержит h :contentReference[oaicite:5]{index=5}
    constexpr bool k_use_h_in_error = false; // <-- поставь true, если нужно строго как в методичке

    long double t = Model->t0();             // текущее время шага интегрирования
    const long double t_end = Model->t1();   // конечное время интегрирования
    long double t_out = t;                  // время выдачи результата (равномерная сетка)

    long double h_done = 0.L;               // фактически выполненный шаг
    long double h_new = Model->step();      // предлагаемый шаг (стартуем с шага выдачи)
    long double eps_step = 0.L;             // относительная вычислительная ошибка ?

    TVector x = Model->X0();                // состояние на начале шага
    TVector x4((int)x.size());              // конец шага (4-й порядок)  x_hat
    TVector x5((int)x.size());              // конец шага (5-й порядок)  x
    TVector x_out((int)x.size());           // состояние для плотной выдачи

    const int system_dim = (int)x.size();

    Model->clearResult();
    Model->prepareResult();

    // Подготовим K
    for (int stage = 0; stage < 7; ++stage)
        K[stage].resize(system_dim);

    // Запишем начальную точку
    Model->addResult(x, t_out);
    t_out += Model->step();

    while (t < t_end)
    {
        // 1) Выбор длины шага (не перелететь конец интервала)
        h_done = h_new;
        if (t + h_done > t_end) h_done = t_end - t;

        // 2) 7 стадий метода Рунге-Кутты (Дорман-Принс) :contentReference[oaicite:6]{index=6}
        K[0] = Model->getRight(x, t);

        K[1] = Model->getRight(
            x + (h_done * a[1][0]) * K[0],
            t + c[1] * h_done
        );

        K[2] = Model->getRight(
            x + (h_done * (a[2][0] * K[0] + a[2][1] * K[1])),
            t + c[2] * h_done
        );

        K[3] = Model->getRight(
            x + (h_done * (a[3][0] * K[0] + a[3][1] * K[1] + a[3][2] * K[2])),
            t + c[3] * h_done
        );

        K[4] = Model->getRight(
            x + (h_done * (a[4][0] * K[0] + a[4][1] * K[1] + a[4][2] * K[2] + a[4][3] * K[3])),
            t + c[4] * h_done
        );

        K[5] = Model->getRight(
            x + (h_done * (a[5][0] * K[0] + a[5][1] * K[1] + a[5][2] * K[2] + a[5][3] * K[3] + a[5][4] * K[4])),
            t + c[5] * h_done
        );

        K[6] = Model->getRight(
            x + (h_done * (a[6][0] * K[0] + a[6][1] * K[1] + a[6][2] * K[2] + a[6][3] * K[3] + a[6][4] * K[4] + a[6][5] * K[5])),
            t + c[6] * h_done
        );

        // 3) Решения 4-го и 5-го порядка на конце шага :contentReference[oaicite:7]{index=7}
        x4 = x + h_done * (b4[0] * K[0] + b4[1] * K[1] + b4[2] * K[2] +
                           b4[3] * K[3] + b4[4] * K[4] + b4[5] * K[5] + b4[6] * K[6]);

        x5 = x + h_done * (b5[0] * K[0] + b5[1] * K[1] + b5[2] * K[2] +
                           b5[3] * K[3] + b5[4] * K[4] + b5[5] * K[5] + b5[6] * K[6]);

        // 4) Оценка относительной вычислительной ошибки ?
        const long double roundoff_term = 2.L * u / this->e_max;

        long double sum_sq = 0.L;
        for (int i = 0; i < system_dim; ++i)
        {
            const long double x0_i    = x[(size_t)i];   // x0i (начало шага)
            const long double x_hat_i = x4[(size_t)i];  // x?i (4-й порядок)
            const long double x_i     = x5[(size_t)i];  // xi (5-й порядок)

            const long double denom_i = max(
                1e-5L,
                max(fabsl(x_hat_i),
                    max(fabsl(x0_i), roundoff_term))
            );

            long double term = (x_hat_i - x_i) / denom_i;

            // В методичке term = h*(x_hat_i - x_i)/denom_i :contentReference[oaicite:10]{index=10}
            if constexpr (k_use_h_in_error)
                term *= h_done;

            sum_sq += term * term;
        }

        eps_step = sqrtl(sum_sq / (long double)system_dim);

        // 5) Подбор нового шага h_new :contentReference[oaicite:11]{index=11}
        // h_new = h / max(0.1, min(5, (eps/eps_max)^(1/5) / 0.9))
        const long double ratio = (eps_step <= 0.L) ? 0.L : powl(eps_step / this->e_max, 1.L/5);
        h_new = h_done / max(0.1L, min(5.L, ratio / 0.9L));

        // Если шаг плохой — отклоняем и пробуем заново с меньшим h_new
        if (eps_step > this->e_max)
            continue;

        // 6) Плотная выдача на отрезке [t, t + h_done] :contentReference[oaicite:12]{index=12}
        // theta = (t_out - t)/h;  x(t_out) ? x(t) + h * ? b_j(theta) * k_j  :contentReference[oaicite:13]{index=13}
        while ((t_out < t + h_done) && (t_out <= t_end))
        {
            const long double theta = (t_out - t) / h_done;

            const long double b1 = b1_theta(theta);
            const long double b2 = b2_theta(theta);
            const long double b3 = b3_theta(theta);
            const long double b4v = b4_theta(theta);
            const long double b5v = b5_theta(theta);
            const long double b6 = b6_theta(theta);

            x_out = x + h_done * (b1 * K[0] + b2 * K[1] + b3 * K[2] +
                                  b4v * K[3] + b5v * K[4] + b6 * K[5]);

            Model->addResult(x_out, t_out);
            t_out += Model->step();
        }

        // 7) Принятие шага: на следующий шаг берём решение 5-го порядка :contentReference[oaicite:14]{index=14}
        x = x5;
        t += h_done;
    }
}
