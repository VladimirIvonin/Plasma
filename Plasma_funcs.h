#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <vector>
#include <map>

#define PI 3.14159265
#define FLENGTH 5000
#define DELTAF 10

using namespace std;

struct Plasma_pars
{
    double lambda; //Длина волны радара, м
    double Ne; //Концентрация=Ne*10^11, 1/м^3
    double Te; //Температура электронов, K
    double Ti; //Температура ионов, K
    double nu_i; //Частота столкновений ионов, Гц
    double nu_e; //Частота столкновений электронов, Гц. Для бесстолкновительной плазмы обе частоты нужно занулить.
    double Con[103]; /* Контейнер для ионного состава, где первый компонент — масса иона,
        второй — его процентное сожержание в плазме. Например, для однородой O+ плазмы: Con[16] = 100
        Для 30% NO+ и 70% O+: Con[29] = 30, Con[16] = 70. Размер контейнера определяет кол-во сортов ионов. 
        Молярные массы сортов частиц, обитающих в ионосфере Земли:
        N: 14
        N2: 28
        NO: 29
        O2: 32
        O: 16
        H: 1
        He: 4
    */
};

complex<double> cintegral(complex<double> a, complex<double> b, unsigned step_count);
void Spectrum(double *x, double *y, Plasma_pars P);
void Set_pars(char* file, Plasma_pars &P);
