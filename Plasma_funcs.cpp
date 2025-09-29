#include "Plasma_funcs.hpp"
#include "FFT.h"

#include <numbers>
#include <iostream>
#include <vector>
#include <fstream>


constexpr double PI = std::numbers::pi;//π
constexpr double E0 = 8.85418781762039E-12;///<Электрическая постоянная, м⁻³·кг⁻¹·с⁴·А²

constexpr size_t FLENGTH = 1250;//Половина длины выходных массивов
constexpr size_t DOUBLEFLENGTH = 2*FLENGTH;//Длина выходных массивов
constexpr unsigned short dF_Hz = 100;//Шаг по частоте для построения спектра, Гц
constexpr size_t dtau_us = 4;//Шаг по задержке для построения АКФ, мкс
constexpr size_t LENGTH = 175;//Длина массива АКФ


std::complex<double> cintegral(std::complex<double> a, std::complex<double> b, unsigned step_count)
{
    std::complex<double> sum{0.0, 0.0};
    std::complex<double> step;
    if(step_count == 0)
        return sum;

    step = (b - a)/double(step_count);
    for (unsigned i = 1; i < step_count; ++i)
        sum += exp(pow(a + double(i)*step, 2.0));

    sum += (exp(a*a) + exp(b*b))/2.0;
    sum *= step;
	return sum;
}

void Spectrum(double *freqs_Hz, double *Amps, Plasma_pars *plasma_pars)
{
    for(size_t i=0; i<DOUBLEFLENGTH; i++)
	{
        freqs_Hz[i] = 0.0;
        Amps[i] = 0.0;
	}
    std::complex<double> im{0.0, 1.0};//Мнимая единица
    std::complex<double> izero{0.0, 0.0};//Комплексный ноль
    std::vector<double> Y(DOUBLEFLENGTH), S(DOUBLEFLENGTH);
    double K = 1.380648813;//Boltzmann constant (*pow(10,-23))
    double E0 = 8.85418781762039;//dielectric constant (*pow(10,-12))
    double eMass = 9.10938291;//mass of electron(kg) (*pow(10,-31))
    double q = 1.602176565;//charge (*pow(10,-19))
    double S_max = 0.0;
    //double k = 2.0*PI/plasma_pars->wavelength_m;
	double alpha, D, R;
    std::complex<double> Rwe, Iwe, Rwi, Iwi, Xe, Xi, Ae, Ai, De, Di, Ce, Ci, SumE, SumI, epsilon;
	double a, b, w, f, Be, Bi;
    R = plasma_pars->Te_K/plasma_pars->Ti_K;
    double k = 4.0*PI/plasma_pars->wavelength_m;
	
    double C = 0.0;
    for (size_t M = 0; M < 103; M++)
        C += plasma_pars->Con[M];
    if (C != 100.0)
        std::cerr<<"Total relative plasma density is more or less than 100%!"<<std::endl;

    for (size_t M = 0; M < 103; M++)
    {
        if (plasma_pars->Con[M] != 0.0)
		{
            double Cm = plasma_pars->Con[M]/C;
            double iMass = 1.67262177774*double(M);//mass of ion(kg) (*pow(10,-27))
            a = sqrt(2.0*K*plasma_pars->Te_K/eMass)*1.E4;
            b = sqrt(2.0*K*plasma_pars->Ti_K/iMass)*1.E2;
            D = sqrt(K*E0*100.0/(plasma_pars->Ne_105cm3*Cm*q*q)*
                (plasma_pars->Te_K*plasma_pars->Ti_K/(plasma_pars->Te_K + plasma_pars->Ti_K)))*1.E-4;
            alpha = 1.0/(D*k);
            for (int i = 0; i < DOUBLEFLENGTH; i++)
			{
                f = double(i - static_cast<int>(FLENGTH))*double(dF_Hz);
                if (fabs(f) < 1.E4)
				{
                    w = 2.0*PI*f;
                    Xe = (w-im*plasma_pars->nu_e_Hz)/(k*a);
                    Xi = (w-im*plasma_pars->nu_i_Hz)/(k*b);
                    SumE = cintegral(izero, Xe, 1000);
                    SumI = cintegral(izero, Xi, 1000);
                    if (!std::isnormal(real(SumE)))
                        SumE = izero;
                    if (!std::isnormal(real(SumI)))
                        SumI = izero;
						
                    Rwe = 1.0 - 2.0*Xe*exp(-Xe*Xe)*SumE;
                    Iwe = sqrt(PI)*Xe*exp(-Xe*Xe);
                    Rwi = 1.0 - 2.0*Xi*exp(-Xi*Xi)*SumI;
                    Iwi = sqrt(PI)*Xi*exp(-Xi*Xi);

                    if (plasma_pars->nu_e_Hz == 0.0 or plasma_pars->nu_i_Hz == 0.0)
					{
                        Ae = std::complex<double> (1.0 + alpha*alpha*R*real(Rwi), alpha*alpha*R*real(Iwi));
                        Ai = std::complex<double> (alpha*alpha*real(Rwe), alpha*alpha*real(Iwe));
                        epsilon = std::complex<double> (1.0 + alpha*alpha*real(Rwe+R*Rwi),
                            alpha*alpha*real(Iwe + R*Iwi));
                        Y[i] = 2.0*sqrt(PI)/(k*a)*(sqrt(iMass*plasma_pars->Te_K/(eMass*plasma_pars->Ti_K))*
                            real(exp(-Xi*Xi))*norm(Ai))/norm(epsilon);
					}
					else
					{
                        De = (im*plasma_pars->nu_e_Hz)/(k*a)*(2.0*exp(-Xe*Xe)*SumE+im*sqrt(PI)*exp(-Xe*Xe));
                        Di = (im*plasma_pars->nu_i_Hz)/(k*b)*(2.0*exp(-Xi*Xi)*SumI+im*sqrt(PI)*exp(-Xi*Xi));

                        // Be = (1.0/(k*a*norm(1.0+De)))*imag(2.0*exp(-Xe*Xe)*SumE + im*sqrt(PI)*exp(-Xe*Xe))\
                        //     -norm(De)/(plasma_pars->nu_e_Hz*norm(1.0+De));
                        Bi = (1.0/(k*b*norm(1.0+Di)))*imag(2.0*exp(-Xi*Xi)*SumI + im*sqrt(PI)*exp(-Xi*Xi))\
                            -norm(Di)/(plasma_pars->nu_i_Hz*norm(1.0+Di));

                        Ce = alpha*alpha*(Rwe - im*Iwe)/(1.0 + De);
                        Ci = R*alpha*alpha*(Rwi - im*Iwi)/(1.0 + Di);
                        epsilon = 1.0 + Ce + Ci;
						// cout<<Rwe+im*Iwe<<"\t"<<Rwi+im*Iwi<<"\t";
						// cout<<De<<"\t"<<Di<<"\t"<<Be<<"\t"<<Bi<<"\t"<<Ce<<"\t"<<Ci<<"\n";
                        Y[i] = norm(Ce/(1.0 + Ce + Ci))*Bi;
					}
				}
				else
                    Y[i] = 0.0;
				// Y[i]=2.0*norm((1.0+Ci)/epsilon)*Be+2.0*norm(Ce/epsilon)*Bi;
			}
            for (int i = 0; i < DOUBLEFLENGTH; i++)
			{
                if (std::isfinite(Y[i]))
                    S[i] += Y[i];
                if(S[i] > S_max)
                    S_max = S[i];
                freqs_Hz[i] = double(i - static_cast<int>(FLENGTH))*double(dF_Hz);
			}
		}
    }
    for (int i = 0; i < DOUBLEFLENGTH; i++)
        Amps[i] = S[i]/S_max;
}

void ACF(double *lags_us, double *Amps, double *spectrum)
{
    double SS[DOUBLEFLENGTH];
    ShortComplex a[DOUBLEFLENGTH];
    for (int i = 0; i < FLENGTH; i++)
    {
        a[i].re = spectrum[i + FLENGTH];
        a[i].im = 0.0;
        a[i + FLENGTH].re = spectrum[i];
        a[i + FLENGTH].im = 0.0;
    }
    universal_fft(a, DOUBLEFLENGTH, true);
    for (int i = 0; i < DOUBLEFLENGTH; i++)
        SS[i] = std::hypot(a[i].re, a[i].im);

    /*fftw_std::complex *Sp, *Cp;
	fftw_plan Plan;
	for(size_t i=0; i<LENGTH; i++)
	{
        x[i] = 0.0;
        y[i] = 0.0;
	}
    Sp = (fftw_std::complex*) fftw_malloc(sizeof(fftw_std::complex) * (2*FLENGTH));
    Cp = (fftw_std::complex*) fftw_malloc(sizeof(fftw_std::complex) * (2*FLENGTH));
    Plan = fftw_plan_dft_1d((2*FLENGTH), Sp, Cp, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (size_t w = 0; w < FLENGTH; w++)
	{
        Sp[w][0] = S[w + FLENGTH];
        Sp[w][1] = 0.0;
        Sp[w + FLENGTH][0] = S[w];
        Sp[w + FLENGTH][1] = 0.0;
	}
    fftw_execute(Plan);
    double Norm = sqrt(Cp[0][0]*Cp[0][0] + Cp[0][1]*Cp[0][1]);*/
    for(size_t tau = 0; tau < LENGTH; tau++)
    {
        lags_us[tau] = double(tau*dtau_us);
        Amps[tau] = /*pow(1.0 - double(tau)/double(LENGTH), 2.0)**/a[tau].re/SS[0];// /Norm;
    }
	
    /*fftw_free(Sp);
	fftw_free(Cp);
    fftw_destroy_plan(Plan);*/
}

void read_pars(char* fullfilename, Plasma_pars *plasma_pars)
{
    std::ifstream inp;
    for(size_t i = 0; i < 103; i++)
        plasma_pars->Con[i] = 0.0;
    std::string line;
    int found;

    inp.open(fullfilename);
    printf("%s\n", fullfilename);
	getline(inp, line);
	getline(inp, line);
    found = line.find_first_of(" = ");
    plasma_pars->wavelength_m = stod(line.substr(found + 3, line.length() - found - 3));
    printf("lambda %f\n", plasma_pars->wavelength_m);
	getline(inp, line);
    found = line.find_first_of(" = ");
    plasma_pars->Ne_105cm3 = stod(line.substr(found + 3, line.length() - found - 3));
    printf("Ne %f\n", plasma_pars->Ne_105cm3);
    getline(inp, line);
    found = line.find_first_of(" = ");
    plasma_pars->Te_K = stod(line.substr(found + 3, line.length() - found - 3));
    printf("Te %f\n", plasma_pars->Te_K);
    getline(inp, line);
    found = line.find_first_of(" = ");
    plasma_pars->Ti_K = stod(line.substr(found + 3, line.length() - found - 3));
    printf("Ti %f\n", plasma_pars->Ti_K);
    getline(inp, line);
    found = line.find_first_of(" = ");
    plasma_pars->nu_i_Hz = stod(line.substr(found + 3, line.length() - found - 3));
    printf("nu_i %f\n", plasma_pars->nu_i_Hz);
    getline(inp, line);
    found = line.find_first_of(" = ");
    plasma_pars->nu_e_Hz = stod(line.substr(found + 3, line.length() - found - 3));
    printf("nu_e %f\n", plasma_pars->nu_e_Hz);
    getline(inp, line);
    found = line.find_first_of("[");
    while(found != -1)
    {
        double M = 0.0, Cm = 0.0;
        M = std::stod(line.substr(found + 1, 3));
        found = line.find(" = ");
        Cm = std::stod(line.substr(found + 3, 3));
        plasma_pars->Con[unsigned(M)] = Cm;
        printf("Con[%3.0f] = %3.0f\n", M, Cm);
        getline(inp, line);
        found = line.find_first_of("[");
    }
}
