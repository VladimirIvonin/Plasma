#include "../Plasma.hpp"
#include "FFT.h"

#include <numbers>
#include <complex>
#include <iostream>
#include <vector>
#include <format>


namespace PlasmaSources
{
    constexpr double PI = std::numbers::pi;//π
    constexpr double E0 = 8.85418781762039E-12;///<Электрическая постоянная, м⁻³·кг⁻¹·с⁴·А²

    constexpr size_t FLENGTH = Plasma::Plasma_pars::spectrum_array_length/2;///<Половина длины выходных массивов для расчёта спектра
    constexpr unsigned short dF_Hz = 100;//Шаг по частоте для построения спектра, Гц
    constexpr size_t dtau_us = 4;//Шаг по задержке для построения АКФ, мкс

    constexpr double error = 1.E-5;//Окрестность для сравнения `double`



    std::complex<double> cintegral(std::complex<double> a, std::complex<double> b, unsigned step_count)
    {
        std::complex<double> sum{0.0, 0.0};
        std::complex<double> step;
        if (step_count == 0)
            return sum;

        step = (b - a)/double(step_count);
        for (unsigned i = 1; i < step_count; ++i)
            sum += exp(pow(a + double(i)*step, 2.0));

        sum += (exp(a*a) + exp(b*b))/2.0;
        sum *= step;
        return sum;
    }
}


void Plasma::Spectrum(double *freqs_Hz, double *normed_spectrum, Plasma_pars *plasma_pars)
{
    for (size_t i = 0; i < Plasma_pars::spectrum_array_length; i++)
	{
        freqs_Hz[i] = 0.0;
        normed_spectrum[i] = 0.0;
	}
    std::complex<double> im{0.0, 1.0};//Мнимая единица
    std::complex<double> izero{0.0, 0.0};//Комплексный ноль
    std::vector<double> Y(Plasma_pars::spectrum_array_length), S(Plasma_pars::spectrum_array_length);
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
    double k = 4.0*PlasmaSources::PI/plasma_pars->wavelength_m;
	
    double C = 0.0;
    for (size_t M = 0; M < 103; M++)
        C += plasma_pars->Con[M];
    if (C < 100.0 - 0.5 || C > 100.0 + 0.5)
    {
        std::cerr<<"Total relative plasma density is more or less than 100%!";
        std::cerr<<std::format(" ({:.{}f})%)", C, 10)<<std::endl;
    }

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
            for (int i = 0; i < Plasma_pars::spectrum_array_length; i++)
			{
                f = double(i - static_cast<int>(PlasmaSources::FLENGTH))*double(PlasmaSources::dF_Hz);
                if (fabs(f) < 1.E4)
				{
                    w = 2.0*PlasmaSources::PI*f;
                    Xe = (w - im*plasma_pars->nu_e_Hz)/(k*a);
                    Xi = (w - im*plasma_pars->nu_i_Hz)/(k*b);
                    SumE = PlasmaSources::cintegral(izero, Xe, 1000);
                    SumI = PlasmaSources::cintegral(izero, Xi, 1000);
                    if (!std::isnormal(real(SumE)))
                        SumE = izero;
                    if (!std::isnormal(real(SumI)))
                        SumI = izero;
						
                    Rwe = 1.0 - 2.0*Xe*exp(-Xe*Xe)*SumE;
                    Iwe = sqrt(PlasmaSources::PI)*Xe*exp(-Xe*Xe);
                    Rwi = 1.0 - 2.0*Xi*exp(-Xi*Xi)*SumI;
                    Iwi = sqrt(PlasmaSources::PI)*Xi*exp(-Xi*Xi);

                    if (plasma_pars->nu_e_Hz == 0.0 or plasma_pars->nu_i_Hz == 0.0)
					{
                        Ae = std::complex<double> (1.0 + alpha*alpha*R*real(Rwi), alpha*alpha*R*real(Iwi));
                        Ai = std::complex<double> (alpha*alpha*real(Rwe), alpha*alpha*real(Iwe));
                        epsilon = std::complex<double> (1.0 + alpha*alpha*real(Rwe+R*Rwi),
                            alpha*alpha*real(Iwe + R*Iwi));
                        Y[i] = 2.0*sqrt(PlasmaSources::PI)/(k*a)*(sqrt(iMass*plasma_pars->Te_K/(eMass*plasma_pars->Ti_K))*
                            real(exp(-Xi*Xi))*norm(Ai))/norm(epsilon);
					}
					else
					{
                        De = (im*plasma_pars->nu_e_Hz)/(k*a)*(2.0*exp(-Xe*Xe)*SumE+im*sqrt(PlasmaSources::PI)*exp(-Xe*Xe));
                        Di = (im*plasma_pars->nu_i_Hz)/(k*b)*(2.0*exp(-Xi*Xi)*SumI+im*sqrt(PlasmaSources::PI)*exp(-Xi*Xi));

                        // Be = (1.0/(k*a*norm(1.0+De)))*imag(2.0*exp(-Xe*Xe)*SumE + im*sqrt(PI)*exp(-Xe*Xe))\
                        //     -norm(De)/(plasma_pars->nu_e_Hz*norm(1.0+De));
                        Bi = (1.0/(k*b*norm(1.0+Di)))*imag(2.0*exp(-Xi*Xi)*SumI + im*sqrt(PlasmaSources::PI)*exp(-Xi*Xi))\
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
            for (int i = 0; i < Plasma_pars::spectrum_array_length; i++)
			{
                if (std::isfinite(Y[i]))
                    S[i] += Y[i];
                if(S[i] > S_max)
                    S_max = S[i];
                freqs_Hz[i] = double(i - static_cast<int>(PlasmaSources::FLENGTH))*double(PlasmaSources::dF_Hz);
			}
		}
    }
    for (int i = 0; i < Plasma_pars::spectrum_array_length; i++)
        normed_spectrum[i] = S[i]/S_max;
}


void Plasma::ACF(double *lags_us, double *Amps, double *spectrum, bool multiply_by_triangle)
{
    double SS[Plasma_pars::spectrum_array_length];
    PlasmaSources::ShortComplex a[Plasma_pars::spectrum_array_length];
    for (int i = 0; i < PlasmaSources::FLENGTH; i++)
    {
        a[i].re = spectrum[i + PlasmaSources::FLENGTH];
        a[i].im = 0.0;
        a[i + PlasmaSources::FLENGTH].re = spectrum[i];
        a[i + PlasmaSources::FLENGTH].im = 0.0;
    }
    PlasmaSources::universal_fft(a, Plasma_pars::spectrum_array_length, true);
    for (int i = 0; i < Plasma_pars::spectrum_array_length; i++)
        SS[i] = std::hypot(a[i].re, a[i].im);

    /*fftw_std::complex *Sp, *Cp;
	fftw_plan Plan;
    for(size_t i=0; i<acf_array_length; i++)
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
    for (size_t tau = 0; tau < Plasma_pars::acf_array_length; tau++)
    {
        lags_us[tau] = double(tau*PlasmaSources::dtau_us);
        Amps[tau] = a[tau].re/SS[0];// /Norm;
        if (multiply_by_triangle)
            // normed_spectrum[tau] *= pow(1.0 - double(tau)/double(acf_array_length), 2.0);
            Amps[tau] *= 1.0 - double(tau)/double(Plasma_pars::acf_array_length);
    }
	
    /*fftw_free(Sp);
	fftw_free(Cp);
    fftw_destroy_plan(Plan);*/
}
