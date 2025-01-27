#include "Plasma_funcs.h"

complex<double> cintegral(complex<double> a, complex<double> b, unsigned step_count) 
{
	complex<double> sum=complex<double> (0.0, 0.0);
	complex<double> step;
	if(step_count==0) return sum;

	step=(b-a)/double(step_count);
	for(unsigned i=1; i<step_count; ++i )
        sum+=exp(pow(a+double(i)*step,2.0));

	sum+=(exp(a*a)+exp(b*b))/2.0;
	sum*=step;
	return sum;
}

void Spectrum(double *x, double *y, Plasma_pars P)
{
	for(size_t i=0; i<2*FLENGTH; i++)
	{
		x[i]=0.0;
		y[i]=0.0;
	}
	complex<double> im=complex<double> (0.0, 1.0);
	complex<double> izero=complex<double> (0.0, 0.0);
	vector<double> Y(2*FLENGTH), S(2*FLENGTH);
	double K=1.380648813;//Boltzmann constant (*pow(10,-23))
	double E0=8.854187817;//dielectric constant (*pow(10,-12))
	double eMass=9.10938291;//mass of electron(kg) (*pow(10,-31))
	double q=1.60217657;//charge (*pow(10,-19))
	double S_max=0.0;
	double k=2.0*PI/P.lambda;
	double alpha, D, R;
	complex<double> Rwe, Iwe, Rwi, Iwi, Xe, Xi, Ae, Ai, De, Di, Ce, Ci, SumE, SumI, epsilon;
	double a, b, w, f, Be, Bi;
	R=P.Te/P.Ti;
	k=4.0*PI/P.lambda;

    for(size_t M=0; M<103; M++)
    {
		if(P.Con[M]!=0.0)
		{
			double Cm=P.Con[M];
			double iMass=1.67262177774*double(M);//mass of ion(kg) (*pow(10,-27))
			a=sqrt(2.0*K*P.Te/eMass)*pow(10.0,4.0);
			b=sqrt(2.0*K*P.Ti/iMass)*pow(10.0,2.0);
			D=sqrt(K*E0*100.0/(P.Ne*Cm*q*q)*(P.Te*P.Ti/(P.Te+P.Ti)))*pow(10.0,-4.0);
			alpha=1.0/(D*k);
			for(int i=0; i<2*FLENGTH; i++)
			{
				f=double(i-FLENGTH)*double(DELTAF);
				w=2.0*PI*f;
				Xe=(w-im*P.nu_e)/(k*a);
				Xi=(w-im*P.nu_i)/(k*b);
				SumE=cintegral(izero, Xe, 1000);
				SumI=cintegral(izero, Xi, 1000);
				if(!isnormal(real(SumE)))
					SumE=complex<double> (0.0, 0.0);
				if(!isnormal(real(SumI)))
					SumI=complex<double> (0.0, 0.0);
					
				Rwe=1.0-2.0*Xe*exp(-Xe*Xe)*SumE;
				Iwe=sqrt(PI)*Xe*exp(-Xe*Xe);
				Rwi=1.0-2.0*Xi*exp(-Xi*Xi)*SumI;
				Iwi=sqrt(PI)*Xi*exp(-Xi*Xi);

				if(P.nu_e==0.0 or P.nu_i==0.0)
				{
					Ae=complex<double> (1.0+alpha*alpha*R*real(Rwi), alpha*alpha*R*real(Iwi));
					Ai=complex<double> (alpha*alpha*real(Rwe), alpha*alpha*real(Iwe));
					epsilon=complex<double> (1.0+alpha*alpha*real(Rwe+R*Rwi), alpha*alpha*real(Iwe+R*Iwi));
					Y[i]=2.0*sqrt(PI)/(k*a)*(sqrt(iMass*P.Te/(eMass*P.Ti))*real(exp(-Xi*Xi))*norm(Ai))/norm(epsilon);
				}
				else
				{
					De=(im*P.nu_e)/(k*a)*(2.0*exp(-Xe*Xe)*SumE+im*sqrt(PI)*exp(-Xe*Xe));
					Di=(im*P.nu_i)/(k*b)*(2.0*exp(-Xi*Xi)*SumI+im*sqrt(PI)*exp(-Xi*Xi));

					Be=(1.0/(k*a*norm(1.0+De)))*imag(2.0*exp(-Xe*Xe)*SumE+im*sqrt(PI)*exp(-Xe*Xe))\
						-norm(De)/(P.nu_e*norm(1.0+De));
					Bi=(1.0/(k*b*norm(1.0+Di)))*imag(2.0*exp(-Xi*Xi)*SumI+im*sqrt(PI)*exp(-Xi*Xi))\
						-norm(Di)/(P.nu_i*norm(1.0+Di));

					Ce=alpha*alpha*(Rwe-im*Iwe)/(1.0+De);
					Ci=R*alpha*alpha*(Rwi-im*Iwi)/(1.0+Di);
					epsilon=1.0+Ce+Ci;
					// cout<<Rwe+im*Iwe<<"\t"<<Rwi+im*Iwi<<"\t";
					// cout<<De<<"\t"<<Di<<"\t"<<Be<<"\t"<<Bi<<"\t"<<Ce<<"\t"<<Ci<<"\n";
					Y[i]=norm(Ce/(1.0+Ce+Ci))*Bi;   
				}
				// Y[i]=2.0*norm((1.0+Ci)/epsilon)*Be+2.0*norm(Ce/epsilon)*Bi;
			}
			for(int i=0; i<2*FLENGTH; i++)
			{
				if(isfinite(Y[i]))
					S[i]+=Y[i];
				if(S[i]>S_max)
					S_max=S[i];
				x[i]=f=double(i-FLENGTH)*double(DELTAF);
			}
		}
    }
	for(int i=0; i<2*FLENGTH; i++)
		y[i]=S[i]/S_max;
}

void Set_pars(char* file, Plasma_pars &P)
{
    ifstream inp;
	for(size_t i=0; i<103; i++)
		P.Con[i]=0.0;
    string line;
    int found;

    inp.open(file);
    printf("%s\n", file);
	getline(inp, line);
	getline(inp, line);
	found=line.find_first_of(" = ");
	P.lambda=stod(line.substr(found+3, line.length()-found-3));
    printf("lambda %f\n", P.lambda);
	getline(inp, line);
    found=line.find_first_of(" = ");
	P.Ne=stod(line.substr(found+3, line.length()-found-3));
	printf("Ne %f\n", P.Ne);
    getline(inp, line);
    found=line.find_first_of(" = ");
	P.Te=stod(line.substr(found+3, line.length()-found-3));
	printf("Te %f\n", P.Te);
    getline(inp, line);
    found=line.find_first_of(" = ");
	P.Ti=stod(line.substr(found+3, line.length()-found-3));
	printf("Ti %f\n", P.Ti);
    getline(inp, line);
    found=line.find_first_of(" = ");
	P.nu_i=stod(line.substr(found+3, line.length()-found-3));
	printf("nu_i %f\n", P.nu_i);
    getline(inp, line);
    found=line.find_first_of(" = ");
	P.nu_e=stod(line.substr(found+3, line.length()-found-3));
	printf("nu_e %f\n", P.nu_e);
    getline(inp, line);
    found=line.find_first_of("[");
    while(found!=-1)
    {
        double M=0.0, Cm=0.0;
        M=stod(line.substr(found+1, 3));
        found=line.find(" = ");
        Cm=stod(line.substr(found+3, 3));
        P.Con[unsigned(M)]=Cm;
        printf("Con[%3.0f] = %3.0f\n", M, Cm);
        getline(inp, line);
        found=line.find_first_of("[");
    }
}