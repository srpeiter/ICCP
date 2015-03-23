#ifndef MONTECARLO_H
#define MONTECARLO_H

class particle
{

friend class observable;
 // defining a friend class so to have acces to class particles


private:		// all variable needed to do computations

double x, a;
int N;

protected:

public:
particle( double x, double a, int N) :  x(x), a(a), N(N) {}

~particle () {};

double wavefunction();	// total wavefunction

void generate_metropolis(double step);



};

class observable //this class is used for calculating energy integral with montecarlo integration
// This class is a friend of class particle
{
// minimize is a friend function of class observable
friend double* minimize(observable& obj, double a0, double a1, double stop_prec);

private:
particle temp;

public:
observable(particle& temp) : temp(temp) {}

double energy();	



void metropolis_walker();

double comp_integral(double a_inp);
};



#endif
