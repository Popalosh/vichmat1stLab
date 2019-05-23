#include <iostream>
#include <math.h>
#include "lib/quanc8.h"
#include "lib/quanc8.cpp"
#include "lib/SPLINES.H"
#include "lib/SPLINES.CPP"

double *x = new double[11]; // Xk
double *z = new double[11]; // f(Xk)
double *y = new double[10]; // Yk
double m = -1.0; // степень для подитегрального выражения

double breakPoint = 1.0/(sin(0.6)); // вычислял аналитически

double f(double x) { // метод для вычисления точного значения
	return 1 - exp(-x);
}

double q_fun(double x) { // подинтегральное выражение для quanc8
	return pow(abs(sin(x)-(0.6)),m);
}

void quanc8Fun(double delta) {
	//----------------------------------
	m = -1.0;
	double abserr = 0.0;
    double relerr = 1.0e-7;
    double *errest= new double;
    double *flag = new double;
    int *nofun = new int;
    double *result = new double[2];
	double *intermediateResult = new double[2];
	//----------------------------------
	printf("\n Quanc8 results in main points with delta = %.10f\n\n", delta);

	for (int i = 0; i < 2; ++i) { // две итерации для разных m
	
		for (int j = 0; j < 2; ++j) { // две итерации для двух промежутков, ибо имеется разрыв
			if (j == 0) {
				double a = 0.5; //нижний предел интеграла
				double b = breakPoint - delta; //верхний предел интеграла
				quanc8(q_fun,a,b,abserr,relerr,&intermediateResult[j],errest,nofun,flag);

			} else {
				double a = breakPoint + delta; //нижний предел интеграла
				double b = 1; //верхний предел интеграла
				quanc8(q_fun,a,b,abserr,relerr,&intermediateResult[j],errest,nofun,flag);
			} 
		}
		result[i] = intermediateResult[0] + intermediateResult[1];

		printf(" m = %.1f\t", m);
		printf("Y = %.8f\n", result[i]);

		m = -0.5;
	}
	printf("\n");
}

void values () { // заполнение массивов входных значений
	for (int i = 0; i <= 10; ++i) {
		x[i] = 0.3*i;
		z[i] = f(x[i]);
		if (i != 10) {
			y[i] = 0.15 + x[i];
		}
	}
}

void lagrangePolynom(double *result) { // вычисление значений с помощью полинома Лагранжа 10й степени
	for (int i = 0; i <= 9; i++) {
		double sum = 0.0;
		for (int j = 0; j <= 10; j++) {
			double intermediate = 1.0;
			for (int k = 0; k <= 10; k++) {
				if ( j != k)
				intermediate *= (y[i] - x[k])/(x[j] - x[k]);
			}
			sum += (intermediate * f(x[j]));
		}
		result[i] = sum;
	}
}

void splineFunction(double *result) { // вычисление значений с помощью сплайн функции
	for (int i = 0; i <= 9; ++i) {
		double *b = new double[12];
		double *c = new double[12];
		double *d = new double[12];
		spline(11,x-1,z-1,b,c,d);

		result[i] = seval(11,&y[i],x-1,z-1,b,c,d);
	}
}

int main() {
	values();
	quanc8Fun(0.1);
	quanc8Fun(0.01);
	quanc8Fun(0.001);
	quanc8Fun(0.0001);
	quanc8Fun(0.00001);
	quanc8Fun(0.000001);
	quanc8Fun(0.0000001);
	quanc8Fun(0.00000001);
	quanc8Fun(0.000000001);
	quanc8Fun(0.0000000001);

	double *lagrangeResult = new double[10];
	double *splineResult = new double[10];

	lagrangePolynom(lagrangeResult);
	splineFunction(splineResult);

	printf("\n");
	printf(" X\t Exact\t\t Lagrange\t Difference\t Spline\t\t Difference\t");
	printf("\n\n");
	for (int i = 0; i <= 9; i++) {
		printf(" %.2f\t", y[i]);
		printf(" %.11f\t", f(y[i]));
		printf(" %.11f\t", lagrangeResult[i]);
		printf(" %.11f\t", abs(f(y[i]) - lagrangeResult[i]));
		printf(" %.11f\t", splineResult[i]);
		printf(" %.11f\n", abs(f(y[i]) - splineResult[i]));
	}
	return 0;
}