#include <iostream>
#include <cmath>
#include <string>
#include <Eigen/Dense>
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;



const double e = 2.7182818284590; //euler's number
const double h = 1e-08; //step size


struct finidif { //finite difference equation for 5-point derivative approximations
				
	VectorXd coeffs; //finite difference coefficients
	VectorXd order; //0 in every row EXCEPT for d+1, which is d! (and d is the order of the derivative)
	MatrixXd samples; //sampled points
	MatrixXd inverse; //inverted samples matrix
	int d;
	double* samps;

};

class func {

public:
	func(double *samps) {
		this->expression = this->define(0);
		this->x = 0;
		this->approx.samps = samps;

	}
	double define(double); //defines the function
	void finitedx(int, int); //generates the coefficients for a finite difference approximation of the nth derivative
	double dxatpoint(double x); //evaluates the derivative of the function at the given point
	void taylor(int); //generates the taylor series
	void sample(int N); //sample points
	void print(int N);

private:
	double x;
	double expression;
	finidif approx;
	

};



int main() {

	cout << "How many points would you like to sample? A taylor polynomial of order (samples-1) will be generated." << endl;
	int N;
	cin >> N;
	double* samps;
	double point;
	samps = new double[N];
	for (int i = 0; i < N; i++) {
		cout << "Please enter a sample point: " << endl;
		cin >> point;
		samps[i] = point;
	}
	func func(samps);
	func.taylor(N-1);
	
	


	

	return 0;
}

double func::define(double x) {
	double ans = pow(e,x);
	return ans;
}

void func::print(int N) {
	cout << "Samples: " << endl;
	for (int i = 0; i < N; i++) {
		cout << this->approx.samps[i] << endl;
	}

	cout << "Sample Matrix: " << endl;
	cout << this->approx.samples << endl;
	cout << "Inverse Matrix: " << endl;
	cout << this->approx.inverse << endl;
	cout << "Order Vector (dx order of " << this->approx.d << "): " << endl;
	cout << this->approx.order << endl;
	cout << "Output Coefficients: " << endl;
	cout << this->approx.coeffs << endl;
	return;
}

	
void func::sample(int N) {
	
	this->approx.samples.resize(N, N);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double num = this->approx.samps[i];
			this->approx.samples(j, i) = pow(num, j);
		}
	}
	return;

}

double func::dxatpoint(double x) {
	cout << "Evaluating dx^" << this->approx.d << endl;
	int size = this->approx.coeffs.size();
	double sum = 0;
	for (int i = 0; i < size; i++) {
		double point = this->approx.samps[i];

		double x = this->x + (point*h);
		cout << "new_x: " << x << endl;

		double ans = this->define(x);
		cout << "Ans = " << ans << endl;
		ans *= this->approx.coeffs[i];
		cout << "Ans after coeff: " << ans << endl;

		sum += ans;

		
	}

	cout << "total: " << sum << endl;
	double check = pow(h, this->approx.d);
	int dif = log10(abs(sum)) - log10(check);
	cout << "dif: " << dif << endl;
	cout << "horder: " << abs(log10(check)) << endl;
	if (dif >= abs(log10(h))) {
		return 0;
	}

	sum /= pow(h, this->approx.d);
	cout << "sum after division: " << sum << endl;
	return sum;
	
}

void func::finitedx(int N, int d) {
	this->sample(N);
	this->approx.inverse = this->approx.samples.inverse();

	VectorXd ord(N);
	for (int i = 0; i < N; i++) {
		ord(i) = 0;
	}

	ord(d) = tgamma(d + 1);
	
	this->approx.coeffs = this->approx.inverse * ord;
	for (int i = 0; i < N; i++) {
		double abso = abs(this->approx.coeffs(i));
		int coorder;
		if (abso == 1) {
			coorder = 0;
		}
		else {
			coorder = log10(abso);
		}
		
		int horder = log10(h);
		cout << "CoOrder: " <<coorder << endl;
		cout << horder << endl;
		if (coorder <= horder) {
			this->approx.coeffs(i) = 0;
		}
	}

	this->approx.d = d;
	this->approx.order = ord;
	cout << "Info about the order " << d << "approximation matrices" << endl;
	this->print(5);
	return;


}

void func::taylor(int terms) {
	double x = this->x;
	double* taylorcos;
	taylorcos = new double[terms];
	taylorcos[0] = this->define(x);
	for (int i = 1; i < terms; i++) {
		this->finitedx(5, i);
		taylorcos[i] = dxatpoint(x);
		taylorcos[i] /= tgamma(i + 1);
	}
	cout << "The taylor expansion of f(x) around x = " << x << " is:" << endl;
	cout << taylorcos[0] << " + ";
	for (int i = 1; i < terms-1; i++) {
		cout << taylorcos[i] << "(x-" << x << ")^" << i << " + ";
	}
	cout << taylorcos[terms - 1] << "x^" << terms-1 << endl;
	return;

}

