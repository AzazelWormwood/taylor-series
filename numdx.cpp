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
	int d = 4; //derivative order
	double* samps; //sampled points, passed by reference to allow it in the constructor

};

class func {

public:
	func() { //default constructor

		this->x = 0;
	}

	func(double* samps) { //constructor with sampled points
		this->x = 0;
		this->approx.samps = samps;

	}

	double define(double); //allows for derived classes to define the function TODO: FIX INHERITANCE
		
	void finitedx(int, int); //generates the coefficients for a finite difference approximation of the nth derivative
	double dxatpoint(double); //evaluates the derivative of the function at the given point
	void taylor(int); //generates the taylor series
	void sample(int); //sample points
	void print(int); //prints info about the finite difference matrices
	double numint(double, double); //performs a numerical approximation of the integral from [a,b] inclusive, using Simpson's rule

private:
	double x;
	finidif approx;
	

};
/*

class etothex: public func{

public:

	etothex() {

	}
	etothex(double* samps) {

	}


	double define(double);
};

class pchempractice : public func {

public:
	pchempractice() {

	}
	pchempractice(double* samps) {

	}

	double define(double);
};
*/

int main() {

	cout << "How many points would you like to sample? A taylor polynomial of order (samples-1) will be generated." << endl;
	int N;
	cin >> N;
	double* samps;
	double point;
	samps = new double[N];
	for (int i = 0; i < N; i++) { //sort these in ascending order TODO
		cout << "Please enter a sample point: " << endl;
		cin >> point;
		samps[i] = point;
	}
	func func(samps);
	func.taylor(N);
	func.print(N);

	
	
	


	

	return 0;
}

/*/double pchempractice::define(double x) {
	int a = 1;
	double square = pow(x, 2);
	double eterm = pow(e, (-a * square));

	double parenth = (2 * a * square) - 1;
	parenth = pow(parenth, 2);
	double ans = parenth * eterm * square;


	return ans;
}*/


/*double etothex::define(double x) {
	double ans = pow(e, x);
	return ans;
}*/

double func::define(double x) {
	double ans = sin(x)+cos(x);
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
			this->approx.samples(j, i) = pow(num, j); //for a square matrix of size NxN, this function takes N sampled points 
		}											  //and for each point raises it to each power up to N-1 to form the columns of the sample matrix
													  //for example, taking points (-3, 2, 5, 6, 9) and derivitave of order 4, the matrix would be:
													  /*                                      | 1   1   1   1     1 |
													                                          |-3   2   5   6     9 |
														                                      | 9   4   25  36    81|    
																							  |-27  8  125 216   729|   
																							  | 81  16 625 1296 6561|                           */
	}
	return;

}


void func::finitedx(int N, int d) {
	this->sample(N);
	this->approx.inverse = this->approx.samples.inverse(); //using the equation derived by Taylor, Cameron at https://web.media.mit.edu/~crtaylor/calculator.html
														   //the inverse matrix of the above sample matrix*a vector based on the derivative order, scaled by 1/h^order
														   //using the above example, the inverse matrix is: 
														   /*                                 |  0.125   -0.122222  0.039120 -0.005093  0.000231   |
																							  |	 1.92857 -0.278571 -0.164286  0.040476 -0.002381   |
																							  |	-3.375    1.5       0.343750 -0.145833  0.010417   |
																							  |	 2.5     -1.19444  -0.231481  0.120370 -0.009259   |
																							  |	-0.178571 0.095238  0.012897 -0.009920  0.000992   |                      */
														   //and the order vector is [0,0,0,0,24] (in general, all entries are 0 except for entry d+1, which is d! 
	VectorXd ord(N);
	for (int i = 0; i < N; i++) {
		ord(i) = 0;
	}

	ord(d) = tgamma(d + 1); //gamma function used as factorial function, off by 1 so input (d+1) for correct factorial

	this->approx.coeffs = this->approx.inverse * ord;

	for (int i = 0; i < N; i++) { //fixes precision limits to avoid floating point errors; sets coefficients of order of magnitude less than h to 0

		double abso = abs(this->approx.coeffs(i));
		int coorder;
		if (abso == 1) {
			coorder = 0;
		}
		else {
			coorder = log10(abso);
		}

		int horder = log10(h);

		if (coorder <= horder) {
			this->approx.coeffs(i) = 0;
		}
	}

	this->approx.d = d;
	this->approx.order = ord;


	return;


}

double func::dxatpoint(double x) { //using the derivative approximation generated by the above function, finds the value of f'(x) 

	cout << "Evaluating dx^" << this->approx.d << "...";
	int size = this->approx.coeffs.size();
	double sum = 0;
	for (int i = 0; i < size; i++) {
		double point = this->approx.samps[i];

		double x = this->x + (point*h);


		double ans = this->define(x);

		ans *= this->approx.coeffs[i];


		sum += ans;

		
	}

	double check = pow(h, this->approx.d);
	int dif = log10(abs(sum)) - log10(check);

	if (dif >= abs(log10(h))) { //precision limit as above
		return 0;
	}

	sum /= pow(h, this->approx.d);

	if (log10(sum) < log10(check)) {
		sum = 0;
	}
	cout << "Done" << endl;
	return sum;
	
}



void func::taylor(int terms) { //runs the derivative approximation for n terms and generates a taylor polynomial
	double x = this->x;		   //around x using the formula T(x) = f(x) + x*dx + ((x^2)*dx^2)/2... + ((x^n)(dx^n))/n!
	double* taylorcos;
	taylorcos = new double[terms];
	taylorcos[0] = this->define(x);
	for (int i = 1; i < terms; i++) {
		this->finitedx(terms, i);
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

double func::numint(double a, double b) {
	int n = 500;
	double dist = b - a;
	double deltax = dist / n;
	double f0 = this->define(a);
	double fn = this->define(b);
	double j0area = 0;
	for (int i = 1; i <= n / 2; i++) {
		double j0 = (2 * i) - 1;
		double x = a + (deltax * j0);
		j0area += this->define(x);
	}
	j0area *= 4;
	double j1area = 0;
	for (int i = 1; i <= (n / 2) - 1; i++) {
		double j1 = (2 * i);
		double x = a + (deltax * j1);
		j1area += this->define(x);
	}
	j1area *= 2;
	double area = f0 + fn + j0area + j1area;
	area *= (deltax / 3);

	return area;
}