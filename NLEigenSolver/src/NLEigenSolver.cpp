#include "nlpch.h"

// Just to test
template <typename T>
void printMatrix(T& A)
{
	std::cout << A << std::endl;
};


int main(int argc, char* argv[])
{
	Log::Init();
	// Run application - NonLinearEigenSolver
	LOG_INFO("Nonlinear Eigenvalue solution!");

	//define 3x3 marix of doubles and set its entries to zero
	Eigen::MatrixXd matrixA(3, 3);
	matrixA.setZero();
	printMatrix(matrixA);

	std::cout << std::endl;
	Eigen::MatrixXd m = Eigen::MatrixXd::Random(5, 5);
	printMatrix(m);

	std::cout << std::endl;
	Eigen::Matrix3d m3d = Eigen::Matrix3d::Random();
	printMatrix(m3d);

	std::cout << "3*m3d = " << std::endl << std::endl;

	double c = 3.0;
	printMatrix(3 * m3d);

	m3d *= 3.0;
	std::cout << "3*m3d = " << std::endl << std::endl;
	printMatrix(m3d);

	std::cout << "inverse of m3d = " << std::endl << std::endl;
	printMatrix(m3d.inverse());

	std::cout << "Checking m3d*inv(m3d) = " << std::endl << std::endl;
	printMatrix(m3d * m3d.inverse());


	Eigen::MatrixXd  matrixB(3, 2);
	matrixB.setOnes();
	matrixB(0, 1) = 10;
	matrixB(1, 1) = 5;
	printMatrix(matrixB);


	std::cout << "A*B = " << std::endl;
	Eigen::MatrixXd D = matrixA * matrixB;
	printMatrix(D);
	return 0;
}