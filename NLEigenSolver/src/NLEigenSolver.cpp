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
	LOG_ERROR("Check error message!");
	LOG_WARN("Check warning message!");
	//define 3x3 marix of doubles and set its entries to zero
	Eigen::MatrixXd matrixA(3, 3);
	matrixA.setZero();
	printMatrix(matrixA);

	std::cout << std::endl;
	Eigen::MatrixXd m = Eigen::MatrixXd::Random(5, 5);

	LOG_ERROR("Matrix m: \n {0}", m);
	LOG_WARN("Block of m1: \n {0}", m.block<3,2>(0,0));
	LOG_WARN("Block of m2: \n {0}", m.block<3, 2>(1, 1));
	LOG_INFO("Sum of m1 + 2*m2: \n {0}", m.block<3, 2>(0, 0) + 2.0*m.block<3, 2>(1, 1));

	LOG_INFO("Get the col3 of m: \n {0}", m.col(2));
	LOG_INFO("Get the m*m.col3 =  \n {0}", m*m.col(2));

	LOG_INFO("Transpose of the col3 of m: \n {0}", m.col(2).transpose());
	LOG_INFO("col3*col3.transpose =  \n {0}", m.col(2)*m.col(2).transpose());
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