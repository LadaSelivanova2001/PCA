#include "Matrix.h"

#include<set>
#include<iomanip>
#include<stdexcept>
#include<cctype>
#include<cmath>

using namespace std;

void swap(vector<double>& a, vector<double>& b) {
	vector<double> tmp;
	tmp = a;
	a = b;
	b = tmp;
}

vector<double> operator+(const vector<double>& lhs, const vector<double>& rhs) {
	vector<double> result(lhs.size());
	for (int i = 0; i < lhs.size(); i++) {
		result[i] = lhs[i] + rhs[i];
	}
	return result;
}

vector<double> operator-(const vector<double>& lhs, const vector<double>& rhs) {
	vector<double> result(lhs.size());
	for (int i = 0; i < lhs.size(); i++) {
		result[i] = lhs[i] - rhs[i];
	}
	return result;
}

vector<double> operator*(const double& n, const vector<double>& rhs) {
	vector<double> result(rhs.size());
	for (int i = 0; i < rhs.size(); i++) {
		result[i] = n * rhs[i];
	}
	return result;
}

std::vector<double> operator*(const std::vector<double>& rhs, const double& n) {
	vector<double> result(rhs.size());
	for (int i = 0; i < rhs.size(); i++) {
		result[i] = n * rhs[i];
	}
	return result;
}

vector<double> operator/(const vector<double>& lhs, const double& n) {
	vector<double> result(lhs.size());
	for (int i = 0; i < lhs.size(); i++) {
		result[i] = lhs[i] / n;
	}
	return result;
}

bool is_zero(vector<double> v) { 
	bool flag = 1;
	for (const auto& elem : v) {
		if (elem != 0) {
			flag = 0;
			break;
		}
	}
	return flag;
}


//дополнение
vector<double>& Matrix::operator[](const int index) {
	return matrix[index];
}

vector<double> Matrix::operator[](const int index) const{
	return matrix[index];
}

Matrix Matrix::get_column(const int index) const {
	Matrix column(x_dimention, 1);
	for (int i = 0; i < x_dimention; i++) {
		column[i][0] = matrix[i][index];
	}
	return column;
}

Matrix Matrix::get_line(const int index) const {
	Matrix line(1, y_dimention);
	for (int i = 0; i < y_dimention; i++) {
		line[0][i] = matrix[index][i];
	}
	return line;
}

Matrix Matrix::add_column(Matrix matr) {
	vector<vector<double>> m = matrix;
	for (int i = 0; i < x_dimention; i++) {
		m[i].push_back(matr[i][0]);
	}
	Matrix result(m);
	return result;
}



void Matrix::read_from_file(ifstream& file) {
	double tmp;
	int tmp_y_dim = 0;
	x_dimention = 0;
	set<int> y_set;

	matrix.push_back({});
	while (true) {
		while (isspace(file.peek())) {
			file.ignore(1);
		}
		if (!isdigit(file.peek()) && file.peek() != '-') {
			throw invalid_argument("file content cannot be represented as a matrix");
		}
		matrix[x_dimention].push_back(0);
		file >> matrix[x_dimention][tmp_y_dim];
		while (isspace(file.peek()) && file.peek() != '\n') {
			file.ignore(1);
		}
		while (file.peek() != '\n' && !file.eof()) {
			tmp_y_dim++;
			matrix[x_dimention].push_back(0);
			if (!isdigit(file.peek()) && file.peek() != '-') {
				throw invalid_argument("file content cannot be represented as a matrix");
			}
			file >> matrix[x_dimention][tmp_y_dim];
			while (isspace(file.peek()) && file.peek() != '\n') {
				file.ignore(1);
			}
		}
		y_set.insert(tmp_y_dim);
		if (y_set.size() > 1) {
			throw invalid_argument("file content cannot be represented as a matrix");
		}
		tmp_y_dim = 0;
		if (file.eof()) { break; }
		else {
			matrix.push_back({});
			x_dimention++;
		}
	}
	x_dimention++;
	y_dimention = (*y_set.begin() + 1); // если множество пустое, знаит в цикл не вошли, значит матрица1х1
}

void Matrix::resize(const int& new_x, const int& new_y) {
	if (new_x < 0 || new_y < 0) {
		throw invalid_argument("matrix dimention connot be negative");
	}
	if (new_x == 0 || new_y == 0) {
		x_dimention = 0;
		y_dimention = 0;
		matrix = {};
	}
	else {
		x_dimention = new_x;
		y_dimention = new_y;
		matrix.resize(new_x);
		for (auto& line : matrix) {
			line.resize(new_y);
		}
	}
}

bool Matrix::is_vector() const {
	return (x_dimention == 1 || y_dimention == 1);
}

Matrix::Matrix() {
	x_dimention = 0;
	y_dimention = 0;
	matrix = {};
}

Matrix::Matrix(const int& x, const int& y) {
	x_dimention = x;
	y_dimention = y;
	matrix = {};

	if (x == 0 || y == 0) {
		throw invalid_argument("cannot create matrix with zero dimention");
	}
	matrix.resize(x);
	for (auto& line : matrix) {
		line.resize(y);
	}
}

Matrix::Matrix(const vector<vector<double>>& m) {
	if (m.empty() || m[0].empty()) {
		throw invalid_argument("cannot create empty matrix");
	}
	y_dimention = m[0].size();
	for (const auto& row : m) {
		if (row.size() != y_dimention) {
			throw invalid_argument("given content cannot be represent as a matrix");
		}
	}
	x_dimention = m.size();
	matrix = m;
}
/*
void Matrix::transpose() {
	if (y_dimention == 0 || x_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	int new_x_dim = y_dimention;
	int new_y_dim = x_dimention;
	vector<vector<double>> new_matrix(y_dimention);
	for (int i = 0; i < new_x_dim; i++) {
		new_matrix[i].resize(new_y_dim);
		for (int j = 0; j < new_y_dim; j++) {
			new_matrix[i][j] = matrix[j][i];
		}
	}
	y_dimention = new_y_dim;
	x_dimention = new_x_dim;
	matrix = new_matrix;
}
*/

Matrix Matrix::transpose() const {
	if (y_dimention == 0 || x_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	int new_x_dim = y_dimention;
	int new_y_dim = x_dimention;
	vector<vector<double>> new_matrix(y_dimention);
	for (int i = 0; i < new_x_dim; i++) {
		new_matrix[i].resize(new_y_dim);
		for (int j = 0; j < new_y_dim; j++) {
			new_matrix[i][j] = matrix[j][i];
		}
	}
	Matrix transposed_matr(new_matrix);
	return transposed_matr;
}

string Matrix::size() const {
	stringstream ss;
	ss << "(" << x_dimention << " x " << y_dimention << ")";
	return ss.str();
}

Matrix& Matrix::operator=(const Matrix& rhs) {
	(*this).resize(rhs.x_dimention, rhs.y_dimention);
	for (int i = 0; i < x_dimention; i++) {
		for (int j = 0; j < y_dimention; j++) {
			matrix[i][j] = rhs.matrix[i][j];
		}
	}
	return *this;
}

ostream& operator<<(ostream& out, const Matrix& m) {
	if (m.x_dimention == 0 || m.y_dimention == 0) {
		throw invalid_argument("connot print empty matrix");
	}
	long long n = 6;
	out << fixed << setprecision(n - 3);       
	for (const auto& line : m.matrix) {
		for (const auto& elem : line) {
			out << setw(n) << setfill(' ') << elem << "  ";
		}
		out << endl;
	}
	return out;
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
	Matrix result;
	if (lhs.x_dimention == 0 || lhs.y_dimention == 0) {
		throw invalid_argument("left matrix is empty");
	}
	if (rhs.x_dimention == 0 || rhs.y_dimention == 0) {
		throw invalid_argument("right matrix is empty");
	}
	if (lhs.x_dimention != rhs.x_dimention || lhs.y_dimention != rhs.y_dimention) {
		stringstream ss;
		ss << "cannot add a matrix of size " << lhs.size() << " to a matrix of size " << rhs.size() << endl;
		throw invalid_argument(ss.str());
	}
	result.resize(lhs.x_dimention, lhs.y_dimention);
	for (int i = 0; i < lhs.x_dimention; i++) {
		for (int j = 0; j < lhs.y_dimention; j++) {
			result.matrix[i][j] = lhs.matrix[i][j] + rhs.matrix[i][j];
		}
	}
	return result;
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs) {
	Matrix result;
	if (lhs.x_dimention == 0 || lhs.y_dimention == 0) {
		throw invalid_argument("left matrix is empty");
	}
	if (rhs.x_dimention == 0 || rhs.y_dimention == 0) {
		throw invalid_argument("right matrix is empty");
	}
	if (lhs.x_dimention != rhs.x_dimention || lhs.y_dimention != rhs.y_dimention) {
		stringstream ss;
		ss << "cannot be subtracted from a matrix of size " << lhs.size() << " matrix size " << rhs.size() << endl;
		throw invalid_argument(ss.str());
	}
	result.resize(lhs.x_dimention, lhs.y_dimention);
	for (int i = 0; i < lhs.x_dimention; i++) {
		for (int j = 0; j < lhs.y_dimention; j++) {
			result.matrix[i][j] = lhs.matrix[i][j] - rhs.matrix[i][j];
		}
	}
	return result;
}

Matrix operator*(const double& n, const Matrix& rhs) {
	Matrix result;
	if (rhs.x_dimention == 0 || rhs.y_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	result.resize(rhs.x_dimention, rhs.y_dimention);
	for (int i = 0; i < rhs.x_dimention; i++) {
		for (int j = 0; j < rhs.y_dimention; j++) {
			result.matrix[i][j] = n * rhs.matrix[i][j];
		}
	}
	for (auto& row : result.matrix) {
		for (auto& elem : row) {
			if (abs(elem) < 1E-9) {
				elem = 0;
			}
		}
	}
	return result;
}

Matrix operator*(const Matrix& lhs, const double& n) {
	Matrix result;
	if (lhs.x_dimention == 0 || lhs.y_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	result.resize(lhs.x_dimention, lhs.y_dimention);
	for (int i = 0; i < lhs.x_dimention; i++) {
		for (int j = 0; j < lhs.y_dimention; j++) {
			result.matrix[i][j] = n * lhs.matrix[i][j];
		}
	}
	for (auto& row : result.matrix) {
		for (auto& elem : row) {
			if (abs(elem) < 1E-9) {
				elem = 0;
			}
		}
	}
	return result;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
	Matrix result;
	if (lhs.x_dimention == 0 || lhs.y_dimention == 0) {
		throw invalid_argument("left matrix is empty");
	}
	if (rhs.x_dimention == 0 || rhs.y_dimention == 0) {
		throw invalid_argument("right matrix is empty");
	}
	if (lhs.y_dimention != rhs.x_dimention) {
		stringstream ss;
		ss << "matrices of size " << lhs.size() << " and size " << rhs.size() << " cannot be multiplied" << endl;
		throw invalid_argument(ss.str());
	}
	result.resize(lhs.x_dimention, rhs.y_dimention);
	for (int i = 0; i < lhs.x_dimention; i++) {
		for (int j = 0; j < rhs.y_dimention; j++) {
			result.matrix[i][j] = 0;
			for (int k = 0; k < lhs.y_dimention; k++) {
				result.matrix[i][j] += lhs.matrix[i][k] * rhs.matrix[k][j];
			}
		}
	}
	for (auto& row : result.matrix) {
		for (auto& elem : row) {
			if (abs(elem) < 1E-8) {
				elem = 0;
			}
		}
	}
	return result;
}

//скал€рное произведение
double operator&(const Matrix& lhs, const Matrix& rhs) {
	if (lhs.x_dimention == 0 || lhs.y_dimention == 0) {
		throw invalid_argument("left matrix is empty");
	}
	if (rhs.x_dimention == 0 || rhs.y_dimention == 0) {
		throw invalid_argument("right matrix is empty");
	}
	if (!lhs.is_vector() && !rhs.is_vector()) {
		stringstream ss;
		ss << "scalar product for " << lhs.size() << " matrix and " << rhs.size()
			<< " matrix " << " is undefined" << endl;
		throw invalid_argument(ss.str());
	}
	double res = 0;
	if (lhs.x_dimention == 1 && rhs.x_dimention == 1) {
		for (int i = 0; i < lhs.y_dimention; i++) {
			res += lhs.matrix[0][i] * rhs.matrix[0][i];
		}
	}
	else {
		Matrix lhs_copy, rhs_copy;
		lhs_copy = lhs;
		rhs_copy = rhs;
		if (lhs.x_dimention > 1) {
			lhs_copy = lhs_copy.transpose();
		}
		if (rhs.x_dimention > 1) {
			rhs_copy = rhs_copy.transpose();
		}
		for (int i = 0; i < lhs_copy.y_dimention; i++) {
			res += lhs_copy.matrix[0][i] * rhs_copy.matrix[0][i];
		}
	}
	if (res < 1E-9) {
		res = 0;
	}
	return res;
}

//произведение јдамара
Matrix operator%(const Matrix& lhs, const Matrix& rhs) {
	Matrix result;
	if (lhs.x_dimention == 0 || lhs.y_dimention == 0) {
		throw invalid_argument("left matrix is empty");
	}
	if (rhs.x_dimention == 0 || rhs.y_dimention == 0) {
		throw invalid_argument("right matrix is empty");
	}
	if (lhs.x_dimention != rhs.x_dimention || lhs.y_dimention != rhs.y_dimention) {
		stringstream ss;
		ss << "Hadamard multiplication for matrices of size " << lhs.size()
			<< " and " << rhs.size() << " is impossible" << endl;
		throw invalid_argument(ss.str());
	}
	result.resize(lhs.x_dimention, lhs.y_dimention);
	for (int i = 0; i < lhs.x_dimention; i++) {
		for (int j = 0; j < lhs.y_dimention; j++) {
			result.matrix[i][j] = lhs.matrix[i][j] * rhs.matrix[i][j];
		}
	}
	return result;
}




/*****************************************
* ќпределени€ методов дл€ Symmetric_matrix
*****************************************/

ifstream& operator>>(ifstream& in, Symmetric_matrix& S_m) {
	if (!in.is_open()) {
		throw runtime_error("file for reading is closed");
	}
	if (in.eof()) {
		throw runtime_error("file is empty");
	}
	S_m.read_from_file(in);
	if (S_m.x_dimention != S_m.y_dimention) {
		throw invalid_argument("symmetric matrix must be square");
	}
	else {
		for (int i = 0; i < S_m.x_dimention; i++) {
			for (int j = 0; j < S_m.y_dimention; j++) {
				if (S_m.matrix[i][j] != S_m.matrix[j][i]) {
					throw invalid_argument("file content is not a symmetric matrix");
				}
			}
		}
	}
	return in;
}


/****************************************
* ќпределени€ методов дл€ Diagonal_matrix
****************************************/

Diagonal_matrix::Diagonal_matrix(const int& size, const vector<double>& content) {
	if (size == 0) {
		throw invalid_argument("cannot create an empty matrix");
	}
	if (content.size() != size) {
		throw invalid_argument("the dimensions of the matrix and the given vector do not match");
	}
	x_dimention = size;
	y_dimention = size;
	matrix.resize(size);
	for (int i = 0; i < size; i++) {
		matrix[i].assign(size, 0.0);
		matrix[i][i] = content[i];
	}
}

ifstream& operator>>(ifstream& in, Diagonal_matrix& D_m) {
	if (!in.is_open()) {
		throw runtime_error("file for reading is closed");
	}
	if (in.eof()) {
		throw runtime_error("file is empty");
	}
	D_m.read_from_file(in);
	if (D_m.x_dimention != D_m.y_dimention) {
		throw invalid_argument("diagonal matrix must be square");
	}
	else {
		for (int i = 0; i < D_m.x_dimention; i++) {
			for (int j = 0; j < D_m.y_dimention; j++) {
				if (i != j && D_m.matrix[i][j] != 0) {
					throw invalid_argument("file content is not a diagonal matrix");
				}
			}
		}
	}
	return in;
}


/****************************************
* ќпределени€ методов дл€ Identity_matrix
****************************************/

Identity_matrix::Identity_matrix(const int& size) {
	x_dimention = size;
	y_dimention = size;
	matrix.resize(size);
	for (int i = 0; i < size; i++) {
		matrix[i].assign(size, 0.0);
		matrix[i][i] = 1.0;
	}
}

ifstream& operator>>(ifstream& in, Identity_matrix& I_m) {
	if (!in.is_open()) {
		throw runtime_error("file for reading is closed");
	}
	if (in.eof()) {
		throw runtime_error("file is empty");
	}
	I_m.read_from_file(in);
	if (I_m.x_dimention != I_m.y_dimention) {
		throw invalid_argument("identity matrix must be square");
	}
	else {
		for (int i = 0; i < I_m.x_dimention; i++) {
			for (int j = 0; j < I_m.y_dimention; j++) {
				if (i == j && I_m.matrix[i][j] != 1) {
					throw invalid_argument("file content is not a identity matrix");
				}
				if (i != j && I_m.matrix[i][j] != 0) {
					throw invalid_argument("file content is not a identity matrix");
				}
			}
		}
	}
	return in;
}


/************************************************
* ќпределени€ методов дл€ Upper_triangular_matrix
************************************************/

Upper_triangular_matrix::Upper_triangular_matrix(const vector<vector<double>>& U_t_m) {
	if (U_t_m.empty() || U_t_m[0].empty()) {
		throw invalid_argument("cannot create empty matrix");
	}
	int y_dimention = U_t_m[0].size();
	for (const auto& row : U_t_m) {
		if (row.size() != y_dimention) {
			throw invalid_argument("given content cannot be represent as a matrix");
		}
	}
	x_dimention = U_t_m.size();
	matrix = U_t_m;
	for (int i = 0; i < y_dimention; i++) {
		for (int j = i + 1; j < x_dimention; j++) {
			if (matrix[j][i] != 0) {
				throw invalid_argument("file content cannot be represented as a Upper triangular matrix");
			}
		}
	}
}

ifstream& operator>>(ifstream& in, Upper_triangular_matrix& U_t_m) {
	if (!in.is_open()) {
		throw runtime_error("file for reading is closed");
	}
	if (in.eof()) {
		throw runtime_error("file is empty");
	}
	U_t_m.read_from_file(in);
	for (int i = 0; i < U_t_m.y_dimention; i++) {
		for (int j = i + 1; j < U_t_m.x_dimention; j++) {
			if (U_t_m.matrix[j][i] != 0) {
				throw invalid_argument("file content cannot be represented as a Upper triangular matrix");
			}
		}
	}
	return in;
}

Lower_triangular_matrix Upper_triangular_matrix::transpose() const {
	if (y_dimention == 0 || x_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	int new_x_dim = y_dimention;
	int new_y_dim = x_dimention;
	vector<vector<double>> new_matrix(y_dimention);
	for (int i = 0; i < new_x_dim; i++) {
		new_matrix[i].resize(new_y_dim);
		for (int j = 0; j < new_y_dim; j++) {
			new_matrix[i][j] = matrix[j][i];
		}
	}
	Lower_triangular_matrix result(new_matrix);
}

/************************************************
* ќпределени€ методов дл€ Lower_triangular_matrix
************************************************/

Lower_triangular_matrix::Lower_triangular_matrix(const vector<vector<double>>& L_t_m) {
	if (L_t_m.empty() || L_t_m[0].empty()) {
		throw invalid_argument("cannot create empty matrix");
	}
	int y_dimention = L_t_m[0].size();
	for (const auto& row : L_t_m) {
		if (row.size() != y_dimention) {
			throw invalid_argument("given content cannot be represent as a matrix");
		}
	}
	x_dimention = L_t_m.size();
	matrix = L_t_m;
	for (int i = 0; i < x_dimention; i++) {
		for (int j = i + 1; j < y_dimention; j++) {
			if (matrix[i][j] != 0) {
				throw invalid_argument("file content cannot be represented as a Lower triangular matrix");
			}
		}
	}
}

ifstream& operator>>(ifstream& in, Lower_triangular_matrix& L_t_m) {
	if (!in.is_open()) {
		throw runtime_error("file for reading is closed");
	}
	if (in.eof()) {
		throw runtime_error("file is empty");
	}
	L_t_m.read_from_file(in);
	for (int i = 0; i < L_t_m.x_dimention; i++) {
		for (int j = i + 1; j < L_t_m.y_dimention; j++) {
			if (L_t_m.matrix[i][j] != 0) {
				throw invalid_argument("file content cannot be represented as a Lower triangular matrix");
			}
		}
	}
	return in;
}

Upper_triangular_matrix Lower_triangular_matrix::transpose() const {
	if (y_dimention == 0 || x_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	int new_x_dim = y_dimention;
	int new_y_dim = x_dimention;
	vector<vector<double>> new_matrix(y_dimention);
	for (int i = 0; i < new_x_dim; i++) {
		new_matrix[i].resize(new_y_dim);
		for (int j = 0; j < new_y_dim; j++) {
			new_matrix[i][j] = matrix[j][i];
		}
	}
	Upper_triangular_matrix result(new_matrix);
}




//лаба 2
double Matrix::trace() const {
	double trace = 0;
	if (x_dimention == 0 || y_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	if (x_dimention != y_dimention) {
		throw invalid_argument("trace operation not defined for rectangular matrix");
	}
	for (int i = 0; i < x_dimention; i++) {
		trace += matrix[i][i];
	}
	return trace;
}

Matrix Matrix::gaussian_algorithm() const {
	int iMax = 0, countSwaps = 1;
	double q = 1;
	Matrix result;
	result = (*this);
	for (int i = 0; i < x_dimention; ++i) {
		iMax = i;
		for (int j = i + 1; j < x_dimention; j++) {
			if (abs(result.matrix[j][i]) > abs(result.matrix[iMax][i])) {
				iMax = j;
			}
		}
		if (abs(result.matrix[iMax][i]) < 1E-9) {
			continue;
		}
		swap(result.matrix[i], result.matrix[iMax]);
		countSwaps = -countSwaps * (i != iMax ? 1 : -1);
		for (int j = i + 1; j < x_dimention; ++j) {
			if (result.matrix[i][i] != 0) {
				q = -result.matrix[j][i] / result.matrix[i][i];
			}
			else {
				q = 1;
			}
			result.matrix[j] = result.matrix[j] + q * result.matrix[i];
		}
	}
	result = countSwaps * result;
	for (auto& row : result.matrix) {
		for (auto& elem : row) {
			if (abs(elem) < 1E-9) {
				elem = 0;
			}
		}
	}
	return result;
}

double Matrix::determinant() const {
	double result = 1;
	if (x_dimention == 0 || y_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	if (x_dimention != y_dimention) {
		throw invalid_argument("determinant operation not defined for rectangular matrix");
	}
	Matrix C;
	C = (*this).gaussian_algorithm();
	for (int i = 0; i < x_dimention; i++) {
		result *= C.matrix[i][i];
	}
	if (result == 0) {
		return 0;
	}
	return result;
}


//это можно использовать и дл€ матриц и дл€ "векторов", это вызовет воответствующий метод
double Matrix::norm() const {
	if (x_dimention == 0 || y_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	double result = 0;
	if ((*this).is_vector()) {
		result = (*this).vector_norm();
	}
	else {
		result = (*this).matrix_norm();
	}
	return result;
}

double Matrix::vector_norm() const {
	double res = 0, sum = 0;
	if (x_dimention > 1) {
		Matrix tmp;
		tmp = (*this);
		tmp = tmp.transpose();
		for (const auto& elem : tmp.matrix[0]) {
			sum += elem * elem;
		}
		res = sqrt(sum);
	}
	else {
		for (int i = 0; i < x_dimention; i++) {
			sum += matrix[i][0] * matrix[i][0];
		}
		res = sqrt(sum);
	}
	return res;
}

double Matrix::matrix_norm() const {
	double res = 0, sum = 0;
	for (const auto& row : matrix) {
		for (const auto& elem : row) {
			sum += elem * elem;
		}
	}
	res = sqrt(sum);
	return res;
}

double Matrix::max_norm() const {
	if (x_dimention == 0 || y_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	if (!(*this).is_vector()) {
		throw invalid_argument("'maximum norm' is undefined for the matrix");
	}
	double result = 0;
	for (const auto& row : matrix) {
		for (const auto& elem : row) {
			if (abs(elem) > result) {
				result = abs(elem);
			}
		}
	}
	return result;
}




//лаба 3
int Matrix::rank() const {
	int result = 1, rank = 0, flag = 0;
	if (x_dimention == 0 || y_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	Matrix C;
	if (x_dimention > y_dimention) {
		Matrix tmp;
		tmp = (*this);
		tmp = tmp.transpose();
		C = tmp.gaussian_algorithm();
		for (const auto& row : C.matrix) {
			if (!is_zero(row)) {
				rank++;
			}
		}
	}
	else {
		C = (*this).gaussian_algorithm();
		for (const auto& row : C.matrix) {
			if (!is_zero(row)) {
				rank++;
			}
		}
	}
	return rank;
}

//угол между векторами (cos)
double operator^(const Matrix& lhs, const Matrix& rhs) {
	if (!lhs.is_vector() || !rhs.is_vector()) {
		throw invalid_argument("operation '^' is undefined for matrix");
	}
	if (lhs.x_dimention == 0 || lhs.y_dimention == 0) {
		throw invalid_argument("left vector is empty");
	}
	if (rhs.x_dimention == 0 || rhs.y_dimention == 0) {
		throw invalid_argument("right vector is empty");
	}
	double result = 0, norm1 = 0, norm2 = 0;
	Matrix tmp;
	norm1 = lhs.norm();
	norm2 = rhs.norm();
	if (norm1 == 0 || norm2 == 0) {
		throw invalid_argument("operation '^' is undefined for zero vectors");
	}
	tmp = lhs * rhs;
	result = tmp.matrix[0][0] / (norm1 * norm2);
	return result;
}

Matrix Matrix::inverse() const {
	Matrix result;
	if (x_dimention == 0 || y_dimention == 0) {
		throw invalid_argument("matrix is empty");
	}
	if (x_dimention != y_dimention) {
		throw invalid_argument("inverse matrix does not exist for rectangular matrices");
	}
	if ((*this).determinant() == 0) {
		throw invalid_argument("inverse matrix does not exist for degenerate matrices");
	}
	else {
		Identity_matrix I(x_dimention);
		result = I;
		int iMax = 0;
		double q = 1;
		Matrix tmp;
		tmp = (*this);
		for (int i = 0; i < x_dimention; i++) {
			if (tmp.matrix[i][i] == 0) {
				for (int k = i + 1; k < x_dimention; k++) {
					if (tmp.matrix[k][i] != 0) {
						swap(result.matrix[i], result.matrix[k]);
						swap(tmp.matrix[i], tmp.matrix[k]);
						break;
					}
				}
			}
			result.matrix[i] = result.matrix[i] / tmp.matrix[i][i];
			tmp.matrix[i] = tmp.matrix[i] / tmp.matrix[i][i];
			for (int k = 0; k < x_dimention; k++) {
				if (k == i) {
					continue;
				}
				result.matrix[k] = result.matrix[k] - tmp.matrix[k][i] * result.matrix[i];
				tmp.matrix[k] = tmp.matrix[k] - tmp.matrix[k][i] * tmp.matrix[i];
			}
		}
		for (auto& row : result.matrix) {
			for (auto& elem : row) {
				if (abs(elem) < 1E-9) {
					elem = 0;
				}
			}
		}
	}
	return result;
}




//лаба 4

ifstream& operator>>(ifstream& in, Matrix& m) {
	if (!in.is_open()) {
		throw runtime_error("file for reading is closed");
	}
	if (in.eof()) {
		throw runtime_error("file is empty");
	}
	m.read_from_file(in);
	return in;
}

ofstream& operator<<(ofstream& out, const Matrix& m) {
	if (!out.is_open()) {
		throw runtime_error("file for writing is closed");
	}
	if (m.x_dimention == 0 || m.y_dimention == 0) {
		throw invalid_argument("connot print empty matrix");
	}
	int n = 6;
	out << fixed << setprecision(n - 3);        //ѕочему переполнение????????????77??7777
	for (const auto& line : m.matrix) {
		for (const auto& elem : line) {
			out << setw(n) << setfill(' ') << elem << "  ";
		}
		out << endl;
	}
	return out;
}

/****************************************
* ¬ начало файла записываютс€ размерности
* x_dimention и y_dimention, только потом
* записываетс€ сама матрица
****************************************/
void Matrix::write_to_binary(const std::string& path) {
	ofstream file(path, ios::out | ios::binary);
	if (!file.is_open()) {
		throw runtime_error("file for writing is closed");
	}
	file.write((const char*)&x_dimention, sizeof(int));
	file.write((const char*)&y_dimention, sizeof(int));
	for (const auto& row : matrix) {
		for (const auto& elem : row) {
			file.write((const char*)&elem, sizeof(double));
		}
	}
	file.close();
}

void Matrix::read_from_binary(const std::string& path) {
	ifstream file(path, ios::in | ios::binary);
	if (!file.is_open()) {
		throw runtime_error("file for reading is closed");
	}
	file.read((char*)&x_dimention, sizeof(int));
	file.read((char*)&y_dimention, sizeof(int));
	(*this).resize(x_dimention, y_dimention);
	for (auto& row : matrix) {
		for (auto& elem : row) {
			file.read((char*)&elem, sizeof(double));
		}
	}
	file.close();
}