#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>

class Matrix 
{
protected:
	int x_dimention;
	int y_dimention;
	std::vector<std::vector<double>> matrix;

	void read_from_file(std::ifstream& file);
	void resize(const int& new_x, const int& new_y);
	bool is_vector() const;

	//лаба 2
	Matrix gaussian_algorithm() const;
	double vector_norm() const;
	double matrix_norm() const;

public:
	Matrix();
	Matrix(const int& x, const int& y);
	Matrix(const std::vector<std::vector<double>>& m);

	std::string size() const;
	Matrix& operator=(const Matrix& rhs);

	friend std::ostream& operator<<(std::ostream& out, const Matrix& m);
	friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
	friend Matrix operator-(const Matrix& lhs, const Matrix& rhs);
	friend Matrix operator*(const double& n, const Matrix& rhs);
	friend Matrix operator*(const Matrix& lhs, const double& n);
	friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);
	friend Matrix operator%(const Matrix& lhs, const Matrix& rhs);      //произведение јдамара
	friend double operator&(const Matrix& lhs, const Matrix& rhs);      //скал€рное произведение

	//лаба 2
	double trace() const;
	double determinant() const;
	double norm() const;
	double max_norm() const;

	//лаба 3
	int rank() const;
	friend double operator^(const Matrix& lhs, const Matrix& rhs);      //угол между векторами (cos)
	Matrix inverse() const;
	//void transpose();
	Matrix transpose() const;

	//лаба 4
	friend std::ifstream& operator>>(std::ifstream& in, Matrix& m);
	friend std::ofstream& operator<<(std::ofstream& out, const Matrix& m);
	void write_to_binary(const std::string& path);
	void read_from_binary(const std::string& path);

	//дополнение
	std::vector<double>& operator[](const int index);
	std::vector<double> operator[](const int index) const;
	Matrix get_column(const int index) const;
	Matrix get_line(const int index) const;
	Matrix add_column(Matrix matr);
};




class Symmetric_matrix : public Matrix {
public:
	Symmetric_matrix() {};
	friend std::ifstream& operator>>(std::ifstream& in, Symmetric_matrix& m);
};

class Diagonal_matrix : public Symmetric_matrix {
public:
	Diagonal_matrix() {};
	Diagonal_matrix(const int& size, const std::vector<double>& content);
	friend std::ifstream& operator>>(std::ifstream& in, Diagonal_matrix& m);
};

class Identity_matrix : public Diagonal_matrix {
public:
	Identity_matrix(const int& size);
	friend std::ifstream& operator>>(std::ifstream& in, Identity_matrix& m);
}; 



class Lower_triangular_matrix;
class Upper_triangular_matrix : public Matrix {
public:
	Upper_triangular_matrix() {};
	Upper_triangular_matrix(const std::vector<std::vector<double>>& U_t_m);
	friend std::ifstream& operator>>(std::ifstream& in, Upper_triangular_matrix& U_t_m);
	Lower_triangular_matrix transpose() const;
};

class Lower_triangular_matrix : public Matrix {
public:
	Lower_triangular_matrix() {};
	Lower_triangular_matrix(const std::vector<std::vector<double>>& L_t_m);
	friend std::ifstream& operator>>(std::ifstream& in, Lower_triangular_matrix& L_t_m);
	Upper_triangular_matrix transpose() const;
};

// из cpp дл€ арифметики с векторами
std::vector<double> operator+(const std::vector<double>& lhs, const std::vector<double>& rhs);

std::vector<double> operator-(const std::vector<double>& lhs, const std::vector<double>& rhs);

std::vector<double> operator*(const double& n, const std::vector<double>& rhs);

std::vector<double> operator*(const double& n, const std::vector<double>& rhs);

std::vector<double> operator*(const std::vector<double>& rhs, const double& n);

std::vector<double> operator/(const std::vector<double>& lhs, const double& n);
