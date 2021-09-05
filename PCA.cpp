#include"PCA.h"
#include"Matrix.h"
#include<algorithm>

using namespace std;


PCA::PCA() 
{
	Matrix data;
	Matrix scores;// матрица счетов
	Matrix loadings;// матрица весов
	Matrix tailings;// матрица остатков
}

PCA::PCA(Matrix a)
{
	data = a;
}

Matrix PCA::centering() //центрирование (вычитание из каждого столбца среднего по столбцу значения)
{ 
	int width = data[0].size();
	int height = data.transpose()[0].size();
	Matrix result;
	result = data.transpose(); //переносим данные в result
	double average_column = 0;
	double sum_column = 0;
	for (int i = 0; i < data[0].size(); i++) {
		for (int j = 0; j < result[i].size(); j++) {
			sum_column += result[i][j]; //сумма элементов каждого столбца
		}
		average_column = sum_column / result[i].size(); //среднее каждого столбца
		for (int j = 0; j < result[i].size(); j++) {
			result[i][j] -= average_column;
		}
		sum_column = 0;
		average_column = 0;
	}
	result = result.transpose();
	data = result;


	return result;
}

Matrix PCA::scaling() //шкалирование (каждый элемент столбца делится на стандартное отклонение столбца)
{ 
	try {
		int width = data[0].size();
		int height = data.transpose()[0].size();
		Matrix result;
		result = data.transpose();
		double sqr_sum_column = 0; 
		double deviation = 0;
		for (int i = 0; i < data[0].size(); i++) {
			for (int j = 0; j < result[i].size(); j++) {
				sqr_sum_column += pow(result[i][j], 2); //сумма квадратов элементов столбца
			}
			deviation = sqrt(sqr_sum_column / result[i].size()); //находим стандартное отклонение
			result[i] = result[i] / deviation;
			sqr_sum_column = 0;
			deviation = 0;
		}
		result = result.transpose();
        data = result;
		return result;
	}
	catch (...) {
		cerr << "Ошибка" << endl;
	}
}

void PCA::NIPALS(const int PC) 
{
	try {
		Matrix E, P, T, t, t0, p, d;
		E = data;
		for (int h = 0; h < PC; h++) {
			t = E.get_column(h);

			do {
				p = (t.transpose() * E) * (1 / (t & t));
				p = p.transpose();
				p = p * (1 / p.norm());
				t0 = t;
				t = (E * p) * (1 / (p & p));
				d = t0 - t;
			} while (d.norm() >= 1e-8); //проверка сходимости

			E = E - t * p.transpose();
			if (h == 0) {
				P = p;
				T = t;
			}
			else {
				P = P.add_column(p); //формирование матриц счетов и весов
				T = T.add_column(t);
			}
		}

		scores = T;
		loadings = P;
		tailings = E;

		return;
	}
	catch (...) {
		cerr << "Ошибка" << endl;
	}
}

Matrix PCA::leverage() const //вычисление размахов
{ 
	try {
		vector<vector<double>> leverage;
		leverage.push_back({});
		Matrix A, t;
		int height = scores.transpose()[0].size();
		double h = 0;

		A = (scores.transpose() * scores).inverse();
		for (int i = 0; i < height; i++) {
			t = scores.get_line(i);
			h = ((t * A) & t);   
			leverage[0].push_back(h);
		}

		Matrix result(leverage);
		result = result.transpose();
		return result;
	}
	catch (...) {
		cerr << "Ошибка" << endl;
	}
}

Matrix PCA::deviation(const int PC) const //вычисление отклонений
{ 
	try {
		Matrix residuals, t, p;
		residuals = data;
		for (int i = 0; i < PC; i++) {
			t = scores.get_column(i);
			p = loadings.get_column(i);
			residuals = residuals - (t * p.transpose()); 
		}
		
		vector<vector<double>> deviation;
		deviation.push_back({});
		Matrix r;
		int width = residuals[0].size();
		int height = residuals.transpose()[0].size();
	    
		double v = 0;

		for (int i = 0; i < height; i++) {
			r = residuals.get_line(i);
			for (int j = 0; j < width; j++) {
				v += pow(r[0][j], 2);
			}
			deviation[0].push_back(v);
			v = 0;
		}

		Matrix result(deviation);
		result = result.transpose();
		return result;
	}
	catch (...) {
		cerr << "Ошибка" << endl;
	}
}

double PCA::TRV(const int PC) const //полная дисперсия
{ 
	try {
		Matrix v;
		v = (*this).deviation(PC);
		double trv = 0, v_0 = 0;

		int I = v.transpose()[0].size(),
			J = data[0].size();

		for (int i = 0; i < I; i++) {
			v_0 += v[i][0];
		}
		v_0 = v_0 / I;

		trv = v_0 / J;

		return trv;
	}
	catch (...) {
		cerr << "Ошибка" << endl;
	}
}

double PCA::ERV(const int PC) const //объясненная дисперсия
{ 
	try {
		double erv = 0;
		erv = 1 - (*this).TRV(PC) / (*this).TRV(0);
		return erv;
	}
	catch (...) {
		cerr << "Ошибка" << endl;
	}
}

Matrix PCA::get_data() const
{
	return data;
}

Matrix PCA::get_scores() const
{
	return scores;
}

Matrix PCA::get_loadings() const
{
	return loadings;
}

Matrix PCA::get_tailings() const
{
	return tailings;
}
