#include"Matrix.h"
#include"PCA.h"
#include<vector>
#include<algorithm>

using namespace std;

int main()
{
	setlocale(LC_ALL, "Russian");
	ifstream f1("data.txt", ios::in); //считываем данные из исходных файлов
	ifstream f2("scores.txt", ios::in);
	ifstream f3("loadings.txt", ios::in);

	Matrix data;
	Matrix check_scores;
	Matrix check_loadings;

	f1 >> data;
	f2 >> check_scores; //check_scores и check_loadings необходимы для проверки полученных scores и loadings 
	f3 >> check_loadings;

	int PC = min(data[0].size(), data.transpose()[0].size());


	PCA pca(data);
	pca.centering(); //центрирование
	pca.scaling(); //шкалирование
	pca.NIPALS(PC); //получение матриц счетов, весов и остатков

	cout << "Начальные данные" << endl;
	cout << data << endl;
	cout << "Данные после центрирования и шкалирования" << endl;
	cout << pca.get_data() << endl;
	cout << "Полученная матрица счётов" << endl;
	cout << pca.get_scores() << endl;
	cout << "Матрица счётов для сравнения" << endl;
	cout << check_scores << endl;
	cout << "Полученная матрица весов" << endl;
	cout << pca.get_loadings() << endl;
	cout << "Матрица весов для сравнения" << endl;
	cout << check_loadings << endl;
	cout << "Размахи" << endl;
	cout << pca.leverage() << endl;
	cout << "Отклонения" << endl;
	cout << pca.deviation(PC) << endl;
	cout << "Полная дисперсия" << endl;
	cout << pca.TRV(PC) << endl;
	cout << "Объясненная дисперсия" << endl;
	cout << pca.ERV(PC) << endl;

	f1.close();
	f2.close();
	f3.close();

	return 0;
}