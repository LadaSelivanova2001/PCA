#include"Matrix.h"
#include"PCA.h"
#include<vector>
#include<algorithm>

using namespace std;

int main()
{
	setlocale(LC_ALL, "Russian");
	ifstream f1("data.txt", ios::in); //��������� ������ �� �������� ������
	ifstream f2("scores.txt", ios::in);
	ifstream f3("loadings.txt", ios::in);

	Matrix data;
	Matrix check_scores;
	Matrix check_loadings;

	f1 >> data;
	f2 >> check_scores; //check_scores � check_loadings ���������� ��� �������� ���������� scores � loadings 
	f3 >> check_loadings;

	int PC = min(data[0].size(), data.transpose()[0].size());


	PCA pca(data);
	pca.centering(); //�������������
	pca.scaling(); //������������
	pca.NIPALS(PC); //��������� ������ ������, ����� � ��������

	cout << "��������� ������" << endl;
	cout << data << endl;
	cout << "������ ����� ������������� � ������������" << endl;
	cout << pca.get_data() << endl;
	cout << "���������� ������� ������" << endl;
	cout << pca.get_scores() << endl;
	cout << "������� ������ ��� ���������" << endl;
	cout << check_scores << endl;
	cout << "���������� ������� �����" << endl;
	cout << pca.get_loadings() << endl;
	cout << "������� ����� ��� ���������" << endl;
	cout << check_loadings << endl;
	cout << "�������" << endl;
	cout << pca.leverage() << endl;
	cout << "����������" << endl;
	cout << pca.deviation(PC) << endl;
	cout << "������ ���������" << endl;
	cout << pca.TRV(PC) << endl;
	cout << "����������� ���������" << endl;
	cout << pca.ERV(PC) << endl;

	f1.close();
	f2.close();
	f3.close();

	return 0;
}