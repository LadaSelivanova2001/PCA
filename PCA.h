#pragma once

#include"Matrix.h"

class PCA 
{
protected:
	Matrix data;
	Matrix scores;
	Matrix loadings;
	Matrix tailings;
public:
	PCA(Matrix matr);
	PCA();
	Matrix get_data() const;
	Matrix get_scores() const;
	Matrix get_loadings() const;
	Matrix get_tailings() const;
	Matrix centering();
	Matrix scaling();
	void NIPALS(const int PC);
	Matrix leverage() const;
	Matrix deviation(const int PC) const;
	double TRV(const int PC) const;
	double ERV(const int PC) const;
};