
#pragma once

namespace Geo
{

struct Solid
{
	double density;
	double thermal_expansion;
	double poisson_ratio;
	double Youngs_modulus;
};

inline void calculateLameConstant(const double nv, const double E, double &Lambda, double &G, double &K)
{
	Lambda = E * nv / ((1. + nv) * (1. - 2. * nv));
	G = 0.5 * E / (1. + nv);
	K = (3.0 * Lambda + 2.0 * G) / 3.0;
}

template<class T_MATRIX>
inline void setElasticConsitutiveTensor(const size_t Dimension, double Lambda, double G, T_MATRIX &D_e)
{
	//D_e *= 0.0;
	D_e(0,0) = Lambda + 2 * G;
	D_e(0,1) = Lambda;
	D_e(0,2) = Lambda;

	D_e(1,0) = Lambda;
	D_e(1,1) = Lambda + 2 * G;
	D_e(1,2) = Lambda;

	D_e(2,0) = Lambda;
	D_e(2,1) = Lambda;
	D_e(2,2) = Lambda + 2 * G;

	D_e(3,3) = G;
	//Plane stress
	// plane stress, only for test
	//(*D_e)(0,0) = (1.0-Mu)*Lambda + 2 * G;
	//(*D_e)(0,1) = Lambda;

	//(*D_e)(1,0) = Lambda;
	// (*D_e)(1,1) = (1.0-Mu)*Lambda + 2 * G;
	// (*D_e)(3,3) = G;

	if(Dimension == 3)
	{
		D_e(4,4) = G;
		D_e(5,5) = G;
	}
}

template <class T_MATRIX>
inline void setB_Matrix(const size_t dim, double dshape_dx, double dshape_dy, double dshape_dz, T_MATRIX &B_matrix)
{
	switch(dim)
	{
	case 2:
		// B_11, dN/dx
		(*B_matrix)(0,0) = dshape_dx;
		// B_12, 0.0
		(*B_matrix)(0,1) = 0.0;
		// B_21, 0.0
		(*B_matrix)(1,0) = 0.0;
		// B_22, dN/dy
		(*B_matrix)(1,1) = dshape_dy;
		// B_31, 0.0
		(*B_matrix)(2,0) = 0.0;
		// B_32, 0.0
		(*B_matrix)(2,1) = 0.0;
		// B_41, dN/dy
		(*B_matrix)(3,0) = dshape_dy;
		// B_42, dN/dx
		(*B_matrix)(3,1) = dshape_dx;

		break;
	case 3:
		// B_11, dN/dx
		(*B_matrix)(0,0) = dshape_dx;
		// B_22, dN/dy
		(*B_matrix)(1,1) = dshape_dy;
		// B_33, dN/dz
		(*B_matrix)(2,2) = dshape_dz;
		//
		// B_41, dN/dy
		(*B_matrix)(3,0) = dshape_dy;
		// B_42, dN/dx
		(*B_matrix)(3,1) = dshape_dx;
		//
		// B_51, dN/dz
		(*B_matrix)(4,0) = dshape_dz;
		// B_53, dN/dx
		(*B_matrix)(4,2) = dshape_dx;
		//
		// B_62, dN/dz
		(*B_matrix)(5,1) = dshape_dz;
		// B_63, dN/dy
		(*B_matrix)(5,2) = dshape_dy;

		break;
	}
}

template <class T_MATRIX>
inline void setB_Matrix4axisymmetry(const size_t dim, double radius, double dshape_dx, double dshape_dy, T_MATRIX &B_matrix)
{
	// B_11, dN/dx
	(*B_matrix)(0,0) = dshape_dx;
	// B_12, 0.0
	(*B_matrix)(0,1) = 0.0;
	// B_21, N/r
	(*B_matrix)(1,0) = dshape_dx / radius;
	// B_22, 0.0
	(*B_matrix)(1,1) = 0.0;
	// B_31, 0.0
	(*B_matrix)(2,0) = 0.0;
	// B_32, dN/dz
	(*B_matrix)(2,1) = dshape_dy;
	// B_41, dN/dy
	(*B_matrix)(3,0) = dshape_dy;
	// B_42, dN/dx
	(*B_matrix)(3,1) = dshape_dx;
}
}
