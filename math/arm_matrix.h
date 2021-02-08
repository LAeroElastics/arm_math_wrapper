/*
 * arm_matrix.h
 *
 *  Created on: Feb 7, 2021
 *      Author:
 */

#ifndef ARM_MATRIX_H_
#define ARM_MATRIX_H_

#define ARM_MATH_CM7

#include "arm_math.h"
#include "arm_const_structs.h"

namespace math{

 template <unsigned int M, unsigned int N>
 class __EXPORT Matrix

 template <unsigned int M, unsigned int N>
 class __EXPORT MatrixBase{
 public:
	 /**
	  * matrix data storage
	  */
	 double data[M][N];

	 /**
	  * struct for usign arm_math
	  */
	 arm_matrix_instance_f64 arm_mat;

	 /**
	  * init to zero
	  */
	 MatrixBase() : data{}, arm_mat{M, N, &data[0][0]}{}
	 ~MatrixBase(){};

	 MatrixBase(const MatrixBase<M,N>& orig) : arm_mat{M, N, &data[0][0]}{
		 memcpy(data, orig.data, sizeof(data));
	 }

	 MatrixBase(const float* scaler) : arm_mat{M, N, &data[0][0]}{
		 memcpy(data, scaler, sizeof(data));
	 }

	 MatrixBase(const float carray[M][N]) : arm_mat{M, N, &data[0][0]}{
		 memcpy(data, carray, sizeof(data));
	 }

	 /**
	  * provide index access
	  */
	 double& operator()(const unsigned int& row, const unsigned int& col){
		 return (*this)[row][col];
	 }

	 unsigned int rows() const{return M;}
	 unsigned int cols() const{return N;}

	 bool operator==(const Matrix<M,N>& rhs) const{
		 for (unsigned int i(0); i < rows(); ++i){
			 for(unsigned int j(0); j < cols(); ++j){
				if(data[i][j] != rhs.data[i][j]){
					return false;
				}
			 }
		 }
		 return true;
	 }

	 bool operator!=(const Matrix<M,N>& rhs) const{
		 return !((*this) == rhs);
	 }

	 Matrix<M,N>& operator=(const Matrix<M,N>& rhs){
		 memcpy(data, rhs.data, sizeof(data));
		 return *this;
	 }

	 Matrix<M,N>& operator+=(const Matrix<M,N>& rhs){
		 for (unsigned int i(0); i < rows(); ++i){
		 		 		for(unsigned int j(0); j < cols(); ++j){
		 		 			data[i][j] += rhs.data[i][j];
		  	  }
		 }
		 return *this;
	 }

	 Matrix<M,N>& operator-=(const Matrix<M,N>& rhs){
		 for (unsigned int i(0); i < rows(); ++i){
		 		 		for(unsigned int j(0); j < cols(); ++j){
		 		 			data[i][j] -= rhs.data[i][j];
		  	  }
		 }
		 return *this;
	 }

	 Matrix<M,N>& operator*=(const double& scaler){
		 for (unsigned int i(0); i < rows(); ++i){
		 		 		for(unsigned int j(0); j < cols(); ++j){
		 		 			data[i][j] *= scaler;
		  	  }
		 }
		 return *this;
	 }

	 Matrix<M,N>& operator/=(const double& scaler){
		 for (unsigned int i(0); i < rows(); ++i){
		 		 		for(unsigned int j(0); j < cols(); ++j){
		 		 			data[i][j] /= scaler;
		  	  }
		 }
		 return *this;
	 }

	 Matrix<M,N> operator+(const Matrix<M,N>& rhs) const{
		 Matrix<M,N> result(*this);
		 result += rhs;
		 return result;
	 }

	 Matrix<M,N> operator-(const Matrix<M,N>& rhs) const{
		 Matrix<M,N> result(*this);
		 result -= rhs;
		 return result;
	 }

	 Matrix<M,N> operator*(const double& scaler) const{
		 Matrix<M,N> result(*this);
		 result *= scaler;
		 return result;
	 }

	 Matrix<M,N> operator/(const double& scaler) const{
		 Matrix<M,N> result(*this);
		 result /= scaler;
		 return result;
	 }

	 /**
	  * Matrix multiply (using CMSIS func)
	  */
	 Matrix<M,N> operator*(const Matrix<M,N>& mat) const{
		 Matrix<M,N> result;
		 arm_mat_mult_f32(&arm_mat, &mat.arm_mat, &result.arm_mat);
		 return result;
	 }

	 /**
	  * Matrix transpose (using CMSIS func)
	  */
	 Matrix<M,N> transpose() const{
		 Matrix<M,N> result;
		 arm_mat_trans_f32(&this->arm_mat, &result.arm_mat);
		 return result;
	 }

	 /**
	  * Matrix inverse (using CMSIS func)
	  */
	 Matrix<M,N> inverse() const{
		 Matrix<M,N> result;
		 arm_mat_inverse_f64(&this->arm_mat, &result.arm_mat);
		 return result;
	 }

	 // @todo: implement simple version(mult, trans, inv)

	 void zero(){
		 memset(data, 0, sizeof(data));
	 }

	 void eye(){
		 this->zero();
		 for(unsigned int i(0); i < M; ++i){
			 data[i][i] = 1.0;
		 }
	 }
 };

template <unsigned int M, unsigned int N>
class __EXPORT Matrix : public MatrixBase<M,N>{
public:
	//using MatrixBase<M,N>::operator*;
	Matrix() : MatrixBase<M,N>(){}
	Matrix(const Matrix<M,N>& orig) : MatrixBase<M,N>(orig){}
	Matrix(const double* scaler) : MatrixBase<M,N>(scaler){}
	Matrix(const double carray[M][N]) : MatrixBase<M,N>(carray){}

	Matrix<M,N>& operator=(const Matrix<M,N>& mat){
		memcpy(this->data, mat.data, sizoef(this->data));
		return *this;
	}

};

template<>
class __EXPORT Matrix<3,3> : public MatrixBase<3,3>{
	Matrix() : MatrixBase<3,3>(){}
	Matrix(const Matrix<3,3>& orig) : MatrixBase<3,3>(orig){}
	Matrix(const double* scaler) : MatrixBase<3,3>(scaler){}
	Matrix(const double carray[3][3]) : MatrixBase<3,3>(carray){}

	Matrix<3,3>& operator=(const Matrix<3,3>& mat){
		memcpy(this->data, mat.data, sizoef(this->data));
		return *this;
	}
};

template<>
class __EXPORT Matrix<4,4> : public MatrixBase<4,4>{
	Matrix() : MatrixBase<4,4>(){}
	Matrix(const Matrix<4,4>& orig) : MatrixBase<4,4>(orig){}
	Matrix(const double* scaler) : MatrixBase<4,4>(scaler){}
	Matrix(const double carray[4][4]) : MatrixBase<4,4>(carray){}

	Matrix<4,4>& operator=(const Matrix<4,4>& mat){
		memcpy(this->data, mat.data, sizoef(this->data));
		return *this;
	}
};

}




#endif /* ARM_MATRIX_H_ */
