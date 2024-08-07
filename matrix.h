#ifndef MATRIX_H__
#define MATRIX_H__

#include <iostream>
//#include <map>
//#include <memory>
//#include <queue>
#include <string>
#include <vector>

using namespace std;

struct matrix
{
    int rows;
    int cols;
    int rank;
    bool REF; // is the matrix in row echelon form?
    vector<vector<double>> data;
    // { {a, b, c},
    //   {d, e, f},
    //   {g, h, i}}
};

// loads the matrix from left to right, top to bottom with the data in d[] following the matrix size rxc
// if vector input is too small this function will fill the rest of the matrix with zeros
// if vector input is too large this function will fill the matrix until capacity and ignore the rest
matrix loadMatrix(int rows, int cols, vector<double> d);

// initializes a zero matrix with size rows x cols
matrix loadZero(int rows, int cols);

// checks whether each element of matrices A and B are equal
bool equal(matrix A, matrix B);

// creates an identity matrix of given size
matrix createIdentity(int size);

// returns the input matrix with each element multiplied by s
matrix scalarMultiplication(const matrix A, int s);

// returns product of imput matricies if dimensions are compatible
matrix matrixMultiplication(const matrix first, const matrix second);

// returns a matrix with swapped rows r1 and r2
matrix rowSwap(const matrix A, int r1, int r2);

//returns the transpose of the input matrix
matrix transpose(const matrix A);
 

// extends matrix A rightwards with vector b if the size is compatible
matrix augmentVector(const matrix A, const vector<double> b);

// extends matrix A rightwards with matrix B if the size is compatible
matrix augmentMatrix(const matrix A, const matrix B);

//performs Gaussian elimination on matrix A to transform it into an
// upper triangular matrix (row echelon form)
// adapted from: https://github.com/akalicki/matrix/blob/master/dist/matrix.cpp#L8 
matrix gaussianElimination(const matrix A);

// from the Gaussian elimination row echelon matrix, this function continues the process to reduce the matrix to reduced row echelon form
matrix reducedRowEchelon(const matrix A);

// augments matrix A with the identity matrix of the same size and performs gaussian elimination to find the inverse of A
matrix inverse(const matrix A);
//no
// Solves a matrix vector equation of form Ax=b using Gaussian elimination
vector<double> matrixVectorEqGauss(const matrix A, const vector<double> b);

// Solves a matrix vector equation of form Ax=b using the inverse
vector<double> matrixVectorEqInv(const matrix A, const vector<double> b);

// Performs LU factorization
vector<matrix> LUfact(const matrix A);


//PRINTING FUNCTIONS
// prints matrix
void printMatrix(const matrix A);

// prints matrix with name
void printMatrix(const matrix A, string name);

// prints matrix with corresponding details
void printMatrixwithdeets(const matrix A);

// prints a vector
void printVector(const vector<double> v);

#endif // MATRIX_H__