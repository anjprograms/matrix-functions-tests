#include "matrix.h"

// initialize tests for each function here
// each one should tell a story of how the function works
void TEST_loadFuncs();
void TEST_equal();
void TEST_multiplication();
void TEST_rowSwap();
void TEST_transpose();
void TEST_augment();
void TEST_gaussElim();
void TEST_rref();
void TEST_inverse();
void TEST_matrixVecSolve();
void TEST_LUfact();
string printBool(bool);

int main()
{
    //uncomment test to be performed

    //TEST_loadFuncs();
    //TEST_equal();
    //TEST_multiplication();
    //TEST_rowSwap();
    //TEST_transpose();
    //TEST_augment();
    //TEST_gaussElim();
    TEST_rref();
    //TEST_inverse();
    //TEST_matrixVecSolve();
    //TEST_LUfact();
}

void TEST_loadFuncs(){
cout << "______________________________________________________" << endl;
    cout << "matrix loadMatrix(int rows, int cols, vector<double> d)" << endl << "matrix loadZero(int rows, int cols)" << endl << "matrix createIdentity(int size)" << endl << endl;

    cout << "loadMatrix should make a matrix and fill it with given elements, it should fill not given elements with zeros and ignore extra elements given" << endl;
    printMatrix(loadMatrix(2, 2, {1, 2, 3, 4}) ,"loadMatrix(2, 2, {1, 2, 3, 4})");
    cout << endl;
    printMatrix(loadMatrix(2, 2, {2, 2, 9.01, 4, 8}), "loadMatrix(2, 2, {2, 2, 9.01, 4, 8})");
    printMatrix(loadMatrix(-1, 2, {2, 2, 9.01, 4, 8}), "loadMatrix(-1, 2, {2, 2, 9.01, 4, 8})");
    cout << endl;
    printMatrix(loadMatrix(2, 2, {1}) ,"loadMatrix(2, 2, {1})");
    cout << endl;
    printMatrix(loadMatrix(0, 0, {}) ,"loadMatrix(0, 0, {})");
    cout << endl << endl;
    printMatrix(loadZero(4, 2), "loadZero(4, 2)");
    cout << endl;
    printMatrix(loadZero(0,0) ,"loadZero(0,0)");
    cout << endl;
    printMatrix(createIdentity(3), "createIdentity(3)");
    cout << endl;
    printMatrix(createIdentity(0), "createIdentity(3)");
    cout << endl;
    cout << "______________________________________________________" << endl;
}

void TEST_equal(){
    cout << "______________________________________________________" << endl;
    cout << "bool equal(matrix A, matrix B)" << endl << endl;

    matrix A = createIdentity(3);
    matrix B = createIdentity(3);

    printMatrix(A ,"A");
    cout << endl;
    printMatrix(B, "B");
    cout << endl;
    cout << "equal(A, B): " << printBool(equal(A, B)) << endl << endl;

    B.data[2][2] = 1.0000000000000001;
    cout << "Set B.data[2][2] = 1.0000000000000001" << endl;
    cout << "equal(A, B): " << printBool(equal(A, B)) << endl << endl;

    B.data[2][2] = 1.0000000001;
    cout << "Set B.data[2][2] = 1.0000000001" << endl;
    cout << "equal(A, B): " << printBool(equal(A, B)) << endl;
    cout << "______________________________________________________" << endl;
}

void TEST_multiplication(){
    cout << "______________________________________________________" << endl;
    cout << "matrix scalarMultiplication(const matrix A, int s)" << endl << "matrix matrixMultiplication(const matrix first, const matrix second)" << endl << endl;

    matrix A = loadMatrix(2,3,{1, 2, 0, 6, 1});
    matrix B = loadMatrix(3,2,{1, 5, 2, 3, 2});

    printMatrix(A, "A");
    cout << endl;
    printMatrix(B, "B");
    cout << endl;

    printMatrix(scalarMultiplication(A, 2) ,"scalarMultiplication(A, 2)");
    cout << endl;
    printMatrix(matrixMultiplication(A, B), "matrixMultiplication(A, B)");
    printMatrix(matrixMultiplication(B, B), "matrixMultiplication(B, B)");
    cout << "______________________________________________________" << endl;
}

void TEST_rowSwap(){
cout << "______________________________________________________" << endl;
    cout << "matrix rowSwap(const matrix A, int r1, int r2)" << endl << endl;

    matrix A = loadMatrix(4,2,{1, 5, 2, 3, 9, 4, 6, 8});

    printMatrix(A, "A");
    cout << endl;

    printMatrix(rowSwap(A, 2, 1) ,"rowSwap(A, 2, 1)");
    cout << endl;
    printMatrix(rowSwap(A, 0, 1), "rowSwap(A, 0, 1)");
    cout << endl;
    printMatrix(rowSwap(A, -1, 1), "rowSwap(A, -1, 1)");
    cout << endl;
    printMatrix(rowSwap(A, 3, 8), "rowSwap(A, 3, 8)");    
    cout << "______________________________________________________" << endl;
}

void TEST_transpose(){
cout << "______________________________________________________" << endl;
    cout << "matrix transpose(const matrix A)" << endl << endl;

    matrix A = loadMatrix(4,2,{1, 5, 2, 3, 9, 4, 6, 8});
    matrix empty = loadZero(0,0);

    printMatrix(A, "A");
    cout << endl;
    printMatrix(empty, "empty");
    cout << endl;

    printMatrix(transpose(A) ,"transpose(A)");
    cout << endl;
    printMatrix(transpose(empty), "transpose(empty)");
    cout << endl;
    cout << "______________________________________________________" << endl;
}

void TEST_augment(){
    cout << "______________________________________________________" << endl;
    cout << "matrix augmentVector(const matrix A, const vector<double> b)" << endl << "matrix augmentMatrix(const matrix A, const matrix B)" << endl << endl;

    matrix A = loadMatrix(3,2,{1, 5, 2, 3, 2});
    matrix B = loadMatrix(3,3,{1, 2, 1 ,4, 5, 2, 1, 0, 2});
    matrix C = loadMatrix(4,2,{1, 5, 2, 3, 9, 4, 6, 8});

    vector<double> v1 = {2,4,5};
    vector<double> v2 = {1};

    printMatrix(A, "A");
    cout << endl;
    printMatrix(B, "B");
    cout << endl;
    printMatrix(C, "C");
    cout << endl;
    cout << "v1: ";
    printVector(v1);
    cout << endl;
    cout << "v2: ";
    printVector(v2);
    cout << endl;

    printMatrix(augmentVector(A, v1),"augmentVector(A, v1)");
    cout << endl;
    printMatrix(augmentVector(A, v2), "augmentVector(A, v2)");
    cout << endl;
    printMatrix(augmentMatrix(A, B), "augmentMatrix(A, B)");
    cout << endl;
    printMatrix(augmentMatrix(A, C), "augmentMatrix(A, C)");    
    cout << "______________________________________________________" << endl;
}

void TEST_gaussElim(){
cout << "______________________________________________________" << endl;
    cout << "matrix gaussianElimination(const matrix A)" << endl << endl;

    matrix tall = loadMatrix(4,2,{1, 5, 2, 3, 9, 4, 6, 8});
    matrix noSol = loadMatrix(3, 4, {1, 1, 1, 2,0,1,-3,1,2, 1, 5, 0});
    matrix infSol = loadMatrix(3, 4, {4,-3,1,-8, -2, 1, -3, -4, 2,-1,3,4});
    matrix yesSol = loadMatrix(3,4, {4, -3, 1, -8, -2, 1, -3, -4, 1, -1, 2, 2});
    matrix empty = loadZero(0,0);

    printMatrix(tall, "tall");
    cout << endl;
    printMatrix(noSol, "noSol");
    cout << endl;
    printMatrix(infSol, "infSol");
    cout << endl;
    printMatrix(yesSol, "yesSol");
    cout << endl;
    printMatrix(empty, "empty");
    cout << endl;

    printMatrix(gaussianElimination(tall) ,"gaussianElimination(tall)");
    cout << endl;
    printMatrix(gaussianElimination(noSol) ,"gaussianElimination(noSol)");
    cout << endl;
    printMatrix(gaussianElimination(infSol) ,"gaussianElimination(infSol)");
    cout << endl;
    printMatrix(gaussianElimination(yesSol) ,"gaussianElimination(yesSol)");
    cout << endl;
    printMatrix(gaussianElimination(empty) ,"gaussianElimination(empty)");
    cout << endl;  
    cout << "______________________________________________________" << endl;
}

void TEST_rref(){
cout << "______________________________________________________" << endl;
    cout << "matrix reducedRowEchelon(const matrix A)" << endl << endl;

    matrix tall = loadMatrix(4,2,{1, 5, 2, 3, 9, 4, 6, 8});
    matrix noSol = loadMatrix(3, 4, {1, 1, 1, 2,0,1,-3,1,2, 1, 5, 0});
    matrix infSol = loadMatrix(3, 4, {4,-3,1,-8, -2, 1, -3, -4, 2,-1,3,4});
    matrix yesSol = loadMatrix(3,4, {4, -3, 1, -8, -2, 1, -3, -4, 1, -1, 2, 2});
    matrix empty = loadZero(0,0);

    printMatrix(tall, "tall");
    cout << endl;
    printMatrix(noSol, "noSol");
    cout << endl;
    printMatrix(infSol, "infSol");
    cout << endl;
    printMatrix(yesSol, "yesSol");
    cout << endl;
    printMatrix(empty, "empty");
    cout << endl;

    printMatrix(reducedRowEchelon(tall) ,"reducedRowEchelon(tall)");
    cout << endl;
    printMatrix(reducedRowEchelon(noSol) ,"reducedRowEchelon(noSol)");
    cout << endl;
    printMatrix(reducedRowEchelon(infSol) ,"reducedRowEchelon(infSol)");
    cout << endl;
    printMatrix(reducedRowEchelon(yesSol) ,"reducedRowEchelon(yesSol)");
    cout << endl;
    printMatrix(reducedRowEchelon(empty) ,"reducedRowEchelon(empty)");
    cout << endl;  
    cout << "______________________________________________________" << endl;
}

void TEST_inverse(){
    cout << "______________________________________________________" << endl;
    cout << "matrix inverse(const matrix A)" << endl << endl;

    matrix notSquare = loadMatrix(2,3,{1, 2, 0, 6, 1});
    matrix noInv = loadMatrix(3, 3, {4,-3,1, -2, 1, -3, 2, -1, 3});
    matrix yesInv = loadMatrix(3, 3, {1, 1, 3,0,1,-3,2, 1, 5});

    printMatrix(notSquare ,"notSquare");
    cout << endl;
    printMatrix(noInv, "noInverse");
    cout << endl;
    printMatrix(yesInv, "yesInverse");
    cout << endl << endl;
    printMatrix(inverse(notSquare), "Inverse of notSquare");
    cout << endl << endl;
    printMatrix(inverse(noInv), "Inverse of noInv");
    cout << endl << endl;
    printMatrix(inverse(yesInv), "Inverse of yesInverse");
    cout << "______________________________________________________" << endl;
}

void TEST_matrixVecSolve(){
    cout << "______________________________________________________" << endl;
    cout << "vector<double> matrixVectorEqGauss(const matrix A, const vector<double> b)" << endl << "vector<double> matrixVectorEqInv(const matrix A, const vector<double> b)" << endl << endl;

    matrix noSol = loadMatrix(3, 4, {1, 1, 1, 2,0,1,-3,1,2, 1, 5, 0});
    matrix infSol = loadMatrix(3, 4, {4,-3,1,-8, -2, 1, -3, -4, 2,-1,3,4});
    matrix yesSol = loadMatrix(3,4, {4, -3, 1, -8, -2, 1, -3, -4, 1, -1, 2, 2});
    
    matrix noSol_A = loadMatrix(3, 3, {1, 1, 1,0,1,-3,2, 1, 5});
    vector<double> noSol_b = {2, 1, 0};
    matrix infSol_A = loadMatrix(3, 3, {4,-3,1, -2, 1, -3,  2,-1,3});
    vector<double> infSol_b = {-8,-4,4};
    matrix yesSol_A = loadMatrix(3,3, {4, -3, 1, -2, 1, -3, 1, -1, 2});
    vector<double> yesSol_b = {-8, -4, 2};

    printMatrix(noSol ,"noSol system");
    cout << endl;
    printMatrix(infSol, "infSol system");
    cout << endl;
    printMatrix(yesSol, "yesSol system");
    cout << endl;

    cout << "matrixVectorEqGauss(noSol): ";
    printVector(matrixVectorEqGauss(noSol_A, noSol_b));
    cout << endl;
    cout << "matrixVectorEqInv(noSol): ";
    printVector(matrixVectorEqInv(noSol_A, noSol_b));
    cout << endl << endl;

    cout << "matrixVectorEqGauss(infSol): ";
    printVector(matrixVectorEqGauss(infSol_A, infSol_b));
    cout << endl;
    cout << "matrixVectorEqInv(infSol): ";
    printVector(matrixVectorEqInv(infSol_A, infSol_b));
    cout << endl << endl;

    cout << "matrixVectorEqGauss(yesSol): ";
    printVector(matrixVectorEqGauss(yesSol_A, yesSol_b));
    cout << endl;
    cout << "matrixVectorEqInv(yesSol): ";
    printVector(matrixVectorEqInv(yesSol_A, yesSol_b));
    cout << "______________________________________________________" << endl;
}

// not done
void TEST_LUfact(){
    cout << "______________________________________________________" << endl;
    cout << "vector<matrix> LUfact(const matrix A)" << endl << endl;

    matrix A = loadMatrix(3, 3, {1, 1, 0,2, 1, 3,3, 1, 1});
    matrix L,U;
    vector<matrix> v;
    v = LUfact(A);
    L = v[0];
    U = v[1];

    printMatrix(A, "A");
    printMatrix(L, "L");
    printMatrix(U, "U");
    printMatrix(matrixMultiplication(L,U), "L*U");

}


string printBool(bool b){
    if (b == true){
        return "true";
    }
    else{
        return "false";
    }
}