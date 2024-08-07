#include "matrix.h"

matrix loadMatrix(int r, int c, vector<double> d){
    if (r < 0 || c < 0){
        return loadZero(0,0);
    }
    matrix A;
    A.rows = r;
    A.cols = c;
    A.REF = false;
    vector<vector<double>> init(r, vector<double>(c, 0.0));
    A.data = init;

    int size = r * c;

    while (d.size() < size)
    {
        d.push_back(0.0);
    }

    int k = 0;
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            A.data[i][j] = d[k];
            k++;
        }
    }

    return A;
}

matrix loadZero(int rows, int cols){
    if (rows < 0 || cols < 0){
        return loadZero(0,0);
    }
    matrix A;
    A.rows = rows;
    A.cols = cols;
    A.REF = false;
    vector<vector<double>> init(rows, vector<double>(cols, 0.0));
    A.data = init;
    return A;
}

void printMatrix(const matrix A){
    if (A.rows == 0 || A.cols == 0)
    {
        cout << "{}" << endl;
    }

    for (int i = 0; i < A.rows; i++)
    {
        if (i == 0)
        {
            cout << "{{";
        }
        else
        {
            cout << " {";
        }
        for (int j = 0; j < A.cols; j++)
        {
            if (j < (A.cols - 1))
            {
                cout << A.data[i][j] << ",\t";
            }
            else if (i == A.rows - 1)
            {
                cout << A.data[i][j] << "\t}}" << endl;
            }
            else
            {
                cout << A.data[i][j] << "\t}," << endl;
            }
        }
    }
}

void printMatrix(const matrix A, string name){
    cout << name << ":" << endl;
    if (A.rows == 0 || A.cols == 0)
    {
        cout << "{}" << endl;
    }

    for (int i = 0; i < A.rows; i++)
    {
        if (i == 0)
        {
            cout << "{{";
        }
        else
        {
            cout << " {";
        }
        for (int j = 0; j < A.cols; j++)
        {
            if (j < (A.cols - 1))
            {
                cout << A.data[i][j] << ",\t";
            }
            else if (i == A.rows - 1)
            {
                cout << A.data[i][j] << "\t}}" << endl;
            }
            else
            {
                cout << A.data[i][j] << "\t}," << endl;
            }
        }
    }

}

void printMatrixwithdeets(const matrix A){
    if (A.rows == 0 || A.cols == 0)
    {
        cout << "{}" << endl;
    }

    for (int i = 0; i < A.rows; i++)
    {
        if (i == 0)
        {
            cout << "{{";
        }
        else
        {
            cout << " {";
        }
        for (int j = 0; j < A.cols; j++)
        {
            if (j < (A.cols - 1))
            {
                cout << A.data[i][j] << ",\t";
            }
            else if (i == A.rows - 1)
            {
                cout << A.data[i][j] << "\t}}" << endl;
            }
            else
            {
                cout << A.data[i][j] << "\t}," << endl;
            }
        }
    }
    cout << "rows: " << A.rows << endl;
    cout << "cols: " << A.cols << endl;
    cout << "ref: " << A.REF << endl;
}

void printVector(const vector<double> v){
    if (v.empty())
    {
        cout << "{}" << endl;
        return;
    }
    cout << "{";
    for (int i = 0; i < v.size(); i++)
    {
        if (i == v.size() - 1)
        {
            cout << v[i] << "\t}" << endl;
        }
        else
        {
            cout << v[i] << ",\t";
        }
    }
}

bool equal(matrix A, matrix B){
    double EPS = 1e-10;
    if (A.rows != B.rows || A.cols != B.cols){
        return false;
    }

    for (int i=0; i<A.rows; i++){
        for (int j=0; j<A.cols; j++){
            double low = A.data[i][j] - EPS;
            double high = A.data[i][j] + EPS;
            // i googled for a function for determinint whether a value is within specific bounds and found clamp
            if (A.data[i][j] != clamp(B.data[i][j], low, high)){
                return false;
            }
        }
    }

    return true;
}

matrix scalarMultiplication(const matrix A, int s)
{
    matrix scaled = loadZero(A.rows, A.cols);
    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            scaled.data[i][j] = A.data[i][j] * s;
        }
    }
    return scaled;
}

matrix matrixMultiplication(const matrix A, const matrix B){
    if (A.cols != B.rows)
    {
        return loadZero(0, 0);
    }

    matrix AB = loadZero(A.rows, B.cols);
    vector<double> v1, v2;
    double dot;
    matrix Bt = transpose(B);
    for (int i = 0; i < AB.rows; i++)
    {
        for (int j = 0; j < AB.cols; j++)
        {
            // calculate and assign dot product
            dot = 0;
            for (int k = 0; k < A.data[i].size(); k++)
            {
                dot += (A.data[i][k] * Bt.data[j][k]);
            }
            AB.data[i][j] = dot;
        }
    }
    return AB;
}

matrix rowSwap(const matrix A, int r1, int r2){
    if (r1 > A.rows || r2 > A.rows || r1 < 1 || r2 < 1)
    {
        return loadZero(0, 0);
    }
    matrix swapped = A;
    vector<double> temp = A.data[r1 - 1];
    swapped.data[r1 - 1] = swapped.data[r2 - 1];
    swapped.data[r2 - 1] = temp;
    return swapped;
}

matrix transpose(const matrix A){
    matrix t = loadZero(A.cols, A.rows);
    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            t.data[j][i] = A.data[i][j];
        }
    }
    return t;
}

// anything not folded/ with a ; is incomplete

//TODO
/*
        - review gaussianelim
        - research on why gauss sucks and other better complexity ways of approaching these problems
        - write all test functions
        - write a script and record a video
        - send to dad

    done
        - rref function finalize
        - test inverse
        - matrix vector solve
        - lu fact


*/

matrix gaussianElimination(const matrix A){    
    double EPS = 1e-10; //stands for epsilon, this is a small constant value used to account for floating-point precision issues in numerical computations
    matrix gauss = A; // copy matrix A into gauss

    int j=0; // current column being processed

    for (int i=0; i<A.rows; i++) // iterate through the rows
    {
        // find a pivot for the row
        bool pivot_found = false;
        while (j < (A.cols - 1) && !pivot_found)
        {
            if (gauss.data[i][j] != 0) { // if pivot is not equal to 0 we found a pivot
                pivot_found = true;
            } else { // check for a possible swap
            // we swap for larger pivots for numerical stability
                int max_row = i;
                double max_val = 0;
                for (int k = i + 1; k < A.rows; ++k) 
                {
                    // cur_abs is the current absolute value of the element we are examining as a potentially better pivot
                    double cur_abs;
                    if (gauss.data[k][j] >= 0){
                        cur_abs = gauss.data[k][j];
                    } 
                    else{
                        cur_abs = -1 * gauss.data[k][j];
                    }

                    if (cur_abs > max_val)
                    {
                        max_row = k;
                        max_val = cur_abs;
                    }
                }
                if (max_row != i) {
                    gauss = rowSwap(gauss, max_row, i);
                    pivot_found = true;
                } else {
                    j++;
                }
            }
        }

        // perform elimination as normal if pivot was found
        if (pivot_found)
        {
            for (int t = i + 1; t < A.rows; ++t) {
                for (int s = j + 1; s < A.cols; ++s) {
                    gauss.data[t][s] = gauss.data[t][s] - gauss.data[i][s] * (gauss.data[t][j] / gauss.data[i][j]);

                    // the following step sets the value to zero if it is close to zero
                    // this helps with numerical stability
                    if (gauss.data[t][s] < EPS && gauss.data[t][s] > -1*EPS)
                        {gauss.data[t][s] = 0;}
                }
                gauss.data[t][j] = 0;
            }
        }

        j++;
    }
    gauss.REF = true;
    return gauss;
}

matrix reducedRowEchelon(const matrix A){
    double EPS = 1e-10;
    matrix rref = A;
    if (rref.REF == false){
        rref = gaussianElimination(rref);
    }
    for (int i=0; i<rref.rows; i++){
        if (rref.data[i][i] != 0){
            double factor = rref.data[i][i];
            for (int j=i; j<rref.cols; ++j){
                rref.data[i][j] = rref.data[i][j]/factor; //normalize pivot row
                if (rref.data[i][j] < EPS && rref.data[i][j] > -1*EPS)
                        {rref.data[i][j] = 0;}
            }
        
            for (int k=0; k<i; k++){
                double factor2 = rref.data[k][i]; 
                for (int h=0; h<rref.cols; h++){
                    rref.data[k][h] = rref.data[k][h] - (factor2*rref.data[i][h]);
                    if (rref.data[k][h] < EPS && rref.data[k][h] > -1*EPS)
                        {rref.data[k][h] = 0;}
                }
            }
        }
    }
    return rref;
}

matrix inverse(const matrix A){
    // check if A is square
    if (A.rows != A.cols){
        return loadZero(0,0);
    }
    int size = A.rows;

    matrix I = createIdentity(size);
    matrix AI = augmentMatrix(A, I);
    matrix U = gaussianElimination(AI);
    matrix IAinv = reducedRowEchelon(U);

    //printMatrix(IAinv, "RREF of A augmented with identity matrix");

    matrix check = loadZero(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            check.data[i][j] = IAinv.data[i][j];
        }
    }
    if (!equal(check, I)){
        return loadZero(0,0);
    }


    matrix Ainv = loadZero(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            Ainv.data[i][j] = IAinv.data[i][j+size];
        }
    }
    return Ainv;
}

matrix createIdentity(int size)
{
    if (size < 0){
        return loadZero(0,0);
    }
    matrix id = loadZero(size, size);
    for (int i = 0; i < id.rows; ++i) {
        for (int j = 0; j < id.cols; ++j) {
            if (i == j) {
                id.data[i][j] = 1;
            } else {
                id.data[i][j] = 0;
            }
        }
    }
    id.rows = size;
    id.cols = size;
    return id;
}

matrix augmentVector(const matrix A, const vector<double> b){
    if (A.rows != b.size()){
        return loadZero(0,0);
    }

    matrix aug = A;
    aug.cols = A.cols+1;
    for (int i=0; i<A.rows; i++){
        aug.data[i].push_back(b[i]);
    }
    return aug;
}

matrix augmentMatrix(const matrix A, const matrix B){
    if (A.rows != B.rows){
        return loadZero(0,0);
    }

    matrix aug = A;
    aug.cols = A.cols+B.cols;
    for (int i=0; i<A.rows; i++){
        for (int j=0; j<B.cols; j++){
            aug.data[i].push_back(B.data[i][j]);
        }
    }
    return aug;
}

vector<double> matrixVectorEqGauss(const matrix A, const vector<double> b){
    matrix aug = augmentVector(A,b);
    aug = reducedRowEchelon(aug);

    // check for inf or no solutions
    bool infSol = true;
    for (int i = 0; i < aug.cols - 1; i++) {
        if (aug.data[aug.rows - 1][i] != 0) {
            infSol = false;
            break;
        }
    }
    
    // check if the last row is all zeros except the last element
    if (infSol && aug.data[aug.rows - 1][aug.cols - 1] != 0) {
        //no solution exists
        return {};
    }

    // check if the last row is all zeros
    if (infSol && aug.data[aug.rows - 1][aug.cols - 1] == 0) {
        //infinite solutions exist
        return {};
    }
    // otherwise get last column
    vector<double> ans;
    int last = aug.cols - 1;
    for (int i=0; i<aug.rows; i++){
        ans.push_back(aug.data[i][last]);
    }

    return ans;
}

vector<double> matrixVectorEqInv(const matrix A, const vector<double> b){
    matrix bmatrix = loadMatrix(b.size(), 1, b);
    bmatrix = matrixMultiplication(inverse(A), bmatrix);
    if (bmatrix.data.empty()){
        return {};
    }
    return transpose(bmatrix).data[0];
}

vector<matrix> LUfact(const matrix A){
    if (A.rows != A.cols){
        return {loadZero(0,0), loadZero(0,0)};
    }
    int size = A.rows;
    matrix L = loadZero(size, size);
    matrix U = loadZero(size, size);

    // Decomposing matrix into Upper and Lower
    // triangular matrix
    for (int i = 0; i < size; i++) 
    {
        // Upper Triangular
        for (int k = i; k < size; k++)
        {
            // Summation of L(i, j) * U(j, k)
            int sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L.data[i][j] * U.data[j][k]);

            // Evaluating U(i, k)
            U.data[i][k] = A.data[i][k] - sum;
        }

        // Lower Triangular
        for (int k = i; k < size; k++) 
        {
            if (i == k)
                L.data[i][i] = 1; // Diagonal as 1
            else 
            {
                // Summation of L(k, j) * U(j, i)
                int sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L.data[k][j] * U.data[j][i]);

                // Evaluating L(k, i)
                L.data[k][i]
                    = (A.data[k][i] - sum) / U.data[i][i];
            }
        }
    }
    return {L, U};
}
