

# Matrix Functions

Anjali Velamala
CSPB 2270-003: Data Structures & Algorithms
Final Project : Implementing Matrix Functionalities


**Compile:** g++ -Wall -Werror -Wpedantic -std=c++17 main.cpp matrix.cpp

**Run:** ./a.out

Video: https://drive.google.com/file/d/1YfAE6K88uH3V6ZOpKPmpg78ddHMexBon/view?usp=share_link


## Motivation
As I prepare to delve deeper into machine learning through upcoming courses in computer science and applied mathematics, I want to solidify my understanding of fundamental linear algebra concepts. Matrices are a crucial component in various machine learning algorithms, from simple linear regression to more complex neural networks. By implementing a matrix class from scratch, I will gain a more comprehensive understanding of matrix operations and their underlying mechanics. This project not only serves as a review but also as a practical exercise to bridge theoretical knowledge with practical implementation. I hope it will also serve as a valuable tool for future assignments.

## Project Scope and Goals
The primary goal of this project is to create a group of functions that can perform a variety of matrix operations on any m x n matrix (given the size is correct for the operation). The functions are:

* Scalar Multiplication
* Matrix Multiplication
* Row Swaps
* Matrix transpose
* Gaussian Elimination and Row Echelon Form
* Matrix Inverse
* Solving a Matrix-Vector Equation
* LU factorization

## Why This Is an Interesting Application
Matrix operations are at the heart of numerous algorithms in machine learning and data science. By implementing these operations from scratch, I will gain a deeper understanding of their computational complexity and efficiency. This understanding is critical as it influences algorithm design and optimization, particularly when working with large datasets. Moreover, this project serves as an intersection between theoretical linear algebra and practical programming, providing a tangible way to apply mathematical concepts.

## Technical Challenges and Learning Opportunities
**Efficiency and Optimization**
I used Gaussian elimination as a basis for all of my other more complicated functions, this algorithm is $O(n^3)$ so my functions are not efficient especially for large matrices.

The idea behind solving using any factorization method is to decouple the factorization phase from the actual solving phase because the factorization phase is usually more complex. The factorization phase only needs matrix $A$, while the actual solving phase uses the factored form of $A$ and the right hand side to solve the linear system. So, once we have the factorization, we can make use of the factored form of $A$, to solve for different right hand sides at a relatively moderate computational cost. The cost of factorizing the matrix 𝐴 into $LU$ is $O(n^3)$. With this factorization, the cost of solving is just $O(n^2)$.


**Numerical Stability**
When applicable I checked for division by zero and I also used a small epsilon value to account for floating-point precision issues in numerical computations performed. Setting a value to zero if val < |epsilon|. 
I also applied pivoting methods when useful. In my Gaussian elimination function the pivot was always the largest. I wasn't able to finish the pivoting for my LU factorization function.


## Resources used

**Reviewing Matrix Operations**

* https://vismor.com/documents/network_analysis/matrix_algorithms/matrix_algorithms.php
* https://math.stackexchange.com/questions/266355/necessity-advantage-of-lu-decomposition-over-gaussian-elimination 

**Help structuring/implementing code**

* My father who is a computer programmer
* My brother who just graduated with a Bachelor's in CS
* https://github.com/akalicki/matrix/blob/master/dist/matrix.cpp#L8 
* https://www.geeksforgeeks.org/doolittle-algorithm-lu-decomposition/