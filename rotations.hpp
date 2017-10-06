#ifndef ROTATIONS_HPP
#define ROTATIONS_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>

//used for storing index of a matrix element which is to be killed at every iteration
struct index
{
    unsigned int row;
    unsigned int column;
} typedef index;

class Jacobi
{
private:
    unsigned int dimension;
    std::vector <double> matrix1;
    std::vector <double> matrix2;
    double mu, alpha, beta;
    double signum (double x)
    {
        return (x > 0) ? 1 : -1;
    }
    unsigned int up_index (unsigned int i, unsigned int j)
    {
        //U-matrixes stored as a linear array
        if (i > j)
        {
            unsigned int tmp = i;
            i = j;
            j = tmp;
        }
        unsigned int result = dimension * i + j - i * (i + 1) / 2;
        //std::cout<<"i = "<<i<<" j = "<<j<<" index is "<<result<<std::endl;
        //WARNING! Be careful using up_index and changing the matrix! 
        //Changing a(i,j) and than a(j,i) will result in an unwanted miscalculation!
        return result; 
    }
    index Find_max ();
    //Find_optimum func should be added here for better performance
public:
    Jacobi (unsigned int dim, std::vector <double> &mtrx);
    void Rotate ();
    void Print ();
};

Jacobi::Jacobi (unsigned int dim, std::vector <double> &mtrx)
{
    dimension = dim;
    if (dim == 0)
    {
        std::cout<<"Dimension must be positive"<<std::endl;
        exit (-1);
    }
    // checking whether the matrix is symmetric
    for (unsigned int i = 0; i < dim; i++)
        for (unsigned int j = 0; j < dim; j++)
            if ((mtrx[i * dimension + j] != mtrx[j * dimension + i]) && (i != j))
            {
                std::cout<<"The matrix is not symmetric"<<std::endl;
                exit (-1);
            }
    // NB: since the matrix is supposed to be symetric, it is enough to store only the upper-triangle and the diagonal elements.
    // The elements are stored in a linear array. 
    matrix1 = std::vector <double> (dimension * (dimension + 1) / 2);
    matrix2 = std::vector <double> (dimension * dimension);
    unsigned int k = 0;
    for (unsigned int i = 0; i < dimension; i ++)
        for (unsigned int j = i; j < dimension; j++)
        {
            matrix1 [k] = mtrx [i * dimension + j];
            k++;
        }
    mu = 0;
    alpha = 0;
    beta = 0;
}

index Jacobi::Find_max ()
{
    index result;
    result.row = 0;
    double cur_element;
    cur_element = (dimension == 1) ? matrix1[0] : matrix1[1];
    result.column = (dimension > 1) ? 1 : 0;
    for (unsigned int i = 0; i < dimension; i++)
        for (unsigned int j = 0; j < dimension; j++)
            if (i != j && matrix1 [up_index (i, j)] > cur_element)
            {
                cur_element = matrix1 [up_index (i, j)];
                //std::cout<<cur_element<<std::endl;
                result.row = i;
                result.column = j;
            }
    std::cout <<"max element is a "<<result.row<<" "<<result.column<<std::endl;
    return result;
}

void Jacobi::Rotate () // row == k, column == l
{
    index kl = Find_max ();
    mu = 2 * matrix1 [up_index (kl.row, kl.column)] / (matrix1 [up_index (kl.row, kl.row)] - matrix1 [up_index (kl.column, kl.column)]);
    alpha = sqrt ((1 + 1 / sqrt (1 + mu * mu)) / 2); 
    beta = signum (mu) * sqrt ((1 - 1 / sqrt (1 + mu * mu)) / 2);
    //std::cout <<"alpha is "<<alpha<<std::endl<<"beta is "<<beta<<std::endl<<"mu is "<<mu<<std::endl;
    for (unsigned int i = 0; i < dimension; i++)
    {
        //if (i <= kl.row)
            matrix2[i * dimension + kl.row]    = matrix1[up_index (i, kl.row)] *  alpha  + matrix1[up_index (i, kl.column)] * beta;
            //std::cout<<i<<" "<<kl.row<<" is "<<matrix2[i * dimension + kl.row]<<std::endl;
        //if (i <= kl.column)
            matrix2[i * dimension + kl.column] = matrix1[up_index (i, kl.row)] * (-beta) + matrix1[up_index (i, kl.column)] * alpha;
            //std::cout<<i<<" "<<kl.column<<" is "<<matrix2[i * dimension + kl.column]<<std::endl;
        for (unsigned int j = 0; j < dimension; j++)
            if (j != kl.row && j != kl.column)
                matrix2[i * dimension + j] = matrix1[up_index (i, j)];
    }
    for (unsigned int i = 0; i < dimension; i++)
    {
        if (i >= kl.row)
            matrix1[up_index (kl.row, i)]    = matrix2[kl.row * dimension + i] *  alpha  + matrix2[kl.column * dimension + i] * beta;
        if (i >= kl.column) 
            matrix1[up_index (kl.column, i)] = matrix2[kl.row * dimension + i] * (-beta) + matrix2[kl.column * dimension + i] * alpha;
        for (unsigned int j = 0; j < dimension; j++)
            if (j != kl.row && j!= kl.column && j <= i)
                matrix1[up_index (j, i)] = matrix2[j * dimension + i];
    }
}

void Jacobi::Print ()
{
    std::cout<<std::endl;
    for (unsigned int i = 0; i < dimension; i++)
    {
        for (unsigned int j = 0; j < dimension; j++)
            std::cout<<matrix1[up_index(i, j)]<<" ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

#endif
