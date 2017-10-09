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
    std::vector <double> eig_pre;
    std::vector <double> eig_post;
    double mu, alpha, beta;
    double signum (double x)
    {
        return (x > 0) ? 1 : -1;
    }
    double abs (double x)
    {
        return (x > 0) ? x : -x;
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
    index Find_optimum ();//Find_optimum func should be added here for better performance
public:
    Jacobi (unsigned int dim, std::vector <double> &mtrx);
    void Rotate ();
    void Print ();
    int eq ();
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
    eig_pre = std::vector <double> (dimension);
    eig_post = std::vector <double> (dimension);
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
    cur_element = (dimension == 1) ? abs(matrix1[0]) : abs(matrix1[1]);
    result.column = (dimension > 1) ? 1 : 0;
    for (unsigned int i = 0; i < dimension; i++)
        for (unsigned int j = 0; j < dimension; j++)
            if (i != j && abs (matrix1 [up_index (i, j)]) > cur_element)
            {
                cur_element = abs (matrix1 [up_index (i, j)]);
                //std::cout<<cur_element<<std::endl;
                result.row = i;
                result.column = j;
            }
    //std::cout <<"max element is a "<<result.row<<" "<<result.column<<std::endl;
    return result;
}

index Jacobi::Find_optimum ()
{
	index result;
	result.row = 0;
	double cur_element;
	unsigned int max_row = 0;
	double max_row_sum = 0;
	double cur_row_sum;
	cur_element = (dimension == 1) ? abs(matrix1[0]) : abs(matrix1[1]);
	result.column = (dimension > 1) ? 1 : 0;
	for (unsigned int i = 0; i < dimension; i++)
	{
		cur_row_sum = 0;
        for (unsigned int j = 0; j < dimension; j++)
		{
			if (i != j)
				cur_row_sum += abs (matrix1 [up_index (i, j)]);
		}
		if (cur_row_sum > max_row_sum)
		{
			max_row_sum = cur_row_sum;
			max_row = i;
		}
	}	
	result.row = max_row;
	for (unsigned int j = 0; j < dimension; j++)
	{
		if (max_row != j && abs (matrix1 [up_index (max_row, j)]) >= cur_element)
		{
			cur_element = abs (matrix1 [up_index (max_row, j)]);
			result.column = j;
		}
	}
	//std::cout <<"opt element is "<<result.row+1<<" "<<result.column+1<<std::endl;
    return result;
}


void Jacobi::Rotate () // row == k, column == l
{
	for (unsigned int i = 0; i < dimension; i++)
    {
		eig_pre [i] = matrix1 [up_index(i, i)];
		//std::cout<<"eig_pre "<<eig_pre[i]<<std::endl;
	}
	
    //index kl = Find_max ();
    index kl = Find_optimum ();
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
    for (unsigned int i = 0; i < dimension; i++)
    {
		eig_post [i] = matrix1 [up_index(i, i)];
		//std::cout<<"eig_post "<<eig_post[i]<<std::endl;
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

int Jacobi::eq()
{
	unsigned int g = 0;
	for (unsigned int i = 0; i < dimension; i++)
	{
		if (eig_pre [i] == eig_post [i])
			g++;
	}
	if (g == dimension)
	{
		for (unsigned int i = 0; i < dimension; i++)
		{
			std::cout<<"eig_post"<<eig_post[i]<<std::endl;
		}
		return 1;
	}
	else
		return 0;
}
		

#endif
