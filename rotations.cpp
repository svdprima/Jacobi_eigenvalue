#include "rotations.hpp"
#include <vector>
#include <fstream>

void Eigen (std::vector <double> &mtrx, unsigned int n, Jacobi &rot)
{
    //a proper stop condition should be added
    unsigned int steps = 9;
    for (unsigned int i = 0; i < steps; i++)
    {
        rot.Rotate ();    
        rot.Print ();
    }
}

int main ()
{
    unsigned int n = 0;
    double delta = 0;
    std::vector <double> A;
    std::ifstream input ("matrix.txt");
    if (input.is_open ())
    {
        input >> n;
        std::cout << n << std::endl;
        A = std::vector <double> (n * n);
        input >> delta;
        for (unsigned int i = 0; (i < n * n) && (!input.eof()); i++)
            input >> A [i];
    }
    else 
        std::cout<<"Unable to open a file"<<std::endl;
    Jacobi rot = Jacobi (n, A);
    rot.Print ();
    Eigen (A, n, rot);
    return 0;
}
