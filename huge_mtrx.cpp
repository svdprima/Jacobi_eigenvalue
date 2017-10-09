#include <vector>
#include "rotations.hpp"
#include <cstdlib>

int main ()
{
	unsigned int n = 100;
	srand(time(NULL));
	std::vector<double> A = std::vector<double> (n * n);
	for (unsigned int i = 0; i < n; i++)
		for (unsigned int j = 0; j < n; j++)
		{
			if (i >= j)
			{
				A[i * n + j] = rand()%100;
				A[j * n + i] = A[i * n + j];
			}
		}
	Jacobi rot = Jacobi (n, A);
	while (1)
	{
		rot.Rotate ();    
        //rot.Print ();
        //std::cout<<"eq() "<<rot.eq()<<std::endl;
        if (rot.eq() == 1)
			break;
	}
	
	return 0;
}
