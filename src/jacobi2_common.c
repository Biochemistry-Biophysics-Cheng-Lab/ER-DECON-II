#include "jacobi2_common.h"

/*
*/




void ApplyComplexTFp ( double *X, const double *TF, long Length )
{
	
	long J;
	double re, im;
	
	for ( J = Length/2L - 1L; J > -1L; J-- )		{
		re    =  *X++;               
		im    =  *X--;                
		*X       =  re*(*TF++);
		*X++    -=  im*(*TF);
		*X       =  re*(*TF--);
		*X++    +=  im*(*TF++);
		TF++;
	}
	
	
}


/*
 */



void ApplyComplexTFpr ( double *X, const double *TF, long Length )
{
	
	long   J;
	double re, im;
	
	for ( J = Length/2L - 1L; J > -1L; J-- )		{
		re    =  *X++;               
		im    =  *X--;                
		*X       =  re*(*TF++);
		*X++    +=  im*(*TF);
		*X       =  -re*(*TF--);
		*X++    +=  im*(*TF++);
		TF++;
	}
	
	
	
}



