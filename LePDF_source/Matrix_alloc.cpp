/*function to create a nxm matrix of zeros*/
double **zeros(int n, int m){
	int i,j;
	double **A;
	A = (double**)malloc(n*sizeof(double *));
	for (i=0;i<n;i++){
		A[i] = (double*)malloc(m*sizeof(double));
	}
	for (i=0;i<n;i++){
		for(j=0;j<m;j++){
		A[i][j]=0;
		}
	}
	return A;
}

/*function to free a matrix of doubles*/
void freemat(double **a, int n){
	int i;
	for(i=0;i<n;i++){
		free(a[i]);
	}
	free(a);
}
