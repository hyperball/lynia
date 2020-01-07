#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MALLOC(pointer, space_element)\
		if(!(pointer = malloc(space_element))){\
			fprintf(stderr, "Insufficient memory");\
			exit(EXIT_FAILURE);\
		}


typedef struct vector{
	int length;
	double *ptr_val;
}vector;

typedef struct matrix{
	int row_length;
	int col_length;
	double **ptr_entry;
}matrix;

 
double infinity_norm(const struct vector *ptr_vector);
struct vector *set_vector(int n, struct vector *ptr_vector);
struct matrix * set_matrix(int r, int c, struct matrix *ptr_matrix);
void get_vandermonde_matrix(struct matrix *ptr_matrix);

void get_vector_b(struct vector *ptr_vector, int n);
struct vector *equal_vector(const struct vector *ptr_vector, struct vector *ptr_res);
struct vector *add_vector(
	struct vector *ptr_vector_a, 
	struct vector *ptr_vector_b, 
	struct vector *ptr_res);
struct vector *subtract_vector(
	struct vector *ptr_vector_a, 
	struct vector *ptr_vector_b,
	struct vector *ptr_res);
struct vector *linear_solver(
	const struct matrix *ptr_M, 
	const struct vector *ptr_RHS, 
	struct vector *ptr_x);

void print_matrix(const struct matrix *ptr_matrix);
void print_vector(const struct vector *ptr_vector);

int main(int argc, char **argv)
{
	int n, r, c;
	int i, j;
	scanf("%d", &n);
	struct vector b;
	struct vector x;
	set_vector(n, &b);
	for(i = 0; i < n; i++){
		scanf("%lf", &(b.ptr_val[i]));	
	}
	set_vector(n, &x);
	struct matrix A;
	scanf("%d %d", &r, &c);
	set_matrix(r, c, &A);
	
	for(i = 0; i < r; i++){
		for(j = 0; j < c; j ++){
			scanf("%lf", &(A.ptr_entry[i][j]));
		}	
	}
	linear_solver(&A, &b, &x);
 	print_vector(&x);
	
	return(0);	
}



double infinity_norm(const struct vector *ptr_vector)
{
	int n = ptr_vector->length;
	double max = 0.0;
	int i;
	for(i = 0; i < n; i++){
		if(max < fabs(ptr_vector->ptr_val[i])){
			max = fabs(ptr_vector->ptr_val[i]);
		}
	}
	return(max);
}

struct vector *set_vector(int n, struct vector *ptr_vector)
{
	ptr_vector->length = n;
	MALLOC(ptr_vector->ptr_val, n * sizeof(*(ptr_vector->ptr_val)))
	return(ptr_vector);
}

struct matrix *set_matrix(int r, int c, struct matrix *ptr_matrix)
{
	ptr_matrix->row_length = r;
	ptr_matrix->col_length = c;
	MALLOC(ptr_matrix ->ptr_entry, c * sizeof(*(ptr_matrix->ptr_entry)))
	int i;
	for(i = 0; i < c; i++){
		MALLOC(ptr_matrix->ptr_entry[i], r * sizeof(*(*(ptr_matrix->ptr_entry))))
	}
	return(ptr_matrix);
}


struct vector *equal_vector(const struct vector * ptr_vector, struct vector *ptr_res){
	int n = ptr_vector->length;
	int i;
	for(i = 0; i < n; i++){
		ptr_res->ptr_val[i] = ptr_vector->ptr_val[i];
	}
	return(ptr_res);
}

struct vector *add_vector(
	struct vector *ptr_vector_a, 
	struct vector *ptr_vector_b, 
	struct vector *ptr_res)
{
	int n = ptr_vector_a->length;
	int i;
	for(i = 0; i < n; i++){
		ptr_res->ptr_val[i] = ptr_vector_a->ptr_val[i] + ptr_vector_b->ptr_val[i];
	}
	return(ptr_res);
}

struct vector *subtract_vector(
	struct vector *ptr_vector_a, 
	struct vector *ptr_vector_b, 
	struct vector *ptr_res)
{
	int n = ptr_vector_a->length;
	int i;
	for(i = 0; i < n; i++){
		ptr_res->ptr_val[i] = ptr_vector_a->ptr_val[i] - ptr_vector_b->ptr_val[i];
	}
	return(ptr_res);
}

void print_vector(const vector *ptr_vector)
{
	int n = ptr_vector->length;
	int i;
	for(i = 0; i < n; i++){
		printf("%-16.8lf\n", ptr_vector->ptr_val[i]);
	}
}


struct vector *linear_solver(
	const struct matrix *ptr_M, 
	const struct vector *ptr_b,
	struct vector *ptr_x)
{
	double tolerance = 1e-12;
	int  dimension;
	dimension = ptr_M->row_length;
	int i, j;

	struct vector xk;
	set_vector(dimension, &xk);
	struct vector tmp_vec;
	set_vector(dimension, &tmp_vec);
	struct vector record_vec;
	set_vector(dimension, &record_vec);
	do{		 
		equal_vector(ptr_x, &record_vec);
		for(i = 0; i < dimension; i++){
			double lx = 0.0;
			double ux = 0.0;
			for(j = i + 1; j < dimension; j++){
				ux += ptr_M->ptr_entry[i][j] * ptr_x->ptr_val[j];
			}
			for(j = 0; j < i; j ++){
				lx += ptr_M->ptr_entry[i][j] * ptr_x->ptr_val[j];	
			}
			xk.ptr_val[i] = (ptr_b->ptr_val[i] - lx - ux)
						  / (ptr_M->ptr_entry[i][i]);
			ptr_x->ptr_val[i] = xk.ptr_val[i];				/*	Update the guess_solution	*/
		}
		subtract_vector(ptr_x, &record_vec, &tmp_vec);
	}while(fabs(infinity_norm(&tmp_vec)) > tolerance);
	return(ptr_x);
}

