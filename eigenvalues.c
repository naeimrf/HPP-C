/* Naeim Rashidfarokhi
   Uppsala university 2019
   Individual HPP project*/

#include <math.h>  // rand()
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>  // calloc()
#include <string.h>  // memset()
#include <time.h>    // time() and clock_gettime()

static inline float** read_matrix_file(char* __restrict file_name, const int n);
static inline float** create_2d_matrix(const int n);
static inline void save_files(float* __restrict max, float* __restrict vector, float** __restrict matrix, const int n);

void power_method(float** __restrict matrix, const int iter, const float tolerance, const int n);
static inline float largest_absolute(float* __restrict vector, const int n);
static inline float* normalize_vector(float* __restrict vector, float max, const int n);

static inline void print_2d_matrix(float** __restrict matrix, const int n);  // could be constant function
static inline void print_vector(float* __restrict vector, const int n);
static inline double measure_time();  //based on clock_gettime()

#define convergence_approach 0

int main(int argc, char* argv[]) {
    double time0 = measure_time();

    if (argc < 5 || argc > 6) {
        printf(">> %d arguments entered!\n", argc - 1);
        printf(">> The program should have 4 or 5 input arguments!\n");
        printf("1. an integer as 'n' for size of square matrix!\n");
        printf("2. an integer for number of iter!\n");
        printf("3. a real number for tolerance!\n");
        printf("4. an integer for number of threads!\n");
        printf("5. an optional text-file for input data\n");
        return 1;
    }
    float** matrix_2d_A = NULL;
    //printf("Available threads on this machine:%d\n", omp_get_max_threads());

    const int n = atoi(argv[1]);
    const int n_itr = atoi(argv[2]);
    const float tol = atof(argv[3]);
    omp_set_num_threads(atoi(argv[4]));

    double time1 = measure_time();
    if (argv[5] != NULL)
        matrix_2d_A = read_matrix_file(argv[5], n);
    else
        matrix_2d_A = create_2d_matrix(n);
    double time2 = measure_time();

    if (n <= 2) {
        print_2d_matrix(matrix_2d_A, n);
    }

    power_method(matrix_2d_A, n_itr, tol, n);
    double time3 = measure_time();

    printf("Time to create matrix A: %.4fs\n", time2 - time1);
    printf("Power_method execution: %.4fs\n", time3 - time2);
    printf("Total execution time: %.4fs\n\n", measure_time() - time0);
    return 0;
}

/* modifed code from Lab1, task3 */
static inline float** create_2d_matrix(const int n) {
    srand(time(NULL));

    /* create and initialize a square matrix with n x n elements */
    float** matrix = (float**)malloc(n * sizeof(float*));
    for (int i = 0; i < n; ++i) {
        matrix[i] = (float*)malloc(n * sizeof(float));  // *matrix[i] same as *(matrix + i)
        for (int j = 0; j < n; j++) {
            matrix[i][j] = rand() % 20 - 10;  // random numbers ranging from -10 to 10.
        }
    }
    printf("Random Matrix: %dx%d created!\n", n, n);
    printf("X*****************************X\n");
    return matrix;
}

static inline void print_2d_matrix(float** __restrict matrix, const int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.0f ", matrix[i][j]);
        putchar('\n');
    }
}

static inline void print_vector(float* __restrict vector, const int n) {
    for (int i = 0; i < n; i++) {
        printf("%5.3f", vector[i]);
        putchar('\n');
    }
}

static inline double measure_time() {
    struct timespec t0;
    clock_gettime(CLOCK_REALTIME, &t0);
    double time_sec = t0.tv_sec + t0.tv_nsec / 1e9;
    return time_sec;
}

static inline void save_files(float* __restrict max, float* __restrict vector, float** __restrict matrix, const int n) {
    FILE* result = NULL;   // to save eigenvalues
    FILE* problem = NULL;  // to save main matrix
    char text0[] = "original matrix:\n";
    char text1[] = "eigenvalue:\n";
    char text2[] = "eigenvector:\n";

    /* to write original matrix to the file */
    if ((problem = fopen("input.txt", "w+")) == NULL) {
        printf("Error: with opening the file in save_files function!\n");
        exit(1);
    }
    if (fprintf(problem, "%s", text0) < 0) {  //fwrite(text0, 1, sizeof(text0), result) != sizeof(text0)
        printf("Error: writing text0 to file in save_files function\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fprintf(problem, "%f ", matrix[i][j]) < 0) {  //fwrite(row, sizeof(float), n, result) != n
                printf("Error: writing original matrix:row[%d] to file in save_files function\n", i);
                exit(1);
            }
        }
        fprintf(problem, "\n");
    }
    fclose(problem);

    /* to write the eigenvalue to the file */
    if ((result = fopen("output.txt", "w+")) == NULL) {
        printf("Error: with opening the file in save_files function!\n");
        exit(1);
    }
    if (fprintf(result, "%s", text1) < 0) {  //fwrite(text1, 1, sizeof(text1), result) != sizeof(text1)
        printf("Error: writing text1 to file in save_files function\n");
        exit(1);
    }
    if (fprintf(result, "%f\n", *max) < 0) {  //fwrite(&max, sizeof(float), 1, result) != 1
        printf("Error: writing eigenvalue to file in save_files function\n");
        exit(1);
    }

    /* to write the eigenvector to the file */
    if (fprintf(result, "%s", text2) < 0) {  //fwrite(text2, 1, sizeof(text2), result) != sizeof(text2)
        printf("Error: writing text2 to file in save_files function\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        if (fprintf(result, "%f ", vector[i]) < 0) {  //fwrite(vector, sizeof(float), n, result) != n
            printf("Error: writing eigenvector to file in save_files function\n");
            exit(1);
        }
        fprintf(result, "\n");
    }
    fclose(result);
    printf("Writing to files: Done!\n");
}

static inline float** read_matrix_file(char* __restrict file_name, const int n) {
    int re;  //used for return types
    FILE* input = fopen(file_name, "r");
    if (input == NULL) {
        printf("File can't be opened or does not exist!\n");
        exit(-1);
    }

    char buff[11 * n];  // n number of floats * 10 number of characterfor for each float + n spaces in between
    float** matrix = (float**)malloc(n * sizeof(float*));

    if (fgets(buff, 50, input) == NULL)  // to skip the first descriptive line
        printf("Error: fgets, error or end of file!\n");

    for (int i = 0; i < n; i++) {  // n number of lines or matrix rows
        matrix[i] = (float*)malloc(n * sizeof(float));

        for (int j = 0; j < n; j++) {
            re = fscanf(input, "%s", buff);
            matrix[i][j] = atof(buff);
        }
        if (re == 0) {
            printf("Error: fscanf, matching failure\n");
            exit(-1);
        }
    }
    printf("FILE: '%s' loaded!\n", file_name);
    printf("X****************************X\n");
    fclose(input);
    return matrix;
}

static inline float largest_absolute(float* __restrict vector, const int n) {
    float max;
    if (vector[0] < 0)
        max = -vector[0];
    else
        max = vector[0];

    for (int i = 1; i < n; i++) {
        if (vector[i] < 0 && -(vector[i]) > max) {
            max = -(vector[i]);
        }
        if (vector[i] > max)
            max = vector[i];
    }
    return max;
}

static inline float* normalize_vector(float* __restrict vector, const float max, const int n) {
    for (int i = 0; i < n; i++)
        vector[i] = vector[i] / max;
    return vector;
}

void power_method(float** __restrict matrix, const int iter, const float tolerance, const int n) {
    float max, check;
    max = check = 0;

    //float old_max = 0; // used in approach 1

    // used to save Product of matrix.vector multiplication P=AX
    float* product = (float*)malloc(n * sizeof(float));
    float* residual = (float*)malloc(n * sizeof(float));
    float* eigenvector = (float*)malloc(n * sizeof(float));

    /* create and initialize a nonzero n-by-1 vector */
    srand(time(NULL));
    for (int i = 0; i < n; ++i)
        eigenvector[i] = rand() % 20 - 10;

    max = largest_absolute(eigenvector, n);
    eigenvector = normalize_vector(eigenvector, max, n);

    int t = 0;
    for (t = 0; t < iter; t++) {
#pragma omp parallel shared(matrix, eigenvector, product)
        {
            float* product_v = (float*)malloc(n * sizeof(float));
            memset(product_v, 0, n * sizeof(float));
/* matrix2d_vector_ji, i for row, j for column */
#pragma omp for schedule(static)  //collapse(2)
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    product_v[j] += matrix[j][i] * eigenvector[i];
                }
            }
#pragma omp critical  //atomic
            {
                for (int k = 0; k < n; k++)
                    product[k] += product_v[k];
            }
            free(product_v);
        }
        //old_max = max;  /* comment out to use convergence method approach 1 */
        /* largest_absolute */
        if (product[0] < 0)
            max = -product[0];
        else
            max = product[0];
        for (int x = 1; x < n; x++) {
            if (product[x] < 0 && -(product[x]) > max)
                max = -(product[x]);

            if (product[x] > max)
                max = product[x];
        }

#if convergence_approach
        //Approach 1: convergence creteria
        check = max - old_max;
        if (check < 0) {
            check *= -1;
        }
        if (check >= 0 && check < tolerance) {
            printf("Tolerance reached: check:%f\n", check);
            break;
        }
#else
/* #pragma omp parallel
        {
            float* residual_v = (float*)malloc(n * sizeof(float)); */
//#pragma omp for schedule(static)
            for (int y = 0; y < n; y++)
                residual[y] = product[y] - max * eigenvector[y];

/* #pragma omp critical  //atomic
            {
                memcpy(residual,residual_v,n * sizeof(float));
            }
            free(residual_v);
        } */

        //Approach 2: convergence creteria
        check = largest_absolute(residual, n);
        if (check < tolerance) {
            printf("Tolerance reached: check:%f\n", check);
            break;
        }
#endif
        /* normalize_vector */
        for (int i = 0; i < n; i++) {
            eigenvector[i] = product[i] / max;
        }
        /* reseting 'product' to hold next matrix product multiplication */
        memset(product, 0, n * sizeof(float));
    }

if (t == iter) {
    printf("NO! convergence after %d iterations with %f as tolerance condition\n", iter, tolerance);
    printf("Value at last iteration: %f\n", max);
} else {
    printf("Eigenvalue:%f, iter:%d, check:%f\n", max, t, check);
    save_files(&max, eigenvector, matrix, n);
}

free(eigenvector);
free(product);
free(residual);

for (int i = 0; i < n; i++)
    free(matrix[i]);
free(matrix);
}

