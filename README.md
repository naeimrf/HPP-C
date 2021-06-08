### High Performance Programming codes in C programming

Eigenvalues, also known as characteristic roots, are numbers that can characterize a matrix and
reveal important behavior of it in a simplified way.
They are associated to matrix eigenvectors. In linear algebra, in a linear system of equation, a linear
transformation can be expressed by a nonzero vector when the direction of that vector remains
unchanged while that transformation is applied to it. This vector is called eigenvector or
characteristic vector.

The determination of eigenvectors and their associated eigenvalues are important in physics and
engineering and it is used in a variety of areas such as stability analysis, physics of rotating bodies,
small oscillations of vibrating systems electric network analysis, quantum mechanics, the
stability analysis of discretization methods for systems of ordinary differential equations, etc.

“For any square matrix M of size m×m, eigenvalues are called lambda λ and associated with an
eigenvector v if:

M.v=λv⟺(M−λI).v=0

when I is the identity matrix (of size m).

Practically, the eigenvalues of M are the roots of its characteristic polynomial P".
In general a square matrix of size m, has m number of eigenvalues but some or all of them may be
complex numbers. A good visualization of these concepts can be found here.

Up to some size of matrix, it is easily possible to use the equation above and find eigenvalues and
eigenvectors. But it gets harder when the size of matrix starts to grow. In engineering applications,
it is usually enough to know the smallest or the biggest eigenvalue. There are a number of
numerical algorithms that can help finding dominant eigenvalue of a square matrix which among
such power method is one of them.

Like Jacobi and Gauss-Seidel methods, the power method for approximating eigenvalues is
iterative. This method is very good at approximating the dominant eigenvalues and especially
for matrices with many zero elements for example sparse matrices.

