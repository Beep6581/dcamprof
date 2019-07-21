/*

  This code is an adaption of Jarno Elonen's TPS demo
  (http://elonen.iki.fi/code/tpsdemo/), which has the
  following copyright:

  Copyright (C) 2003, 2004 by Jarno Elonen

  Permission to use, copy, modify, distribute and sell this
  software and its documentation for any purpose is hereby
  granted without fee, provided that the above copyright
  notice appear in all copies and that both that copyright
  notice and this permission notice appear in supporting
  documentation. The authors make no representations about
  the suitability of this software for any purpose. It is
  provided "as is" without express or implied warranty.

*/
#include <stdio.h>
#include <math.h>

#include <matmath.h>
#include <elog.h>
#include <tps.h>

// Solve a linear equation system a*x=b using inplace LU decomposition.
//
// Stores x in 'b' and overwrites 'a' (with a pivotted LUD).
//
// Matrix 'b' may have any (>0) number of columns but
// must contain as many rows as 'a'.
//
// Possible return values:
//  0=success
//  1=singular matrix
//  2=a.rows != b.rows
static int
LU_Solve(matrix *a,
         matrix *b)
{
    // This routine is originally based on the public domain draft for JAMA,
    // Java matrix package available at http://math.nist.gov/javanumerics/jama/

//    typedef boost::numeric::ublas::matrix<T> Matrix;
//    typedef boost::numeric::ublas::matrix_row<Matrix> Matrix_Row;
//    typedef boost::numeric::ublas::matrix_column<Matrix> Matrix_Col;

    if (a->rows != b->rows) {
        return 2;
    }

    int m = a->rows, n = a->cols;
    int pivsign = 0;
    int piv[m];

    // PART 1: DECOMPOSITION
    //
    // For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
    // unit lower triangular matrix L, an n-by-n upper triangular matrix U,
    // and a permutation vector piv of length m so that A(piv,:) = L*U.
    // If m < n, then L is m-by-m and U is m-by-n.
    {
        // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
        for (int i = 0; i < m; ++i) {
            piv[i] = i;
        }
        pivsign = 1;

        // Outer loop.
        for (int j=0; j<n; ++j) {
            // Make a copy of the j-th column to localize references.
            // Apply previous transformations.
            for (int i = 0; i < m; ++i) {
                // This dot product is very expensive.
                // Optimize for SSE2?
                int kmax = (i<=j)?i:j;
                int k = 0;
                double sum = 0;
                while( kmax-- > 0 ) {
                    sum += (*rc(a, i, k)) * (*rc(a, k, j));
                    k++;
                }
                *rc(a, i, j) -= sum;
            }

            // Find pivot and exchange if necessary.
            //
            // Slightly optimized version of:
            //  for (int i = j+1; i < m; ++i)
            //    if ( fabs(LUcolj[i]) > fabs(LUcolj[p]) )
            //      p = i;
            int p = j;
            double coljp_abs = fabs(*rc(a,p,j));
            for (int i = j+1; i < m; i++) {
                if (fabs(*rc(a,i,j)) > coljp_abs) {
                    p = i;
                    coljp_abs = fabs(*rc(a,p,j));
                }
            }

            if (p != j) {
                swap_mrow(a, p, j);

                int tmp = piv[p];
                piv[p] = piv[j];
                piv[j] = tmp;
                pivsign = -pivsign;
            }

            // Compute multipliers.
            if (j < m && *rc(a,j,j) != 0.0) {
                for (int i = j+1; i < m; ++i) {
                    *rc(a,i,j) /= *rc(a,j,j);
                }
            }
        }
    }

    // PART 2: SOLVE

    // Check singluarity
    for (int j = 0; j < n; ++j) {
        if (*rc(a,j,j) == 0.0) {
            return 1;
        }
    }

    // Reorder b according to pivotting
    for (int i=0; i<m; ++i) {
        if ( piv[i] != i ) {
            swap_mrow(b, piv[i], i);
            for ( int j=i; j<m; ++j ) {
                if ( piv[j] == i ) {
                    piv[j] = piv[i];
                    break;
                }
            }
        }
    }

    // Solve L*Y = B(piv,:)
    for (int k=0; k<n; ++k) {
        for (int i = k+1; i < n; i++) {
            const double aik = *rc(a,i,k);
            for (int col = 0; col < b->cols; col++) {
                *rc(b, i, col) -= *rc(b, k, col) * aik;
            }
        }
    }

    // Solve U*X = Y;
    for (int k=n-1; k>=0; k--) {
        const double akk = 1.0 / *rc(a,k,k);
        for (int col = 0; col < b->cols; col++) {
            *rc(b, k, col) *= akk;
        }
        for (int i = 0; i < k; i++) {
            const double aik = *rc(a,i,k);
            for (int col = 0; col < b->cols; col++) {
                *rc(b, i, col) -= *rc(b, k, col) * aik;
            }
        }
    }

    return 0;
}

static inline double
tps_base_func(double r)
{
    if ( r == 0.0 )
        return 0.0;
    else
        return r*r * log(r);
}

static inline double
difflen(const v3 c1,
        const v3 c2)
{
    v3 c = {{ c1.v[0]-c2.v[0], c1.v[1]-c2.v[1], c1.v[2]-c2.v[2] }};

    return sqrt(c.v[0]*c.v[0] + c.v[1]*c.v[1] + c.v[2]*c.v[2]);
}

struct tps_t_ {
    matrix *mtx_v;
    matrix *mtx_orig_k;
    v3 *control_points;
    size_t control_points_len;
};

void
tps_delete(tps_t *tps)
{
    if (tps == NULL) {
        return;
    }
    free(tps->mtx_v);
    free(tps->mtx_orig_k);
    free(tps->control_points);
    free(tps);
}

/*
 *  Calculate Thin Plate Spline (TPS) weights from
 *  control points and build a new height grid by
 *  interpolating with them.
 */
tps_t *
tps_new(const v3 control_points[],
        size_t control_points_len,
        const double regularization)
{
    // You We need at least 3 points to define a plane
    if (control_points_len < 3 ) {
        elog("TPS: too few control points! Aborting.\n");
        abort();
    }

    const unsigned p = control_points_len;

    // Allocate the matrix and vector
    matrix *mtx_l = alloc_matrix(p+3, p+3);
    matrix *mtx_v = alloc_matrix(p+3, 1);
    matrix *mtx_orig_k = alloc_matrix(p, p);

    // Fill K (p x p, upper left of L) and calculate
    // mean edge length from control points
    //
    // K is symmetrical so we really have to
    // calculate only about half of the coefficients.
    double a = 0.0;
    for ( unsigned i=0; i<p; ++i ) {
        for ( unsigned j=i+1; j<p; ++j ) {
            v3 pt_i = control_points[i];
            v3 pt_j = control_points[j];
            pt_i.v[1] = pt_j.v[1] = 0;
            double elen = difflen(pt_i, pt_j);
            *rc(mtx_l,i,j) = *rc(mtx_l,j,i) = *rc(mtx_orig_k,i,j) = *rc(mtx_orig_k,j,i) = tps_base_func(elen);
            a += elen * 2; // same for upper & lower tri
        }
    }
    a /= (double)(p*p);

    // Fill the rest of L
    for ( unsigned i=0; i<p; ++i ) {
        // diagonal: reqularization parameters (lambda * a^2)
        *rc(mtx_l, i,i) = *rc(mtx_orig_k,i,i) = regularization * (a*a);

        // P (p x 3, upper right)
        *rc(mtx_l, i, p+0) = 1.0;
        *rc(mtx_l, i, p+1) = control_points[i].v[0];
        *rc(mtx_l, i, p+2) = control_points[i].v[2];

        // P transposed (3 x p, bottom left)
        *rc(mtx_l, p+0, i) = 1.0;
        *rc(mtx_l, p+1, i) = control_points[i].v[0];
        *rc(mtx_l, p+2, i) = control_points[i].v[2];
    }
    // O (3 x 3, lower right)
    for ( unsigned i=p; i<p+3; ++i ) {
        for ( unsigned j=p; j<p+3; ++j ) {
            *rc(mtx_l, i,j) = 0.0;
        }
    }

    // Fill the right hand vector V
    for ( unsigned i=0; i<p; ++i ) {
        *rc(mtx_v, i,0) = control_points[i].v[1];
    }
    *rc(mtx_v, p+0, 0) = *rc(mtx_v, p+1, 0) = *rc(mtx_v, p+2, 0) = 0.0;

    // Solve the linear system "inplace"
    if (LU_Solve(mtx_l, mtx_v) != 0) {
        // Singular matrix! No solution
        free(mtx_l);
        free(mtx_v);
        free(mtx_orig_k);
        return NULL;
    }

    // Interpolate grid heights
    /*
#define GRID_W 9
#define GRID_H 9
    float grid[GRID_W][GRID_H];

    for ( int x=-GRID_W/2; x<GRID_W/2; ++x ) {
        for ( int z=-GRID_H/2; z<GRID_H/2; ++z ) {
            double h = *rc(mtx_v, p+0, 0) + *rc(mtx_v, p+1, 0)*x + *rc(mtx_v, p+2, 0)*z;
            v3 pt_i, pt_cur = {{x,0,z}};
            for ( unsigned i=0; i<p; ++i ) {
                pt_i = control_points[i];
                pt_i.v[1] = 0;
                h += *rc(mtx_v, i,0) * tps_base_func(difflen(pt_i, pt_cur));
            }
            grid[x+GRID_W/2][z+GRID_H/2] = h;
            elog("%.2f ", h);
        }
        elog("\n");
    }
    elog("\n");

    printf("mtx_orig_k:\n");
    for (int i = 0; i < mtx_orig_k->rows; i++) {
        for (int j = 0; j < mtx_orig_k->cols; j++) {
            printf("%.2f ", *rc(mtx_orig_k, i, j));
        }
        printf("\n");
    }
    printf("\n");

    printf("mtx_v:\n");
    for (int i = 0; i < mtx_v->rows; i++) {
        for (int j = 0; j < mtx_v->cols; j++) {
            printf("%.2f ", *rc(mtx_v, i, j));
        }
        printf("\n");
    }
    */

    free(mtx_l);
    tps_t *tps = calloc(1, sizeof(*tps));
    tps->mtx_v = mtx_v;
    tps->mtx_orig_k = mtx_orig_k;
    tps->control_points = malloc(p*sizeof(control_points[0]));
    memcpy(tps->control_points, control_points, p*sizeof(control_points[0]));
    tps->control_points_len = p;

    return tps;
}

double
tps_bending_energy(tps_t *tps)
{
    // Calc bending energy
    const unsigned p = tps->control_points_len;
    matrix *w = alloc_matrix(p, 1);
    matrix *w_trans = alloc_matrix(1, p);
    for (unsigned i = 0; i < p; i++) {
        *rc(w,i,0) = *rc(tps->mtx_v,i,0);
        *rc(w_trans,0,i) = *rc(tps->mtx_v,i,0);
    }
    matrix *prod = matrix_mul(w_trans, tps->mtx_orig_k);
    free(w_trans);
    matrix *be = matrix_mul(prod, w);
    free(prod);
    free(w);
    double bending_energy = *rc(be,0,0);
    free(be);
    return bending_energy;
}

double
tps_interpolate(tps_t *tps,
                double x,
                double z)
{
    const unsigned p = tps->control_points_len;
    double h = *rc(tps->mtx_v, p+0, 0) + *rc(tps->mtx_v, p+1, 0)*x + *rc(tps->mtx_v, p+2, 0)*z;
    v3 pt_i, pt_cur = {{x,0,z}};
    for (unsigned i = 0; i < p; i++) {
        pt_i = tps->control_points[i];
        pt_i.v[1] = 0;
        h += *rc(tps->mtx_v, i,0) * tps_base_func(difflen(pt_i, pt_cur));
    }
    return h;
}
