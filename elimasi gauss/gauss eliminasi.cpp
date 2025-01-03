#include <iostream>
#include <iomanip>
#include <cmath>

#define MAX 10
#define MAX_ITER 1000
#define EPSILON 1e-9

using namespace std;

double absolute(double num) {
    return num < 0 ? -num : num;
}

void displayMatrix(double mat[MAX][MAX + 1], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            cout << setw(10) << fixed << setprecision(4) << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int gaussElimination(double mat[MAX][MAX + 1], double sol[MAX], int n) {
    int iter = 0;

    for (int i = 0; i < n; i++) {
        cout << "\nLangkah " << i + 1 << ":\n";

        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (absolute(mat[k][i]) > absolute(mat[maxRow][i])) {
                maxRow = k;
            }
        }

        if (absolute(mat[maxRow][i]) < EPSILON) {
            cout << "Error: Matrix singular atau tidak memiliki solusi unik.\n";
            return -1;
        }

        if (maxRow != i) {
            for (int k = i; k <= n; k++) {
                swap(mat[i][k], mat[maxRow][k]);
            }
        }

        for (int k = i + 1; k < n; k++) {
            double factor = mat[k][i] / mat[i][i];
            for (int j = i; j <= n; j++) {
                mat[k][j] -= factor * mat[i][j];
            }
            iter++;
            if (iter > MAX_ITER) {
                cout << "Error: Iterasi melebihi batas maksimal.\n";
                return -1;
            }
        }

        displayMatrix(mat, n);
    }

    cout << "\nMatriks setelah eliminasi gauss:\n";
    displayMatrix(mat, n);

    for (int i = n - 1; i >= 0; i--) {
        sol[i] = mat[i][n] / mat[i][i];
        for (int j = i - 1; j >= 0; j--) {
            mat[j][n] -= mat[j][i] * sol[i];
            iter++;
            if (iter > MAX_ITER) {
                cout << "Error: Iterasi melebihi batas maksimal.\n";
                return -1;
            }
        }
    }

    return 0;
}

int main() {
    int n;
    double mat[MAX][MAX + 1], sol[MAX];

    cout << "Masukkan ukuran matriks (n x n): ";
    cin >> n;

    if (n <= 0 || n > MAX) {
        cout << "Error: Ukuran matriks tidak valid.\n";
        return -1;
    }

    cout << "Masukkan elemen matriks augmented:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            cout << "mat[" << i << "][" << j << "]: ";
            cin >> mat[i][j];
        }
    }

    cout << "\nMatriks sebelum eliminasi gauss:\n";
    displayMatrix(mat, n);

    if (gaussElimination(mat, sol, n) == 0) {
        cout << "\nSolusi:\n";
        for (int i = 0; i < n; i++) {
            cout << "x[" << i << "] = " << setw(10) << fixed << setprecision(4) << sol[i] << endl;
        }
    } else {
        cout << "Proses eliminasi gagal.\n";
    }

    return 0;
}
