package mathcomp.oletsky.mathhelper;

import java.util.Arrays;

public class VectMatr {

    public static double[] getRowSums(double[][] matr) {
        int nRows = matr.length;
        int nCols = matr[0].length;
        double[] res = new double[nRows];
        for (int i = 0; i < nRows; i++) {
            res[i] = VectMatr.calculateSumOfComponents(matr[i]);
        }
        return res;
    }

    public static double[] getColumnSums(double[][] matr) {
        int nRows = matr.length;
        int nCols = matr[0].length;
        double[] res = new double[nCols];
        for (int j = 0; j < nCols; j++) {
            double s = 0;
            for (int i = 0; i < nRows; i++) {
                s += matr[i][j];
            }
            res[j] = s;
        }
        return res;
    }

    public static double calculateAverage(double[] v) {
        return VectMatr.calculateSumOfComponents(v) / v.length;
    }

    public static double calculateSumOfComponents(double[] v) {
        double s = 0.;
        for (double d : v) s += d;
        return s;
    }

    public static double calculateProdOfComponents(double[] v) {
        double s = 1.;
        for (double d : v) s *= d;
        return s;
    }

    public static double calculateEvklVectorNorm(double[] v) {
        double s = 0;
        for (double d : v) s += (d * d);
        return Math.sqrt(s);
    }

    public static double[] calculateAverageVector(double[] a, double[] b) {
        if (a.length != b.length) throw new RuntimeException("Vectors are of different length!");
        int n = a.length;
        double[] res = new double[n];
        for (int i = 0; i < n; i++) {
            res[i] = (a[i] + b[i]) / 2.;
        }
        return res;
    }


    public static void defaultOutputVector(double[] vect) {
        for (int j = 0; j < vect.length; j++) {
            System.out.printf("%11.4f", vect[j]);
        }
        System.out.print("\n");
    }

    public static void defaultOutputMatrix(double[][] matr) {
        for (int j = 0; j < matr.length; j++) {
            defaultOutputVector(matr[j]);
        }

    }

    public static double[] normalizeVectorBySum(double[] v) {
        int n = v.length;
        double[] res = new double[n];
        double s = VectMatr.calculateSumOfComponents(v);
        for (int i = 0; i < n; i++) {
            res[i] = v[i] / s;
        }

        return res;
    }

    public static double[] normalizeVectorByEvkl(double[] v) {
        int n = v.length;
        double[] res = new double[n];
        double s = VectMatr.calculateEvklVectorNorm(v);
        for (int i = 0; i < n; i++) {
            res[i] = v[i] / s;
        }

        return res;
    }

    public static double calculateScalarProduct(double[] a, double[] b) {
        int n = a.length;

        if (n != b.length) throw new RuntimeException("Vectors have different length!");
        double d = 0.;
        for (int i = 0; i < n; i++) {
            d += a[i] * b[i];
        }

        return d;
    }

    public static double[] leftMultiply(double[] v, double[][] m) {
        int n = v.length;
        int nMatr = m[0].length;
        int rows = m.length;
        if (n != rows) throw new RuntimeException("Wrong dimensions!");
        double res[] = new double[nMatr];
        for (int j = 0; j < nMatr; j++) {
            double s = 0.;
            for (int k = 0; k < rows; k++) {
                s += v[k] * m[k][j];
            }
            res[j] = s;
        }

        return res;
    }

    public static double[] rightMultiply(double[][] m, double[] v) {
        int n = v.length;
        int rows = m.length;
        double res[] = new double[rows];
        for (int i = 0; i < rows; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) {
                s += m[i][j] * v[j];
            }
            res[i] = s;
        }

        return res;
    }


    public static double[] getLeftEigenVector(double[][] m) {

        final double EPS = 1.E-15;
        int n = m.length;
        double[] prev = new double[n];
        double[] next = new double[n];

        for (int i = 0; i < n; i++) {
            prev[i] = 1. / n;
            next[i] = 1. / n;
        }

        do {
            prev = next;
            next = leftMultiply(prev, m);
        } while (calculateDistance(next, prev) > EPS);
        return next;
    }

    public static double[] getMainEigenVector(double[][] m) {

        final double EPS = 1.E-15;
        int n = m.length;
        double[] prev = new double[n];
        double[] next = new double[n];

        for (int i = 0; i < n; i++) {
            prev[i] = 1.;
            next[i] = 1.;
        }

        do {
            prev = next;
            next = rightMultiply(m, prev);
            next = VectMatr.normalizeVectorByEvkl(next);
        } while (calculateDistance(next, prev) > EPS);
        return next;
    }

    public static double calculateDistance(double[] a, double[] b) {
        int alen = a.length;
        int blen = b.length;
        if (alen != blen)
            throw new RuntimeException("Different vector dimensions");
        double s = 0.;
        for (int i = 0; i < alen; i++) {
            s += ((a[i] - b[i]) * (a[i] - b[i]));
        }
        return Math.sqrt(s);
    }

    public static double[][] transposeSquareMatrix(double[][] matr) {
        int n = matr.length;
        double[][] transposed = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                transposed[i][j] = matr[j][i];
            }
        return transposed;
    }

    public static double[][] transposeMatrix(double[][] matr) {
        int m = matr.length;
        int n = matr[0].length;
        double[][] transposed = new double[n][m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                transposed[i][j] = matr[j][i];
            }
        return transposed;
    }

    public static double[][] multiplyMatrices(double[][] a, double[][] b) {
        int m = a.length;
        int r = a[0].length;
        int r1 = b.length;
        int n = b[0].length;
        if (r != r1) throw new RuntimeException("Invalid dimensions");
        double[][] res = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                double s = 0;
                for (int k = 0; k < r; k++) {
                    s += a[i][k] * b[k][j];
                }
                res[i][j] = s;
            }
        return res;
    }

    public static double[] getMatrixRow(double[][] m, int i) {
        int n = m[i].length;
        double[] res = Arrays.copyOf(m[i], n);
        return res;
    }

    public static double[] getMatrixColumn(double[][] m, int j) {
        double[][] transp = VectMatr.transposeMatrix(m);
        return VectMatr.getMatrixRow(transp, j);
    }

    public static double calculateGeomAverage(double[] v) {
        double pw = 1. / v.length;
        double prod = VectMatr.calculateProdOfComponents(v);
        return Math.pow(prod, pw);
    }

    public static double[][] normalizeMatrixByRows(double[][] matr) {
        int m = matr.length;
        int n = matr[0].length;
        double[][] res = new double[m][n];
        for (int i = 0; i < m; i++) {
            res[i] = VectMatr.normalizeVectorBySum(matr[i]);
        }
        return res;
    }

    public static double[] calculatePercentage(double[] v) {
        int n = v.length;
        double[] res = new double[n];
        double s = VectMatr.calculateSumOfComponents(v);
        for (int i = 0; i < n; i++) {
            res[i] = v[i] / s;
        }
        return res;
    }

    public static double getMinValue(double[] v) {
        int n = v.length;
        double min = v[0];
        if (n == 1) return min;
        for (int i = 1; i < n; i++) {
            if (v[i] < min) min = v[i];
        }
        return min;
    }

    public static double getMaxValue(double[] v) {
        int n = v.length;
        double max = v[0];
        if (n == 1) return max;
        for (int i = 1; i < n; i++) {
            if (v[i] > max) max = v[i];
        }
        return max;
    }

    public static double getMaxValueInMatrix(double[][] matr) {
        int m = matr.length;
        int n = matr[0].length;
        double max = matr[0][0];
        if ((n == 1) && (m == 1)) return max;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (matr[i][j] > max) max = matr[i][j];
            }
        }
        return max;
    }

    public static double[][] invertByValue(double[][] matr,
                                           double value) {
        int m = matr.length;
        int n = matr[0].length;
        double[][] res = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res[i][j] = value - matr[i][j];
            }

        }
        return res;
    }

    public static int getMaxIndex(double[] v) {
        int n = v.length;
        double max = v[0];
        int maxIndex = 0;
        if (n == 1) return maxIndex;
        for (int i = 1; i < n; i++) {
            if (v[i] > max) {
                max = v[i];
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    public static double[] normalizeByMax(double[] v) {
        int n = v.length;
        double[] res = new double[n];
        double max = VectMatr.getMaxValue(v);
        for (int i = 0; i < n; i++) {
            res[i] = v[i] / max;
        }
        return res;
    }

    public static double[] normalizeMaxMin(double[] v) {
        int n = v.length;
        double[] res = new double[n];
        double max = VectMatr.getMaxValue(v);
        double min = VectMatr.getMinValue(v);
        double diff = max - min;
        for (int i = 0; i < n; i++) {
            res[i] = (v[i] - min) / diff;

        }
        return res;
    }

    public static double[] selectColumn(double[][] matr, int j) {
        int m = matr.length;
        int n = matr[0].length;
        if (j >= n) throw new RuntimeException("Index out of range");
        double[] res = new double[m];
        for (int i = 0; i < m; i++) {
            res[i] = matr[i][j];
        }
        return res;
    }

    public static double[] getUnityVector(int n) {
        double[] arr = new double[n];
        for (int i = 0; i < n; i++) {
            arr[i] = 1.;
        }
        return arr;
    }

    public static double getMinElement(double[][] matr) {
        double min = matr[0][0];
        for (int i = 0; i < matr.length; i++) {
            for (int j = 0; j < matr[0].length; j++) {
                if (matr[i][j] < min) {
                    min = matr[i][j];
                }

            }


        }
        return min;
    }

    public static void shiftMatrixValues(double[][] matr, double beta) {
        for (int i = 0; i < matr.length; i++) {
            for (int j = 0; j < matr[0].length; j++) {
                matr[i][j] += beta;

            }

        }


    }

    public static double[][] normalizeMatrixColumns(
            double[][] matr) {
        int m = matr.length;
        int n = matr[0].length;
        double[][] res = new double[m][n];
        for (int j = 0; j < n; j++) {
            double[] col = VectMatr.getMatrixColumn(matr, j);
            double[] normCol = VectMatr.normalizeVectorBySum(col);
            for (int i = 0; i < m; i++) {
                res[i][j] = normCol[i];
            }
        }
        return res;

    }

    public static double[][] multMatrixScalar(
            double[][] matr,
            double scal
    ) {
        int m = matr.length;
        int n = matr[0].length;
        double[][] res = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                res[i][j] = matr[i][j] * scal;

            }

        }
        return res;
    }

    public static int findIndexByValue(double[] vect,
                                       double value) {
        final double EPS=1.E-7;
        int n = vect.length;
        int index = -1;
        for (int i=0; i<n; i++) {
            if (Math.abs(vect[i]-value)<EPS) {
                index=i;
                break;
            }
        }
        return index;
    }

    public static double[] calculateConvexHull(
            double[] alpha,
            double[]... vects
    ){
        int m=vects[0].length;
        int kolOpor= vects.length;
        double[] conv = new double[m];

        for (int i = 0; i <m ; i++) {
            double sConv=0.;
            for (int k = 0; k <kolOpor ; k++) {
                sConv+=alpha[k]*vects[k][i];
            }
            conv[i]=sConv;
        }
        return conv;
    }

}



