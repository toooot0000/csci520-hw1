#include "Mat.h"
#include <math.h>
static double det(const double m[4][4]) {
    double r0 = m[2][3] * m[3][2] * m[0][1] * m[1][0] - m[2][2] * m[3][3] * m[0][1] * m[1][0] - m[2][3] * m[3][1] * m[0][2] * m[1][0] + m[2][1] * m[3][3] * m[0][2] * m[1][0],
        r1 = m[2][2] * m[3][1] * m[0][3] * m[1][0] - m[2][1] * m[3][2] * m[0][3] * m[1][0] - m[0][0] * m[2][3] * m[3][2] * m[1][1] + m[0][0] * m[2][2] * m[3][3] * m[1][1],
        r2 = m[2][3] * m[3][0] * m[0][2] * m[1][1] - m[2][2] * m[3][0] * m[0][3] * m[1][1] + m[0][0] * m[2][3] * m[3][1] * m[1][2] - m[0][0] * m[2][1] * m[3][3] * m[1][2],
        r3 = m[2][3] * m[3][0] * m[0][1] * m[1][2] + m[2][1] * m[3][0] * m[0][3] * m[1][2] - m[0][0] * m[2][2] * m[3][1] * m[1][3] + m[0][0] * m[2][1] * m[3][2] * m[1][3],
        r4 = m[2][2] * m[3][0] * m[0][1] * m[1][3] - m[2][1] * m[3][0] * m[0][2] * m[1][3] - m[3][3] * m[0][2] * m[1][1] * m[2][0] + m[3][2] * m[0][3] * m[1][1] * m[2][0],
        r5 = m[3][3] * m[0][1] * m[1][2] * m[2][0] - m[3][1] * m[0][3] * m[1][2] * m[2][0] - m[3][2] * m[0][1] * m[1][3] * m[2][0] + m[3][1] * m[0][2] * m[1][3] * m[2][0];
    return r0 + r1 + r2 - r3 + r4 + r5;
}

static double det(const double m[3][3]) {
    return m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1]
        - m[0][2] * m[1][1] * m[2][0] - m[0][1] * m[1][0] * m[2][2] - m[0][0] * m[1][2] * m[2][1];
}

static void getMinor(const double m[4][4], int ii, int jj, double ret[3][3]) {
    int iHit = 0, jHit = 0;
    for (int i = 0; i < 4; i++) {
        if (i == ii) {
            iHit = 1;
            continue;
        }
        jHit = 0;
        for (int j = 0; j < 4; j++) {
            if (j == jj) {
                jHit = 1;
                continue;
            }
            ret[i - iHit][j - jHit] = m[i][j];
        }
    }
}


void inverse(const double a[4][4], double ret[4][4]) {
    double temp[3][3];
    double detA = det(a), detTemp;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            getMinor(a, j, i, temp);
            detTemp = det(temp);
            ret[i][j] = (1 - (((i + j) & 1) << 1)) * detTemp /detA;
        }
    }
}

void multiply(const double lhs[4][4], const double rhs[4][4], double ret[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ret[i][j] = lhs[i][0] * rhs[0][j] + lhs[i][1] * rhs[1][j] + lhs[i][2] * rhs[2][j] + lhs[i][3] * rhs[3][j];
        }
    }

}

void transform(const double a[4][4], const Vector3& point, Vector3& ret) {
    double w = 0;
    for (int i = 0; i < 4; i++) {
        if (i == 3) {
            //w += point.x() * a[i][0] + point.y() * a[i][1] + point.z() * a[i][2] + a[i][3];
            w += point.x() * a[0][i] + point.y() * a[1][i] + point.z() * a[2][i] + a[3][i];
        }
        else {
            //ret.set(i, point.x() * a[i][0] + point.y() * a[i][1] + point.z() * a[i][2] + a[i][3]);
            ret.set(i, point.x() * a[0][i] + point.y() * a[1][i] + point.z() * a[2][i] + a[3][i]);
        }
    }
    if (fabs(w) < 0.00000000001) return;
    ret.x(ret.x() / w);
    ret.y(ret.y() / w);
    ret.z(ret.z() / w);
}