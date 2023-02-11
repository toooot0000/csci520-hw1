#pragma once
#include "Vector3.h"


//class Mat4 {
//private:
//	union {
//		double a[16];
//		double b[4][4];
//	}m;
//public:
//	Mat4() {
//		for (int i = 0; i < 16; i++) m.a[i] = 0;
//	}
//
//	inline double operator[](size_t i) const { return m.a[i]; }
//	inline double& operator[](size_t i) { return m.a[i]; }
//
//	//Mat4 inversed() const;
//
//};
extern void inverse(const double a[4][4], double ret[4][4]);

extern void multiply(const double lhs[4][4], const double rhs[4][4], double ret[4][4]);

extern void transform(const double a[4][4], const Vector3& point, Vector3& ret);

inline void print(const double m[4][4]) {
	for (int i = 0; i < 4; i++)
		printf("%.3f %.3f %.3f %.3f\n", m[0][i], m[1][i], m[2][i], m[3][i]);
}