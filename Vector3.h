#pragma once
#include "jello.h"
#include <math.h>


class Vector3
{
private:
 	struct point &_p;
	bool selfManaged = false;
	static Vector3 _zero;

public:
	inline static const Vector3& zero() { 
		return _zero;
	}

#define X _p.x
#define Y _p.y
#define Z _p.z

	Vector3() : _p(*(new struct point())), selfManaged(true){};

	~Vector3() {
		if (!selfManaged) return;
		delete &_p;
	}

	inline Vector3(double x, double y, double z): _p(*(new struct point())), selfManaged(true){
		X = x; Y = y; Z = z;
	}

	inline Vector3(struct point& p):_p(p){}

	inline Vector3(const struct point& other) : _p(*(new struct point())), selfManaged(true) {
		X = other.x; Y = other.y; Z = other.z;
	}

	inline Vector3(const Vector3& other):_p(other._p){
		X = other.X;
		Y = other.Y;
		Z = other.Z;
	}

	inline Vector3(Vector3&& other) noexcept : _p(other._p) {
		other.selfManaged = false;
	}

	inline void operator=(const Vector3& other) = delete;

	inline void operator=(Vector3&& other) noexcept {
		_p = other._p;
		other.selfManaged = false;
	}
	
	inline double x() const { return X; }
	inline double y() const { return Y; }
	inline double z() const { return Z; }

	inline const struct point& point() const { return _p; }
	
	inline Vector3 operator*(double other) const{
		return { X * other, Y * other, Z * other };
	}

	inline double operator*(const Vector3& other) const{
		return X * other.X + Y * other.Y + Z * other.Z ;
	}

	inline void copyFrom(const Vector3& other) {
		X = other.X;
		Y = other.Y;
		Z = other.Z;
	}

	inline void copyFrom(const struct point& other) {
		X = other.x;
		Y = other.y;
		Z = other.z;
	}

	inline double lengthSquare() const { return X * X + Y * Y + Z * Z; }
	
	inline double length() const { return sqrt(lengthSquare()); }

	inline void normalize(){ 
		auto l = length();
		X /= l;
		Y /= l;
		Z /= l;
	}

	inline Vector3 normalized() const {
		auto l = length();
		return { X / l, Y / l, Z / l };
	}

	inline Vector3 operator-() const {
		return { -X, -Y, -Z };
	}

	inline void operator+=(const Vector3& other) {
		X += other.X;
		Y += other.Y;
		Z += other.Z;
	}
	
	inline void operator-=(const Vector3& other) {
		operator+=(-other);
	}
	
	inline static friend Vector3 operator+(const Vector3& lhs, const Vector3& rhs) {
		return { lhs.X + rhs.X, lhs.Y + rhs.Y, lhs.Z + rhs.Z };
	}

	inline static friend Vector3 operator-(const Vector3& lhs, const Vector3& rhs) {
		return { lhs.X - rhs.X, lhs.Y - rhs.Y, lhs.Z - rhs.Z };
	}
	
	inline static friend Vector3 operator*(float lhs, const Vector3& rhs) {
		return rhs * lhs;
	}

	inline void reset() { X = 0; Y = 0; Z = 0; }

	inline void operator*=(float other) {
		X *= other;
		Y *= other;
		Z *= other;
	}

	inline void clampAll(double low, double high) {
#define clamp(x, lo, hi) ((x)<(lo)) ? (lo) : (((x) > (hi)) ? (hi) : (x) )
		X = clamp(X, low, high);
		Y = clamp(Y, low, high);
		Z = clamp(Z, low, high);
#undef clamp
	}

#undef X
#undef Y
#undef Z
};

