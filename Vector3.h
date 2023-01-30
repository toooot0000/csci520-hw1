#pragma once
#include "jello.h"
#include <math.h>


class Vector3
{
private:
 	struct point *_p;
	static Vector3 _zero;

public:
	inline static const Vector3& zero() { 
		return _zero;
	}

#define X _p->x
#define Y _p->y
#define Z _p->z

	/// <summary>
	/// will create a point
	/// </summary>
	Vector3() : _p(new struct point()){};

	~Vector3() {
		if (!_p) return;
		delete _p;
	}
	/// <summary>
	/// Self managed
	/// </summary>
	/// <param name="x"></param>
	/// <param name="y"></param>
	/// <param name="z"></param>
	inline Vector3(double x, double y, double z): _p(new struct point()){
		X = x; Y = y; Z = z;
	}

	/// <summary>
	/// Copied from other, create new
	/// </summary>
	/// <param name="other"></param>
	inline Vector3(const struct point& other) : _p((new struct point())) {
		X = other.x; Y = other.y; Z = other.z;
	}

	/// <summary>
	/// use other's point structure
	/// </summary>
	inline Vector3(const Vector3& other): _p((new struct point())) {
		X = other.X; Y = other.Y; Z = other.Z;
	}


	/// <summary>
	/// Other will not be self managed.
	/// </summary>
	/// <param name="newVal"></param>
	/// <returns></returns>
	inline Vector3(Vector3&& other) noexcept : _p(other._p){
		other._p = nullptr;
	}

	/// <summary>
	/// No copy
	/// </summary>
	inline void operator=(const Vector3& other) {
		X = other.X; Y = other.Y; Z = other.Z;
	};


	/// <summary>
	/// Only Move
	/// </summary>
	inline void operator=(Vector3&& other) noexcept {
		_p = other._p;
		other._p = nullptr;
	}
	
	inline double x() const { return X; }
	inline double y() const { return Y; }
	inline double z() const { return Z; }

	inline void x(double newVal) { X = newVal; }
	inline void y(double newVal) { Y = newVal; }
	inline void z(double newVal) { Z = newVal; }


	inline const struct point& point() const { return *_p; }
	
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

