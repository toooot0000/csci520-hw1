/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include "Vector3.h"
#include <cmath>

static inline void computeHookForceOnA(
    const Vector3 &a, 
    const Vector3 &b, 
    double rest, 
    double kHook,
    Vector3& accRet
) {
    auto l = a - b;
    accRet += -kHook * (l.length() - rest) * l.normalized();
}

static inline void computDampingForce(
    const Vector3& a,
    const Vector3& b,
    const Vector3& va,
    const Vector3& vb,
    double kd,
    Vector3& accRet
) {
    auto l = a - b;
    accRet += -kd * ((va - vb).dot(l)) / l.length() * l.normalized();
}

/// <summary>
/// Check i, j, k boundry. if inside the index boundry, ret the position and velocity.
/// </summary>
/// <param name="jello"></param>
/// <param name="i"></param>
/// <param name="j"></param>
/// <param name="k"></param>
/// <param name="retPos"></param>
/// <param name="retVel"></param>
/// <returns></returns>
static inline bool getPosAndVel(
    struct world* jello,
    int i, int j, int k,
    Vector3& retPos,
    Vector3& retVel
) {
#define guard(x) do{ if(x < 0 || x >= 8) return false; } while(0);
    guard(i);
    guard(j);
    guard(k);
    retPos.copyFrom(jello->p[i][j][k]);
    retVel.copyFrom(jello->v[i][j][k]);
    return true;
#undef guard
}




static void computeInnerForce(
    struct world* jello,
    int i, int j, int k,
    Vector3 &ret
) {
    Vector3 a, b, va, vb;
    getPosAndVel(jello, i, j, k, a, va);

#define compute_spring(i, j, k) do{ if (getPosAndVel(jello, i, j, k, b, vb)) {\
    computeHookForceOnA(a, b, rest, jello->kElastic, ret);\
    computDampingForce(a, b, va, vb, jello->dElastic, ret);\
}} while(0)
    // structual force
    auto rest = 1.0 / 7.0;
    compute_spring(i - 1, j, k);
    compute_spring(i + 1, j, k);
    compute_spring(i, j - 1, k);
    compute_spring(i, j + 1, k);
    compute_spring(i, j, k - 1);
    compute_spring(i, j, k + 1);

    // shear force
    rest = sqrt(2.0) / 7.0;
    compute_spring(i - 1, j - 1, k);
    compute_spring(i - 1, j + 1, k);
    compute_spring(i - 1, j, k - 1);
    compute_spring(i - 1, j, k + 1);
    compute_spring(i + 1, j - 1, k);
    compute_spring(i + 1, j + 1, k);
    compute_spring(i + 1, j, k - 1);
    compute_spring(i + 1, j, k + 1);
    compute_spring(i, j - 1, k - 1);
    compute_spring(i, j - 1, k + 1);
    compute_spring(i, j + 1, k - 1);
    compute_spring(i, j + 1, k + 1);

    // diagnal force
    rest = sqrt(3.0) / 7.0;
    compute_spring(i - 1, j - 1, k - 1);
    compute_spring(i - 1, j - 1, k + 1);
    compute_spring(i - 1, j + 1, k - 1);
    compute_spring(i - 1, j + 1, k + 1);
    compute_spring(i + 1, j - 1, k - 1);
    compute_spring(i + 1, j - 1, k + 1);
    compute_spring(i + 1, j + 1, k - 1);
    compute_spring(i + 1, j + 1, k + 1);


    // bend force
    rest = 2.0/7.0;
    compute_spring(i - 2, j, k);
    compute_spring(i + 2, j, k);
    compute_spring(i, j - 2, k);
    compute_spring(i, j + 2, k);
    compute_spring(i, j, k - 2);
    compute_spring(i, j, k + 2);

#undef compute_neighbor
}

/*
* Force field ret index
*           011---------111
*           /|          /|
*          / |         / |
*         /  |        /  |
*       001---------101  |
*        |  010----- | -110
*        |  /        |  /
*        | /         | /  
*        |/          |/
*       000---------100
* 
* params is the uniformed distance to (000)
* 
*/     

static inline int getIndex(size_t resolution, double pos, double &diff) {
    const double interval = 4.0 / (resolution - 1);
    auto ret = (size_t)floor((pos + 2.0) / interval);
    if(ret >= resolution - 1) return -1;
    if (ret < 0) return -1;
    diff = (pos + 2.0 - ret * interval)/interval;
    diff = min(1, diff);
    diff = max(0, diff);
    return ret;
}

static void computeForceField(
    struct world* jello,
    int ii, int jj, int kk,
    Vector3& accRet
) {
#define index(i, j, k) (jello->forceField)[(i)*jello->resolution*jello->resolution + (j)*jello->resolution + (k)]
    const Vector3& pos = { jello->p[ii][jj][kk] };
    int i, j, k;
    double xFac, yFac, zFac, curCoef;
    struct point curForce;
    Vector3 force;
    if ((i = getIndex(jello->resolution, pos.x(), xFac)) == -1) return;
    if ((j = getIndex(jello->resolution, pos.y(), yFac)) == -1) return;
    if ((k = getIndex(jello->resolution, pos.z(), zFac)) == -1) return;

    curForce = index(i, j, k);
    curCoef = (1 - xFac) * (1 - yFac) * (1 - zFac);
    curForce *= curCoef;
    accRet += curForce;

    curForce = index(i + 1, j, k);
    curCoef = (xFac) * (1 - yFac) * (1 - zFac);
    curForce *= curCoef;
    accRet += curForce;

    curForce = index(i, j + 1, k);
    curCoef = (1 - xFac) * (yFac) * (1 - zFac);
    curForce *= curCoef;
    accRet += curForce;

    curForce = index(i + 1, j + 1, k);
    curCoef = (xFac) * (yFac) * (1 - zFac);
    curForce *= curCoef;
    accRet += curForce;

    curForce = index(i, j, k + 1);
    curCoef = (1 - xFac) * (1 - yFac) * (zFac);
    curForce *= curCoef;
    accRet += curForce;

    curForce = index(i + 1, j, k + 1);
    curCoef = (xFac) * (1 - yFac) * (zFac);
    curForce *= curCoef;
    accRet += curForce;

    curForce = index(i, j + 1, k + 1);
    curCoef = (1 - xFac) * (yFac) * (zFac);
    curForce *= curCoef;
    accRet += curForce;

    curForce = index(i + 1, j + 1, k + 1);
    curCoef = (xFac) * (yFac) * (zFac);
    curForce *= curCoef;
    accRet += curForce;

#undef index()
}

static inline bool getPlaneCollision(
    double a, double b, double c, double d,
    const Vector3& p,
    Vector3& ret
) {
    auto btm = a * a + b * b + c * c;
    if (fabs(btm) < 0.0000001) return false;
    auto point = p.point();
    auto t = -(a * point.x + b * point.y + c * point.z + d) / (btm);
    if (t > 0) {
        ret = {point.x + a * t, point.y + b * t, point.z + c * t};
        return true;
    }
    else return false;
}

/// <summary>
/// Calculate if hitting bound, if hit, assign the bounding point to ret
/// </summary>
/// <param name="pos">tested position</param>
/// <param name="ret">hitting point on the surface of the bound</param>
/// <returns></returns>
static inline bool getBoundCollision(
    const Vector3& pos,
    Vector3& ret
) {
    if (pos.x() < -2) {
        ret = { -2, pos.y(), pos.z() };
        return true;
    }
    else if (pos.x() > 2) {
        ret = { 2, pos.y(), pos.z() };
        return true;
    } else if (pos.y() < -2) {
        ret = { pos.x(), -2, pos.z()};
        return true;
    }
    else if (pos.y() > 2) {
        ret = { pos.x(), 2, pos.z() };
        return true;
    }
    else if (pos.z() < -2) {
        ret = { pos.x(), pos.y(), -2};
        return true;
    }
    else if (pos.z() > 2) {
        ret = { pos.x(), pos.y(), 2};
        return true;
    }
    return false;
}

static void computeCollision(
    struct world* jello,
    int i, int j, int k,
    Vector3& accRet
) {
    Vector3 a, b, va;
    if (!getPosAndVel(jello, i, j, k, a, va)) return;
    if (!getBoundCollision(a, b)) return;
    auto &vb = Vector3::zero();
    computeHookForceOnA(a, b, 0, jello->kCollision, accRet);
    computDampingForce(a, b, va, vb, jello->dCollision, accRet);
}

static void computePlaneCollision(
    struct world* jello,
    int i, int j, int k,
    Vector3& accRet
) {
    Vector3 a, b, va;
    if (!getPosAndVel(jello, i, j, k, a, va)) return;
    if (!getPlaneCollision(jello->a, jello->b, jello->c, jello->d, a, b)) return;
    auto& vb = Vector3::zero();
    computeHookForceOnA(a, b, 0, jello->kCollision, accRet);
    computDampingForce(a, b, va, vb, jello->dCollision, accRet);
}


/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'acc'. */
void computeAcceleration(struct world * jello, struct point acc[8][8][8])
{
  /* for you to implement ... */
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                Vector3 cur;
                computeForceField(jello, i, j, k, cur);
                computeInnerForce(jello, i, j, k, cur);
                computeCollision(jello, i, j, k, cur);
                if (jello->incPlanePresent) computePlaneCollision(jello, i, j, k, cur);
                if (i == hitI && j == hitJ && k == hitK) {
                    cur += dragForce;
                }
                cur *= 1.0/jello->mass; 
                acc[i][j][k] = cur.point();
            }
        }
    }
    
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k], 0.5, buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k], 0.5, buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a); 

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
