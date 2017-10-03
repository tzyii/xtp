#include <votca/xtp/gpu/GPUCTypes.cuh>

namespace votca { namespace xtp { namespace gpu {

__device__ double3 operator+(const double3& a, const double3& b){
    double3 c;
    c.x = a.x+b.x;
    c.y = a.y+b.y;
    c.z = c.z+b.z;
    return c;
}

__device__ double3 operator*(const double& s, const double3& a){
    double3 b;
    b.x = s*a.x;
    b.y = s*a.y;
    b.z = s*a.z;
    return b;
}

__device__ double3 operator*(const double3& a, const double s){
    double3 b;
    b.x = s*a.x;
    b.y = s*a.y;
    b.z = s*a.z;
    return b;
}

__device__ double3 operator/(const double3& a, const double s){
    return (1/s) * a; 
}

__device__ double3 operator-(const double3& a, const double3& b){
    double3 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c; 
}

__device__ void operator+=(double3& a, const double3& b){
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}

__device__ void operator-=(double3& a, const double3& b){
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z; 
}

__device__ double3 R3NIndex(const R3Nptrs& a, const size_t& i){
    return make_double3(a.x[i], a.y[i], a.z[i]);
}

__device__ double dot(const double3& a, const double3& b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

}}}
