/*******************************************************************
*   KAero_CUDA.cpp
*   KPS
*
*	Author: Kareem Omar
*	kareem.omar@uah.edu
*	https://github.com/komrad36
*
*	Last updated Feb 12, 2016
*   This application is entirely my own work.
*******************************************************************/
//
// KAero uses CUDA or CPU and operates on a satellite or body defined as a series of polygons
// in 3-D space. When supplied with air density and velocity (in the Body frame),
// it approximates the drag force and torque on the body by simulating
// collisions and accumulating impulses and angular impulses per unit time.
//
// Note that net force and torque are returned in the Body frame.
//
// This module is intended for use by a 6 DoF orbital/attitude propagator which
// calls the aer() function from its integrators to obtain forces and/or torques,
// such as KPS.
//

#include "KAero.h"

#define PAD						(0.000001)

// placeholder x-value to signal no collision occurred
#define NONE					(-9999999.0)

// threads in one dimension
//
//
#define THREADS					(16)

// wrap around a CUDA function call to check and report any errors
bool checkForCUDAError(cudaError_t cuda_status) {
	if (cuda_status != cudaSuccess) {
		std::cerr << "ERROR: CUDA reports " << cudaGetErrorName(cuda_status) << ". Aborting." << std::endl
			<< "More info: " << cudaGetErrorString(cuda_status) << std::endl;
		return true;
	}
	return false;
}

// calculate the number of compute cores in a CUDA device
// based on its major and minor version numbers
//
// used in finding the fastest available GPU on the system
int getNumCores(const int major, const int minor) {
	const std::vector<std::pair<int, int>> mapping {
		{ 20, 32 }, // Fermi GF100
		{ 21, 48 }, // Fermi GF10x
		{ 30, 192 }, // Kepler GK10x
		{ 32, 192 }, // Kepler GK10x
		{ 35, 192 }, // Kepler GK11x
		{ 37, 192 }, // Kepler GK21x
	};

	const int combined_version = major * 10 + minor;

	// all the Maxwells have 128 so far, so making that the default
	int num_cores = 128;

	for (auto&& pair : mapping) {
		if (pair.first == combined_version) {
			num_cores = pair.second;
			break;
		}
	}

	return num_cores;
}

// approximate the GFLOPS of compute-enabled CUDA devices
// and return the fastest one
bool getHighestGFLOPSDevice(int& ID) {
	cudaDeviceProp deviceProp;
	int num_devices = 0;
	if (checkForCUDAError(cudaGetDeviceCount(&num_devices))) return false;

	if (num_devices == 0) {
		std::cerr << "ERROR: CUDA reports no capable devices on this system. Aborting." << std::endl;
		return false;
	}

	int devices_prohibited = 0;
	int best_SM_arch = 0;
	for (int i = 0; i < num_devices; ++i) {
		if (checkForCUDAError(cudaGetDeviceProperties(&deviceProp, i))) return false;
		if (deviceProp.computeMode == cudaComputeModeProhibited) {
			++devices_prohibited;
		}
		else {
			if (deviceProp.major > 0 && deviceProp.major < 9999) {
				best_SM_arch = std::max(best_SM_arch, deviceProp.major);
			}
		}
	}

	if (devices_prohibited == num_devices) {
		std::cerr << "ERROR: CUDA reports that all CUDA devices have compute mode prohibited. Aborting." << std::endl;
		return false;
	}

	size_t best_perf = 0;
	int best_i = 0;
	for (int i = 0; i < num_devices; ++i) {
		cudaGetDeviceProperties(&deviceProp, i);
		if (deviceProp.computeMode != cudaComputeModeProhibited) {
			size_t num_SMs = (deviceProp.major == 9999 && deviceProp.minor == 9999) ? 1 : getNumCores(deviceProp.major, deviceProp.minor);
			size_t perf = static_cast<size_t>(deviceProp.multiProcessorCount) * num_SMs * deviceProp.clockRate;
			if (perf  > best_perf) {
				if (best_SM_arch > 2) {
					if (deviceProp.major == best_SM_arch) {
						best_perf = perf;
						best_i = i;
					}
				}
				else {
					best_perf = perf;
					best_i = i;
				}
			}
		}
	}

	std::cout << "Autoselecting CUDA device with highest GFLOPS: Device " << best_i << std::endl;

	ID = best_i;
	return true;
}

bool KAero_CUDA::init(const int cuda_device_ID) {
	if (cuda_device_ID == AUTO_SELECT) {
		if (!getHighestGFLOPSDevice(cuda_device)) return false;
	}
	else {
		cuda_device = cuda_device_ID;
	}

	if (checkForCUDAError(cudaSetDevice(cuda_device))) return false;

	if (checkForCUDAError(cudaMalloc(&d_N, num_poly*sizeof(vec3)))) return false;
	if (checkForCUDAError(cudaMalloc(&d_precomp, num_poly*sizeof(double)))) return false;
	if (checkForCUDAError(cudaMalloc(&d_min_y, num_poly*sizeof(double)))) return false;
	if (checkForCUDAError(cudaMalloc(&d_max_y, num_poly*sizeof(double)))) return false;
	if (checkForCUDAError(cudaMalloc(&d_min_z, num_poly*sizeof(double)))) return false;
	if (checkForCUDAError(cudaMalloc(&d_max_z, num_poly*sizeof(double)))) return false;
	if (checkForCUDAError(cudaMalloc(&d_totals, 4 * sizeof(double)))) return false;
	if (checkForCUDAError(cudaMalloc(&d_P_rot, total_pts * sizeof(vec3)))) return false;
	if (checkForCUDAError(cudaMalloc(&d_CM_R, sizeof(vec3)))) return false;

	cudaDeviceProp device_prop;
	if (checkForCUDAError(cudaGetDeviceProperties(&device_prop, cuda_device))) return false;

	std::cout << "CUDA KAero initialized on device " << cuda_device << ": " << device_prop.name << std::endl;

	return true;
}

KAero_CUDA::KAero_CUDA(const double linear_pitch, const int num_polygons, const vec3* const poly, const vec3 sat_CM) :
pitch(linear_pitch),
f_scalar(2.0*linear_pitch*linear_pitch),
num_poly(num_polygons),
total_pts(num_polygons * NUM_VTX),
CM(sat_CM),
P_s(new vec3[num_polygons * NUM_VTX]),
P_rot(new vec3[num_polygons * NUM_VTX]) {

	// copy polygon data into internal storage
	for (int i = 0; i < total_pts; ++i) {
		P_s[i] = poly[i];
	}

}

KAero_CUDA::~KAero_CUDA() {
	cudaDeviceReset();
}

// kernel to precompute normals and portion of collision location finding
template<int block_size>
__global__ void precompute(vec3* __restrict__ const d_P_rot,
	vec3* __restrict__ const d_N, double* __restrict__ const d_precomp,
	double* __restrict__ const d_min_y, double* __restrict__ const d_max_y,
	double* __restrict__ const d_min_z, double* __restrict__ const d_max_z,
	double* __restrict__ const d_totals, int num_poly) {

	int i = threadIdx.x;

	// partition shared memory into 4 chunks of num_poly doublesm
	// for storing the min y's, max y's, min z's, and max z's, respectively
	extern __shared__ double s_min_y[];
	double* s_max_y = s_min_y + num_poly;
	double* s_min_z = s_min_y + 2 * num_poly;
	double* s_max_z = s_min_y + 3 * num_poly;

	double l_min_y = DBL_MAX;
	double l_max_y = -DBL_MAX;
	double l_min_z = DBL_MAX;
	double l_max_z = -DBL_MAX;

	if (i < num_poly) {
		// precompute some info for speed,
		// including panel normals and some of collision location math
		vec3* P = d_P_rot + i*NUM_VTX;
		d_N[i] = glm::normalize(glm::cross(P[2] - P[1], P[2] - P[3]));
		d_precomp[i] = glm::dot(d_N[i], P[0]);

		l_min_y = l_max_y = P[0].y;
		l_min_z = l_max_z = P[0].z;

		// find mins and maxes of y and z of each poly
#pragma unroll
		for (int j = 1; j < NUM_VTX; ++j) {
			l_min_y = (P[j].y < l_min_y) ? P[j].y : l_min_y;
			l_max_y = (P[j].y > l_max_y) ? P[j].y : l_max_y;
			l_min_z = (P[j].z < l_min_z) ? P[j].z : l_min_z;
			l_max_z = (P[j].z > l_max_z) ? P[j].z : l_max_z;
		}

		d_min_y[i] = l_min_y;
		d_max_y[i] = l_max_y;
		d_min_z[i] = l_min_z;
		d_max_z[i] = l_max_z;

	}

	s_min_y[i] = l_min_y;
	s_max_y[i] = l_max_y;
	s_min_z[i] = l_min_z;
	s_max_z[i] = l_max_z;
	__syncthreads();

	// fully unrolled reduction above 64
	// see CUDA reduction documentation or KPS reseach paper for more info

	if (block_size >= 512 && i < 256) {
		s_min_y[i] = l_min_y = s_min_y[i + 256] < l_min_y ? s_min_y[i + 256] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 256] > l_max_y ? s_max_y[i + 256] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 256] < l_min_z ? s_min_z[i + 256] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 256] > l_max_z ? s_max_z[i + 256] : l_max_z;
	}
	__syncthreads();

	if (block_size >= 256 && i < 128) {
		s_min_y[i] = l_min_y = s_min_y[i + 128] < l_min_y ? s_min_y[i + 128] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 128] > l_max_y ? s_max_y[i + 128] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 128] < l_min_z ? s_min_z[i + 128] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 128] > l_max_z ? s_max_z[i + 128] : l_max_z;
	}
	__syncthreads();

	if (block_size >= 128 && i < 64) {
		s_min_y[i] = l_min_y = s_min_y[i + 64] < l_min_y ? s_min_y[i + 64] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 64] > l_max_y ? s_max_y[i + 64] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 64] < l_min_z ? s_min_z[i + 64] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 64] > l_max_z ? s_max_z[i + 64] : l_max_z;
	}
	__syncthreads();

	// at this point 64 remain
	// CUDA architectures above 3.0 support the warp shuffle operation to do the last warp
	// very efficiently.
	//
	// see CUDA warp shuffle documentation for more info
	//
	// if available, warp shuffle is used.
	// if not, the remaining steps are manually unrolled just like above

#if (__CUDA_ARCH__ >= 300)
	if (i < 32) {
		if (block_size >= 64) {
			l_min_y = s_min_y[i + 32] < l_min_y ? s_min_y[i + 32] : l_min_y;
			l_max_y = s_max_y[i + 32] > l_max_y ? s_max_y[i + 32] : l_max_y;
			l_min_z = s_min_z[i + 32] < l_min_z ? s_min_z[i + 32] : l_min_z;
			l_max_z = s_max_z[i + 32] > l_max_z ? s_max_z[i + 32] : l_max_z;
		}

		for (int offset = warpSize / 2; offset > 0; offset /= 2) {
			double shuffled_min_y = __shfl_down(l_min_y, offset);
			l_min_y = shuffled_min_y < l_min_y ? shuffled_min_y : l_min_y;

			double shuffled_max_y = __shfl_down(l_max_y, offset);
			l_max_y = shuffled_max_y > l_max_y ? shuffled_max_y : l_max_y;

			double shuffled_min_z = __shfl_down(l_min_z, offset);
			l_min_z = shuffled_min_z < l_min_z ? shuffled_min_z : l_min_z;

			double shuffled_max_z = __shfl_down(l_max_z, offset);
			l_max_z = shuffled_max_z > l_max_z ? shuffled_max_z : l_max_z;
		}
	}

#else
	// fully unrolled reduction within single warp
	if (block_size >= 64 && i < 32)
	{
		s_min_y[i] = l_min_y = s_min_y[i + 32] < l_min_y ? s_min_y[i + 32] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 32] > l_max_y ? s_max_y[i + 32] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 32] < l_min_z ? s_min_z[i + 32] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 32] > l_max_z ? s_max_z[i + 32] : l_max_z;
	}
	__syncthreads();

	if (block_size >= 32 && i < 16)
	{
		s_min_y[i] = l_min_y = s_min_y[i + 16] < l_min_y ? s_min_y[i + 16] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 16] > l_max_y ? s_max_y[i + 16] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 16] < l_min_z ? s_min_z[i + 16] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 16] > l_max_z ? s_max_z[i + 16] : l_max_z;
	}
	__syncthreads();

	if (block_size >= 16 && i < 8)
	{
		s_min_y[i] = l_min_y = s_min_y[i + 8] < l_min_y ? s_min_y[i + 8] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 8] > l_max_y ? s_max_y[i + 8] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 8] < l_min_z ? s_min_z[i + 8] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 8] > l_max_z ? s_max_z[i + 8] : l_max_z;
	}
	__syncthreads();

	if (block_size >= 8 && i < 4)
	{
		s_min_y[i] = l_min_y = s_min_y[i + 4] < l_min_y ? s_min_y[i + 4] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 4] > l_max_y ? s_max_y[i + 4] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 4] < l_min_z ? s_min_z[i + 4] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 4] > l_max_z ? s_max_z[i + 4] : l_max_z;
	}
	__syncthreads();

	if (block_size >= 4 && i < 2)
	{
		s_min_y[i] = l_min_y = s_min_y[i + 2] < l_min_y ? s_min_y[i + 2] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 2] > l_max_y ? s_max_y[i + 2] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 2] < l_min_z ? s_min_z[i + 2] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 2] > l_max_z ? s_max_z[i + 2] : l_max_z;
	}
	__syncthreads();

	if (block_size >= 2 && i < 1)
	{
		s_min_y[i] = l_min_y = s_min_y[i + 1] < l_min_y ? s_min_y[i + 1] : l_min_y;
		s_max_y[i] = l_max_y = s_max_y[i + 1] > l_max_y ? s_max_y[i + 1] : l_max_y;
		s_min_z[i] = l_min_z = s_min_z[i + 1] < l_min_z ? s_min_z[i + 1] : l_min_z;
		s_max_z[i] = l_max_z = s_max_z[i + 1] > l_max_z ? s_max_z[i + 1] : l_max_z;
	}
	__syncthreads();
#endif

	// thread 0 contains this block's final result
	// store it out to global mem
	if (i == 0) {
		d_totals[0] = l_min_y;
		d_totals[1] = l_max_y;
		d_totals[2] = l_min_z;
		d_totals[3] = l_max_z;
	}

}

// kernel to collide particles and determine resultant forces and/or torques
__global__ void collide(vec3* __restrict__ const d_P_rot,
	vec3* __restrict__ const d_N, double* __restrict__ const d_precomp,
	double* __restrict__ const d_min_y, double* __restrict__ const d_max_y,
	double* __restrict__ const d_min_z, double* __restrict__ const d_max_z,
	int num_k, int num_m, int min_k, int min_m, double rho, double v_mag2,
	vec3* __restrict__ const d_CM_R, vec3* __restrict__ const d_block_sums_f,
	vec3* __restrict__ const d_block_sums_t, double pitch,
	int num_poly, double f_scalar) {

	int tidx = THREADS*blockIdx.x + threadIdx.x;
	int tidy = THREADS*blockIdx.y + threadIdx.y;

	// partition shared memory into forces and torques
	extern __shared__ vec3 sdata_f[];
	vec3* sdata_t = sdata_f + THREADS*THREADS;

	int sid = threadIdx.y*THREADS + threadIdx.x;

	vec3 force, torque;

	if (tidx < num_k && tidy < num_m) {
		double y = (min_k + tidx) * pitch;
		double z = (min_m + tidy) * pitch;

		double best_x = NONE;
		vec3 best_N;

		// check for collisions against each poly:
		for (int p = 0; p < num_poly; ++p) {
			// bail early if outside bounding box of this particular poly...
			if (y < d_min_y[p] || y > d_max_y[p] || z < d_min_z[p] || z > d_max_z[p]) continue;

			// ...otherwise, perform point-in-polygon anaylsis.
			// Looks scary but isn't too bad. Pretend coordinate system
			// is such that test point is at origin. Send the point to the right (+x)
			// and flip a flag each time it crosses a segment of the polygon. Flag will
			// end up flipped if inside polygon because it will cross an odd number of
			// segments. This works even for concave polygons.
			// For each segment of the polygon:
			//	if both ends are above or below x axis (i.e. share same sign), not a crossing
			//		otherwise, if both ends are left of the y axis, not a crossing
			//			otherwise, we do the slow bit (but usually don't have to due to the above):
			//			see if the segment intersects the x axis right of 0. if it does, crossing!
			// See KPS research paper for more.
			vec3* P = d_P_rot + p*NUM_VTX;
			int j = NUM_VTX - 1;
			int odd_nodes = 0;

#pragma unroll
			for (int i = 0; i < NUM_VTX; ++i) {
				odd_nodes ^= (((((P[i].z < z && P[j].z >= z) || (P[j].z < z && P[i].z >= z))
					&& (P[i].y <= y || P[j].y <= y)))
					&& ((P[i].y + (z - P[i].z) / (P[j].z - P[i].z)*(P[j].y - P[i].y) < y)));
				j = i;
			}

			// if inside polygon, compute collision location
			if (odd_nodes && d_N[p].x) {
				double x = (d_precomp[p] - (d_N[p].y*y + d_N[p].z*z)) / d_N[p].x;

				// and if it's the best (first) collision, update best_x
				if (x > best_x) {
					best_x = x;
					best_N = d_N[p];
				}
			}
		}

		// if a collision occurred
		if (best_x > NONE + PAD) {
			// see Equation 62 in KPS Research paper
			force = f_scalar*rho*best_N*(-v_mag2*best_N.x*fabs(best_N.x));

			// see Equation 63 in KPS Research paper
			torque = glm::cross(vec3{ best_x, y, z } -(*d_CM_R), force);
		}
	}

	sdata_f[sid] = force;
	sdata_t[sid] = torque;
	__syncthreads();

	// sum reduction in shared mem (block size is known to be 256)
	if (sid < 128) {
		sdata_f[sid] = force = force + sdata_f[sid + 128];
		sdata_t[sid] = torque = torque + sdata_t[sid + 128];
	}
	__syncthreads();

	if (sid < 64) {
		sdata_f[sid] = force = force + sdata_f[sid + 64];
		sdata_t[sid] = torque = torque + sdata_t[sid + 64];
	}
	__syncthreads();

	// at this point 64 remain
	// CUDA architectures above 3.0 support the warp shuffle operation to do the last warp
	// very efficiently.
	//
	// see CUDA warp shuffle documentation for more info
	//
	// if available, warp shuffle is used.
	// if not, the remaining steps are manually unrolled just like above

#if (__CUDA_ARCH__ >= 300)
	if (sid < 32) {
		force += sdata_f[sid + 32];
		torque += sdata_t[sid + 32];

		for (int offset = warpSize / 2; offset > 0; offset /= 2) {
			force.x += __shfl_down(force.x, offset);
			force.y += __shfl_down(force.y, offset);
			force.z += __shfl_down(force.z, offset);

			torque.x += __shfl_down(torque.x, offset);
			torque.y += __shfl_down(torque.y, offset);
			torque.z += __shfl_down(torque.z, offset);
		}

	}
#else
	// fully unrolled reduction within single warp
	if (sid < 32) {
		sdata_f[sid] = force = force + sdata_f[sid + 32];
		sdata_t[sid] = torque = torque + sdata_t[sid + 32];
	}
	__syncthreads();

	if (sid < 16) {
		sdata_f[sid] = force = force + sdata_f[sid + 16];
		sdata_t[sid] = torque = torque + sdata_t[sid + 16];
	}
	__syncthreads();

	if (sid < 8) {
		sdata_f[sid] = force = force + sdata_f[sid + 8];
		sdata_t[sid] = torque = torque + sdata_t[sid + 8];
	}
	__syncthreads();

	if (sid < 4) {
		sdata_f[sid] = force = force + sdata_f[sid + 4];
		sdata_t[sid] = torque = torque + sdata_t[sid + 4];
	}
	__syncthreads();

	if (sid < 2) {
		sdata_f[sid] = force = force + sdata_f[sid + 2];
		sdata_t[sid] = torque = torque + sdata_t[sid + 2];
	}

	__syncthreads();

	if (sid < 1) {
		sdata_f[sid] = force = force + sdata_f[sid + 1];
		sdata_t[sid] = torque = torque + sdata_t[sid + 1];
	}
	__syncthreads();
#endif

	// thread 0 contains this block's final result
	// store it out to global mem
	if (sid == 0) {
		d_block_sums_f[blockIdx.x*gridDim.y + blockIdx.y] = force;
		d_block_sums_t[blockIdx.x*gridDim.y + blockIdx.y] = torque;
	}

}

void KAero_CUDA::aer(vec3& f, vec3& t, const double rho, const vec3& v) {

	// zero out 'f' and 't'
	f = t = vec3();

	double v_mag2 = glm::length2(v);
	double v_mag = sqrt(v_mag2);

	// --- ROTATOR SETUP ---
	// rotates entire satellite such that velocity is in +x in new frame, i.e. relative wind
	// direction is -x
	//
	// Uses Rodrigues' rotation formula for speed
	// the ternary handles cases where velocity is entirely in +x or entirely in -x already
	//
	// does NOT handle zero velocity, but if that's the case in an *orbital simulation*
	// you have bigger problems to worry about
	double cos_theta = v.x / v_mag;
	double sin_theta = sin(acos(cos_theta));
	vec3 k = sin_theta ? vec3{ 0.0, v.z, -v.y } / (v_mag * sin_theta) : vec3();
	vec3 k_times_1_minus_cos_theta = k*(1 - cos_theta);
	// --- /ROTATOR SETUP ---

	// Use of ROTATOR (Rodrigues' formula) to rotate satellite CM into new frame
	vec3 CM_R = cos_theta*CM + sin_theta*glm::cross(k, CM) + k_times_1_minus_cos_theta*glm::dot(k, CM);

	// push CM to GPU
	// NOTE that no error checks are being done on CUDA returns;
	// they are during initialization, but not in the hot path
	cudaMemcpy(d_CM_R, &CM_R, sizeof(vec3), cudaMemcpyHostToDevice);

	for (int i = 0; i < total_pts; ++i) {
		// Use of ROTATOR (Rodrigues' formula) to rotate polygons into new frame
		P_rot[i] = cos_theta*P_s[i] + sin_theta*glm::cross(k, P_s[i]) + k_times_1_minus_cos_theta*glm::dot(k, P_s[i]);
	}

	// push rotated polygons to GPU (proved faster than doing the rotation on GPU)
	cudaMemcpy(d_P_rot, P_rot, total_pts*sizeof(vec3), cudaMemcpyHostToDevice);


	// --- PRECOMPUTE ---
	// precomp includes a reduction for total satellite bounding box,
	// so launch a power of 2 number of threads
	//
	// see KPS reseach paper and/or CUDA documentation on reductions for more info
	int precomp_threads = nextHigherPow2(num_poly);

	// increase kernel performance by templating on thread count
	// polygon number doesn't change that much... currently supporting < 256
	//
	// Blocks: 1
	// Threads: next power of two above num_poly
	// Shared memory: 4 doubles (min_y, max_y, min_z, max_z) for each thread,
	//				EXCEPT not less than 64 total groups of 4 doubles
	//				because 2 warps (of 32) are combined before shuffling begins
	//				within each warp
	//
	// These switch statements are out of order intentionally; they are arranged
	// in order of most likelihood to assist the branch predictor
	// (usually ~10 polygons, so 16 threads... sometimes ~6 polygons, so 8, etc...)
	//
	switch (precomp_threads) {
	case 16:
		precompute<16> << <1, 16, 4 * 64 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 8:
		precompute<8> << <1, 8, 4 * 64 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 4:
		precompute<4> << <1, 4, 4 * 64 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 32:
		precompute<32> << <1, 32, 4 * 64 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 2:
		precompute<2> << <1, 2, 4 * 64 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 64:
		precompute<64> << <1, 64, 4 * 64 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 1:
		precompute<1> << <1, 1, 4 * 64 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 128:
		precompute<128> << <1, 128, 4 * 128 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 256:
		precompute<256> << <1, 256, 4 * 256 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	case 512:
		precompute<512> << <1, 512, 4 * 512 * sizeof(double) >> >(d_P_rot, d_N,
			d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, d_totals, num_poly);
		break;
	default:
		std::cerr
			<< "ERROR: Too many polygons! This should never occur, as it should be caught" << std::endl
			<< "in the polygon input stage." << std::endl;
		break;
	}

	cudaMemcpy(totals, d_totals, 4 * sizeof(double), cudaMemcpyDeviceToHost);
	// --- /PRECOMPUTE ---


	// CPU: use total bounding box information to get rectangular area to simulate
	const int min_k = static_cast<int>((totals[0] + PAD) / pitch);
	const int max_k = static_cast<int>((totals[1] - PAD) / pitch);
	const int min_m = static_cast<int>((totals[2] + PAD) / pitch);
	const int max_m = static_cast<int>((totals[3] - PAD) / pitch);

	// determine number of particles required in each direction
	const int num_k = max_k - min_k + 1;
	const int num_m = max_m - min_m + 1;

	// determine required number of collide blocks in 2-D grid
	const int collide_blks_x = (num_k - 1) / THREADS + 1;
	const int collide_blks_y = (num_m - 1) / THREADS + 1;

	const int collide_blks_total = collide_blks_x * collide_blks_y;

	// dynamically allocate device memory for the force and torque sums from each block
	cudaMalloc(&d_block_sums_f, collide_blks_total*sizeof(vec3));
	cudaMalloc(&d_block_sums_t, collide_blks_total*sizeof(vec3));

	// dynamically allocate host memory for the force and torque sums from each block
	block_sums_f = new vec3[collide_blks_total];
	block_sums_t = new vec3[collide_blks_total];

	// collide!
	// shared memory is 2 vec3's (one for force, one for torque) per thread
	// NOTE that collide also performs a reduction, but the thread count is fixed at THREADS, which is 16,
	// already a multiple of 2.
	collide<<<dim3(collide_blks_x, collide_blks_y), dim3(THREADS, THREADS), 2*THREADS*THREADS*sizeof(vec3)>>>(d_P_rot,
		d_N, d_precomp, d_min_y, d_max_y, d_min_z, d_max_z, num_k, num_m, min_k, min_m, rho, v_mag2,
		d_CM_R, d_block_sums_f, d_block_sums_t, pitch, num_poly, f_scalar);

	// retrieve force block sums
	cudaMemcpy(block_sums_f, d_block_sums_f, collide_blks_total*sizeof(vec3), cudaMemcpyDeviceToHost);

	// free force block sum device dynamic allocation
	cudaFree(d_block_sums_f);

	// retrieve torque block sums
	cudaMemcpy(block_sums_t, d_block_sums_t, collide_blks_total*sizeof(vec3), cudaMemcpyDeviceToHost);

	// free torque block sum device dynamic allocation
	cudaFree(d_block_sums_t);

	// accumulate the block sums
	// there won't be too many; Kahan summation is not required
	for (int i = 0; i < collide_blks_total; ++i) {
		f += block_sums_f[i];
		t += block_sums_t[i];
	}

	// Use of ROTATOR (Rodrigues' formula) to rotate summed force and torque BACK TO BODY FRAME
	// (note the negative sign)
	f = cos_theta*f - sin_theta*glm::cross(k, f) + k_times_1_minus_cos_theta*glm::dot(k, f);
	t = cos_theta*t - sin_theta*glm::cross(k, t) + k_times_1_minus_cos_theta*glm::dot(k, t);

	// free force block sum host dynamic allocation
	delete[] block_sums_f;

	// free torque block sum host dynamic allocation
	delete[] block_sums_t;

}
