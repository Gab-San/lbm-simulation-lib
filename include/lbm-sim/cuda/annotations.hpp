#ifndef __LBM_SIM_BACKEND_CUDA_METADATA_HPP
#define __LBM_SIM_BACKEND_CUDA_METADATA_HPP

#ifdef __CUDACC__
#define LBM_HD_FUNC __host__ __device__
#else
#define LBM_HD_FUNC
#endif

#endif // __LBM_SIM_BACKEND_CUDA_METADATA_HPP
