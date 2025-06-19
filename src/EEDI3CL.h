#pragma once

#include <cfloat>

#include <avisynth_c.h>

#define CL_TARGET_OPENCL_VERSION 300

#define BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION
#define BOOST_COMPUTE_HAVE_THREAD_LOCAL
#define BOOST_COMPUTE_THREAD_SAFE
#define BOOST_COMPUTE_USE_OFFLINE_CACHE

#include <boost/compute/core.hpp>
#include <boost/compute/utility/dim.hpp>

struct EEDI3CLData
{
    AVS_Clip* sclip;
    int field;
    int mdis;
    int vcheck;
    bool dh;
    bool dw;
    bool process[4];
    bool ucubic;
    float gamma;
    float vthresh2;
    int peak;
    int vectorSize;
    int tpitch;
    int mdisVector;
    int tpitchVector;
    float rcpVthresh0;
    float rcpVthresh1;
    float rcpVthresh2;
    cl_image_format clImageFormat;
    boost::compute::kernel applyInterpolationKernel;
    boost::compute::kernel vCheckKernel;
    boost::compute::kernel transpose_u8_kernel;
    boost::compute::kernel transpose_u16_kernel;
    boost::compute::kernel transpose_f32_kernel;
    boost::compute::command_queue queue;
    boost::compute::kernel calculateConnectionCosts;
    boost::compute::image2d src;
    boost::compute::image2d dst;
    boost::compute::image2d vcheck_tmp;
    boost::compute::buffer ccosts;
    boost::compute::buffer fpath_gpu;
    boost::compute::buffer dmap_gpu;
    float* pcosts;
    int* pbackt;
    int* fpath;
    int* dmap;
    std::string err;

    void (*filter)(const AVS_VideoFrame* __restrict, const AVS_VideoFrame* __restrict, AVS_VideoFrame* __restrict, const int, bool, EEDI3CLData* __restrict, const AVS_FilterInfo* __restrict);
};

template<typename T>
void copyPlane(void* __restrict dstp_, const int dstStride, const void* __restrict srcp_, const int srcStride, const int width, const int height) noexcept;
template<typename T>
void copyPad(const AVS_VideoFrame* __restrict src, AVS_VideoFrame* __restrict dst, const int plane, const int plane0, const int off, const bool dh, AVS_ScriptEnvironment* __restrict env) noexcept;

template<typename T>
void filterCL_sse2(const AVS_VideoFrame* __restrict src, const AVS_VideoFrame* __restrict scp, AVS_VideoFrame* __restrict dst, const int field_n, bool use_dh, EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
template<typename T>
void filterCL_avx2(const AVS_VideoFrame* __restrict src, const AVS_VideoFrame* __restrict scp, AVS_VideoFrame* __restrict dst, const int field_n, bool use_dh, EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
template<typename T>
void filterCL_avx512(const AVS_VideoFrame* __restrict src, const AVS_VideoFrame* __restrict scp, AVS_VideoFrame* __restrict dst, const int field_n, bool use_dh, EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
