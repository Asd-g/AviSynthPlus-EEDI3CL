#pragma once

#include "avisynth_c.h"
#include "VCL2/vectorclass.h"

#define BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION
#define BOOST_COMPUTE_HAVE_THREAD_LOCAL
#define BOOST_COMPUTE_THREAD_SAFE
#define BOOST_COMPUTE_USE_OFFLINE_CACHE
#include <boost/compute/core.hpp>
#include <boost/compute/utility/dim.hpp>

struct EEDI3CLData
{
    AVS_FilterInfo* fi;
    AVS_VideoInfo vi_pad;
    AVS_Clip* sclip;
    int field;
    int mdis;
    int vcheck;
    bool dh;
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
    boost::compute::command_queue queue;
    boost::compute::kernel calculateConnectionCosts;
    boost::compute::image2d src;
    boost::compute::buffer ccosts;
    float* pcosts;
    int* pbackt;
    int* fpath;
    int* dmap;
    int* tline;
    std::unique_ptr<char[]> err;

    void (*filter)(const AVS_VideoFrame*, const AVS_VideoFrame*, AVS_VideoFrame*, AVS_VideoFrame**, const int, const EEDI3CLData* const __restrict);
};

template<typename T>
void copyPlane(void* __restrict dstp_, const int dstStride, const void* srcp_, const int srcStride, const int width, const int height) noexcept;
template<typename T>
void copyPad(const AVS_VideoFrame* src, AVS_VideoFrame* dst, const int plane, const int plane0, const int off, const bool dh, AVS_ScriptEnvironment* env) noexcept;
template<typename T>
inline void interpolate(const T* src3p, const T* src1p, const T* src1n, const T* src3n, const bool* bmask, const int* fpath, int* __restrict dmap, T* __restrict dstp, const int width, const bool ucubic, const int peak) noexcept;
template<typename T>
void vCheck(const T* srcp, const T* scpp, T* __restrict dstp, const int* dmap, void* _tline, const int field_n, const int dstWidth, const int srcHeight, const int srcStride, const int dstStride, const int vcheck,
    const float vthresh2, const float rcpVthresh0, const float rcpVthresh1, const float rcpVthresh2, const int peak) noexcept;

template<typename T>
void filterCL_sse2(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, AVS_VideoFrame** pad, const int field_n, const EEDI3CLData* const __restrict d);
