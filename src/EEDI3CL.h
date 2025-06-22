#pragma once

#include <cfloat>

#define CL_TARGET_OPENCL_VERSION 300

#define BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION
#define BOOST_COMPUTE_HAVE_THREAD_LOCAL
#define BOOST_COMPUTE_THREAD_SAFE
#define BOOST_COMPUTE_USE_OFFLINE_CACHE

#include <boost/compute/core.hpp>
#include <boost/compute/utility/dim.hpp>

#include "avs_c_api_loader.hpp"

AVS_FORCEINLINE void* aligned_malloc(size_t size, size_t align)
{
    void* result = [&]() {
#ifdef _WIN32
        return _aligned_malloc(size, align);
#else
        if (posix_memalign(&result, align, size))
            return result = nullptr;
        else
            return result;
#endif
    }();

    return result;
}

AVS_FORCEINLINE void aligned_free(void* ptr)
{
#ifdef _WIN32
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

template<typename T>
struct aligned_array_deleter
{
    void operator()(T* ptr) const noexcept
    {
        if (ptr)
            aligned_free(ptr);
    }
};

template<typename T>
using aligned_unique_ptr = std::unique_ptr<T[], aligned_array_deleter<T>>;

template<typename T>
inline aligned_unique_ptr<T> make_unique_aligned_array(size_t num_elements, size_t alignment)
{
    if (num_elements == 0)
        return aligned_unique_ptr<T>(nullptr);

    T* ptr = static_cast<T*>(aligned_malloc(num_elements * sizeof(T), alignment));

    return aligned_unique_ptr<T>(ptr);
}

struct EEDI3CLData
{
    avs_helpers::avs_clip_ptr sclip;
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
    aligned_unique_ptr<float> pcosts;
    aligned_unique_ptr<int> pbackt;
    std::unique_ptr<int[]> fpath;
    std::unique_ptr<int[]> dmap;
    std::string err;

    void (*filter)(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, const int field_n, bool use_dh,
        EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
};

template<typename T>
void copyPlane(void* dstp_, const int dstStride, const void* srcp_, const int srcStride, const int width, const int height) noexcept;
template<typename T>
void copyPad(const AVS_VideoFrame* src, AVS_VideoFrame* dst, const int plane, const int plane0, const int off, const bool dh,
    AVS_ScriptEnvironment* __restrict env) noexcept;

template<typename T>
void filterCL_sse2(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, const int field_n, bool use_dh,
    EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
template<typename T>
void filterCL_avx2(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, const int field_n, bool use_dh,
    EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
template<typename T>
void filterCL_avx512(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, const int field_n, bool use_dh,
    EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
