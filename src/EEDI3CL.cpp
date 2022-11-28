/*
**   Another VapourSynth port by HolyWu
**
**   eedi3 (enhanced edge directed interpolation 3). Works by finding the
**   best non-decreasing (non-crossing) warping between two lines according to
**   a cost functional. Doesn't really have anything to do with eedi2 aside
**   from doing edge-directed interpolation (they use different techniques).
**
**   Copyright (C) 2010 Kevin Stone - some part by Laurent de Soras, 2013
**
**   This program is free software; you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation; either version 2 of the License, or
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License
**   along with this program; if not, write to the Free Software
**   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <locale>
#include <memory>
#include <sstream>
#include <string>

#include "EEDI3CL.cl"
#include "EEDI3CL.h"


AVS_FORCEINLINE void* aligned_malloc(size_t size, size_t align)
{
    void* result = [&]() {
#ifdef _MSC_VER 
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
#ifdef _MSC_VER 
    _aligned_free(ptr);
#else 
    free(ptr);
#endif
}

template<typename T>
void copyPlane(void* __restrict dstp_, const int dstStride, const void* srcp_, const int srcStride, const int width, const int height) noexcept
{
    const T* srcp = reinterpret_cast<const T*>(srcp_);
    T* __restrict dstp = reinterpret_cast<T*>(dstp_);

    for (int y{ 0 }; y < height; ++y)
    {
        for (int x{ 0 }; x < width; ++x)
            dstp[x] = srcp[x];

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
void copyPad(const AVS_VideoFrame* src, AVS_VideoFrame* dst, const int plane, const int plane0, const int off, const bool dh, AVS_ScriptEnvironment* env) noexcept
{
    const int srcWidth{ static_cast<int>(avs_get_row_size_p(src, plane) / sizeof(T)) };
    const int dstWidth{ static_cast<int>(avs_get_row_size_p(dst, plane0) / sizeof(T)) };
    const int srcHeight{ avs_get_height_p(src, plane) };
    const int dstHeight{ avs_get_height_p(dst, plane0) };
    const int srcStride{ static_cast<int>(avs_get_pitch_p(src, plane) / sizeof(T)) };
    const int dstStride{ static_cast<int>(avs_get_pitch_p(dst, plane0) / sizeof(T)) };
    const T* srcp{ reinterpret_cast<const T*>(avs_get_read_ptr_p(src, plane)) };
    T* __restrict dstp{ reinterpret_cast<T*>(avs_get_write_ptr_p(dst, plane0)) };

    if (!dh)
        copyPlane<T>(dstp + dstStride * (4 + off) + 12, dstStride * 2, srcp + srcStride * off, srcStride * 2, srcWidth, srcHeight / 2);
    else
        copyPlane<T>(dstp + dstStride * (4 + off) + 12, dstStride * 2, srcp, srcStride, srcWidth, srcHeight);

    dstp += dstStride * (4 + off);

    for (int y{ 4 + off }; y < dstHeight - 4; y += 2)
    {
        for (int x{ 0 }; x < 12; ++x)
            dstp[x] = dstp[24 - x];

        for (int x{ dstWidth - 12 }, c{ 2 }; x < dstWidth; ++x, c += 2)
            dstp[x] = dstp[x - c];

        dstp += dstStride * 2;
    }

    dstp = reinterpret_cast<T*>(avs_get_write_ptr_p(dst, plane0));

    for (int y{ off }; y < 4; y += 2)
        memcpy(dstp + dstStride * y, dstp + dstStride * (8 - y), dstWidth * sizeof(T));

    for (int y{ dstHeight - 4 + off }, c{ 2 + 2 * off }; y < dstHeight; y += 2, c += 4)
        memcpy(dstp + dstStride * y, dstp + dstStride * (y - c), dstWidth * sizeof(T));
}

template<typename T>
inline void interpolate(const T* src3p, const T* src1p, const T* src1n, const T* src3n, const bool* bmask, const int* fpath, int* __restrict dmap, T* __restrict dstp, const int width, const bool ucubic, const int peak) noexcept
{
    if constexpr (std::is_integral_v<T>)
    {
        for (int x{ 0 }; x < width; ++x)
        {
            if (bmask && !bmask[x])
            {
                dmap[x] = 0;

                if (ucubic)
                    dstp[x] = std::min(std::max((9 * (src1p[x] + src1n[x]) - (src3p[x] + src3n[x]) + 8) >> 4, 0), peak);
                else
                    dstp[x] = (src1p[x] + src1n[x] + 1) >> 1;
            }
            else
            {
                const int dir{ fpath[x] };
                const int dir3{ dir * 3 };
                const int absDir3{ std::abs(dir3) };

                dmap[x] = dir;

                if (ucubic && x >= absDir3 && x <= width - 1 - absDir3)
                    dstp[x] = std::min(std::max((9 * (src1p[x + dir] + src1n[x - dir]) - (src3p[x + dir3] + src3n[x - dir3]) + 8) >> 4, 0), peak);
                else
                    dstp[x] = (src1p[x + dir] + src1n[x - dir] + 1) >> 1;
            }
        }
    }
    else
    {
        for (int x{ 0 }; x < width; ++x)
        {
            if (bmask && !bmask[x])
            {
                dmap[x] = 0;

                if (ucubic)
                    dstp[x] = 0.5625f * (src1p[x] + src1n[x]) - 0.0625f * (src3p[x] + src3n[x]);
                else
                    dstp[x] = (src1p[x] + src1n[x]) / 2.0f;
            }
            else
            {
                const int dir{ fpath[x] };
                const int dir3{ dir * 3 };
                const int absDir3{ std::abs(dir3) };

                dmap[x] = dir;

                if (ucubic && x >= absDir3 && x <= width - 1 - absDir3)
                    dstp[x] = 0.5625f * (src1p[x + dir] + src1n[x - dir]) - 0.0625f * (src3p[x + dir3] + src3n[x - dir3]);
                else
                    dstp[x] = (src1p[x + dir] + src1n[x - dir]) / 2.0f;
            }
        }
    }
}

template<typename T>
void vCheck(const T* srcp, const T* scpp, T* __restrict dstp, const int* dmap, void* _tline, const int field_n, const int dstWidth, const int srcHeight, const int srcStride, const int dstStride, const int vcheck,
    const float vthresh2, const float rcpVthresh0, const float rcpVthresh1, const float rcpVthresh2, const int peak) noexcept
{
    if constexpr (std::is_integral_v<T>)
    {
        for (int y{ 4 + field_n }; y < srcHeight - 4; y += 2)
        {
            if (y >= 6 && y < srcHeight - 6)
            {
                const T* dst3p{ srcp - srcStride * 3 };
                const T* dst2p{ dstp - dstStride * 2 };
                const T* dst1p{ dstp - dstStride };
                const T* dst1n{ dstp + dstStride };
                const T* dst2n{ dstp + dstStride * 2 };
                const T* dst3n{ srcp + srcStride * 3 };
                T* __restrict tline{ reinterpret_cast<T*>(_tline) };

                for (int x{ 0 }; x < dstWidth; ++x)
                {
                    const int dirc{ dmap[x] };
                    const T cint{ static_cast<T>(scpp ? scpp[x] : std::min(std::max((9 * (dst1p[x] + dst1n[x]) - (dst3p[x] + dst3n[x]) + 8) >> 4, 0), peak)) };

                    if (dirc == 0)
                    {
                        tline[x] = cint;
                        continue;
                    }

                    const int dirt{ dmap[x - dstWidth] };
                    const int dirb{ dmap[x + dstWidth] };

                    if (std::max(dirc * dirt, dirc * dirb) < 0 || (dirt == dirb && dirt == 0))
                    {
                        tline[x] = cint;
                        continue;
                    }

                    const int it{ (dst2p[x + dirc] + dstp[x - dirc] + 1) >> 1 };
                    const int vt{ std::abs(dst2p[x + dirc] - dst1p[x + dirc]) + std::abs(dstp[x + dirc] - dst1p[x + dirc]) };
                    const int ib{ (dstp[x + dirc] + dst2n[x - dirc] + 1) >> 1 };
                    const int vb{ std::abs(dst2n[x - dirc] - dst1n[x - dirc]) + std::abs(dstp[x - dirc] - dst1n[x - dirc]) };
                    const int vc{ std::abs(dstp[x] - dst1p[x]) + std::abs(dstp[x] - dst1n[x]) };

                    const int d0{ std::abs(it - dst1p[x]) };
                    const int d1{ std::abs(ib - dst1n[x]) };
                    const int d2{ std::abs(vt - vc) };
                    const int d3{ std::abs(vb - vc) };

                    const int mdiff0{ (vcheck == 1) ? std::min(d0, d1) : (vcheck == 2 ? (d0 + d1 + 1) >> 1 : std::max(d0, d1)) };
                    const int mdiff1{ (vcheck == 1) ? std::min(d2, d3) : (vcheck == 2 ? (d2 + d3 + 1) >> 1 : std::max(d2, d3)) };

                    const float a0{ mdiff0 * rcpVthresh0 };
                    const float a1{ mdiff1 * rcpVthresh1 };
                    const float a2{ std::max((vthresh2 - std::abs(dirc)) * rcpVthresh2, 0.0f) };
                    const float a{ std::min(std::max({ a0, a1, a2 }), 1.0f) };

                    tline[x] = static_cast<T>((1.0f - a) * dstp[x] + a * cint);
                }

                memcpy(dstp, tline, dstWidth * sizeof(T));
            }

            srcp += srcStride * 2;
            if (scpp)
                scpp += dstStride * 2;
            dstp += dstStride * 2;
            dmap += dstWidth;
        }
    }
    else
    {
        for (int y{ 4 + field_n }; y < srcHeight - 4; y += 2)
        {
            if (y >= 6 && y < srcHeight - 6)
            {
                const float* dst3p{ srcp - srcStride * 3 };
                const float* dst2p{ dstp - dstStride * 2 };
                const float* dst1p{ dstp - dstStride };
                const float* dst1n{ dstp + dstStride };
                const float* dst2n{ dstp + dstStride * 2 };
                const float* dst3n{ srcp + srcStride * 3 };
                float* __restrict tline{ reinterpret_cast<float*>(_tline) };

                for (int x{ 0 }; x < dstWidth; ++x)
                {
                    const int dirc{ dmap[x] };
                    const float cint{ scpp ? scpp[x] : 0.5625f * (dst1p[x] + dst1n[x]) - 0.0625f * (dst3p[x] + dst3n[x]) };

                    if (dirc == 0)
                    {
                        tline[x] = cint;
                        continue;
                    }

                    const int dirt{ dmap[x - dstWidth] };
                    const int dirb{ dmap[x + dstWidth] };

                    if (std::max(dirc * dirt, dirc * dirb) < 0 || (dirt == dirb && dirt == 0))
                    {
                        tline[x] = cint;
                        continue;
                    }

                    const float it{ (dst2p[x + dirc] + dstp[x - dirc]) / 2.0f };
                    const float vt{ std::abs(dst2p[x + dirc] - dst1p[x + dirc]) + std::abs(dstp[x + dirc] - dst1p[x + dirc]) };
                    const float ib{ (dstp[x + dirc] + dst2n[x - dirc]) / 2.0f };
                    const float vb{ std::abs(dst2n[x - dirc] - dst1n[x - dirc]) + std::abs(dstp[x - dirc] - dst1n[x - dirc]) };
                    const float vc{ std::abs(dstp[x] - dst1p[x]) + std::abs(dstp[x] - dst1n[x]) };

                    const float d0{ std::abs(it - dst1p[x]) };
                    const float d1{ std::abs(ib - dst1n[x]) };
                    const float d2{ std::abs(vt - vc) };
                    const float d3{ std::abs(vb - vc) };

                    const float mdiff0{ (vcheck == 1) ? std::min(d0, d1) : (vcheck == 2 ? (d0 + d1) / 2.0f : std::max(d0, d1)) };
                    const float mdiff1{ (vcheck == 1) ? std::min(d2, d3) : (vcheck == 2 ? (d2 + d3) / 2.0f : std::max(d2, d3)) };

                    const float a0{ mdiff0 * rcpVthresh0 };
                    const float a1{ mdiff1 * rcpVthresh1 };
                    const float a2{ std::max((vthresh2 - std::abs(dirc)) * rcpVthresh2, 0.0f) };
                    const float a{ std::min(std::max({ a0, a1, a2 }), 1.0f) };

                    tline[x] = (1.0f - a) * dstp[x] + a * cint;
                }

                memcpy(dstp, tline, dstWidth * sizeof(float));
            }

            srcp += srcStride * 2;
            if (scpp)
                scpp += dstStride * 2;
            dstp += dstStride * 2;
            dmap += dstWidth;
        }
    }
}

template<typename T>
static void filterCL_c(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, AVS_VideoFrame** pad, const int field_n, const EEDI3CLData* const __restrict d)
{
    const int planes_y[4]{ AVS_PLANAR_Y, AVS_PLANAR_U, AVS_PLANAR_V, AVS_PLANAR_A };
    const int planes_r[4]{ AVS_PLANAR_R, AVS_PLANAR_G, AVS_PLANAR_B, AVS_PLANAR_A };
    const int* planes{ (avs_is_rgb(&d->fi->vi) ? planes_r : planes_y) };

    for (int plane{ 0 }; plane < avs_num_components(&d->fi->vi); ++plane)
    {
        if (d->process[plane])
        {
            copyPad<T>(src, pad[plane], planes[plane], planes[0], 1 - field_n, d->dh, d->fi->env);

            const int srcWidth{ static_cast<int>(avs_get_row_size_p(pad[plane], planes[0]) / sizeof(T)) };
            const int dstWidth{ static_cast<int>(avs_get_row_size_p(dst, planes[plane]) / sizeof(T)) };
            const int srcHeight{ avs_get_height_p(pad[plane], planes[0]) };
            const int dstHeight{ avs_get_height_p(dst, planes[plane]) };
            const int srcStride{ static_cast<int>(avs_get_pitch_p(pad[plane], planes[0]) / sizeof(T)) };
            const int dstStride{ static_cast<int>(avs_get_pitch_p(dst, planes[plane]) / sizeof(T)) };
            const T* _srcp{ reinterpret_cast<const T*>(avs_get_read_ptr_p(pad[plane], planes[0])) };
            T* __restrict _dstp{ reinterpret_cast<T*>(avs_get_write_ptr_p(dst, planes[plane])) };

            auto queue{ d->queue };
            auto calculateConnectionCosts{ d->calculateConnectionCosts };
            auto srcImage{ d->src };
            auto _ccosts{ d->ccosts };
            float* pcosts{ d->pcosts + d->mdis };
            int* pbackt{ d->pbackt + d->mdis };
            int* fpath{ d->fpath };
            int* _dmap{ d->dmap };
            int* tline{ d->tline };

            const size_t globalWorkSize[] = { static_cast<size_t>((dstWidth + 63) & -64), 1 };
            constexpr size_t localWorkSize[] = { 64, 1 };
            const int bufferSize{ static_cast<int>(dstWidth * d->tpitch * sizeof(cl_float)) };

            copyPlane<T>(_dstp + dstStride * (1 - field_n), dstStride * 2, _srcp + srcStride * (4 + 1 - field_n) + 12, srcStride * 2, dstWidth, dstHeight / 2);

            queue.enqueue_write_image(srcImage, boost::compute::dim(0, 0), boost::compute::dim(srcWidth, srcHeight), _srcp, avs_get_pitch_p(pad[plane], planes[0]));

            for (int y{ 4 + field_n }; y < srcHeight - 4; y += 2)
            {
                const int off{ (y - 4 - field_n) >> 1 };
                const T* srcp{ _srcp + srcStride * (4 + field_n + 2 * off) + 12 };
                T* dstp{ _dstp + dstStride * (field_n + 2 * off) };
                int* dmap{ _dmap + dstWidth * off };

                const T* src3p{ srcp - srcStride * 3 };
                const T* src1p{ srcp - srcStride };
                const T* src1n{ srcp + srcStride };
                const T* src3n{ srcp + srcStride * 3 };

                calculateConnectionCosts.set_args(srcImage, _ccosts, dstWidth, srcHeight - 4, y);
                queue.enqueue_nd_range_kernel(calculateConnectionCosts, 2, nullptr, globalWorkSize, localWorkSize);

                float* ccosts{ reinterpret_cast<float*>(queue.enqueue_map_buffer(_ccosts, CL_MAP_READ, 0, bufferSize)) + d->mdis };

                // calculate path costs
                *pcosts = *ccosts;
                for (int x{ 1 }; x < dstWidth; ++x)
                {
                    const float* tT{ ccosts + d->tpitch * x };
                    const float* ppT{ pcosts + d->tpitch * (x - 1) };
                    float* pT{ pcosts + d->tpitch * x };
                    int* piT{ pbackt + d->tpitch * (x - 1) };

                    const int umax{ std::min({ x, dstWidth - 1 - x, d->mdis }) };
                    const int umax2{ std::min({ x - 1, dstWidth - x, d->mdis }) };
                    for (int u{ -umax }; u <= umax; ++u)
                    {
                        int idx{ 0 };
                        float bval{ FLT_MAX };

                        for (int v{ std::max(-umax2, u - 1) }; v <= std::min(umax2, u + 1); ++v)
                        {
                            const double z{ ppT[v] + d->gamma * std::abs(u - v) };
                            const float ccost{ static_cast<float>(std::min(z, FLT_MAX * 0.9)) };
                            if (ccost < bval)
                            {
                                bval = ccost;
                                idx = v;
                            }
                        }

                        const double z{ bval + tT[u] };
                        pT[u] = static_cast<float>(std::min(z, FLT_MAX * 0.9));
                        piT[u] = idx;
                    }
                }

                // backtrack
                fpath[dstWidth - 1] = 0;
                for (int x{ dstWidth - 2 }; x >= 0; --x)
                    fpath[x] = pbackt[d->tpitch * x + fpath[x + 1]];

                interpolate<T>(src3p, src1p, src1n, src3n, nullptr, fpath, dmap, dstp, dstWidth, d->ucubic, d->peak);

                queue.enqueue_unmap_buffer(_ccosts, ccosts - d->mdis);
            }

            if (d->vcheck)
            {
                const T* srcp{ _srcp + srcStride * (4 + field_n) + 12 };
                const T* scpp{ nullptr };
                if (d->sclip)
                    scpp = reinterpret_cast<const T*>(avs_get_read_ptr_p(scp, planes[plane])) + dstStride * field_n;
                T* dstp{ _dstp + dstStride * field_n };

                vCheck<T>(srcp, scpp, dstp, _dmap, tline, field_n, dstWidth, srcHeight, srcStride, dstStride, d->vcheck, d->vthresh2, d->rcpVthresh0, d->rcpVthresh1, d->rcpVthresh2, d->peak);
            }
        }
    }
}

/* multiplies and divides a rational number, such as a frame duration, in place and reduces the result */
AVS_FORCEINLINE void muldivRational(unsigned* num, unsigned* den, int64_t mul, int64_t div)
{
    /* do nothing if the rational number is invalid */
    if (!*den)
        return;

    int64_t a;
    int64_t b;
    *num *= static_cast<unsigned>(mul);
    *den *= static_cast<unsigned>(div);
    a = *num;
    b = *den;

    while (b != 0)
    {
        int64_t t{ a };
        a = b;
        b = t % b;
    }

    if (a < 0)
        a = -a;

    *num /= static_cast<unsigned>(a);
    *den /= static_cast<unsigned>(a);
}

static AVS_VideoFrame* AVSC_CC get_frame_EEDI3CL(AVS_FilterInfo* fi, int n)
{
    EEDI3CLData* d{ static_cast<EEDI3CLData*>(fi->user_data) };

    const int field_no_prop = [&]()
    {
        if (d->field == -1)
            return (avs_get_parity(fi->child, n)) ? 1 : 0;
        else if (d->field == -2)
            return (avs_get_parity(fi->child, n >> 1)) ? 3 : 2;
        else
            return -1;
    }();

    int field{ (d->field > -1) ? d->field : field_no_prop };

    AVS_VideoFrame* src{ avs_get_frame(fi->child, (field > 1) ? (n >> 1) : n) };
    if (!src)
        return nullptr;

    AVS_VideoFrame* dst{ avs_new_video_frame_p(fi->env, &fi->vi, src) };
    AVS_VideoFrame* scp{ (d->vcheck && d->sclip) ? avs_get_frame(d->sclip, n) : nullptr };
    AVS_VideoFrame* pad[4]{};

    const int planes_y[4]{ AVS_PLANAR_Y, AVS_PLANAR_U, AVS_PLANAR_V, AVS_PLANAR_A };
    const int planes_r[4]{ AVS_PLANAR_R, AVS_PLANAR_G, AVS_PLANAR_B, AVS_PLANAR_A };
    const int* planes{ (avs_is_rgb(&fi->vi) ? planes_r : planes_y) };

    for (int plane{ 0 }; plane < avs_num_components(&fi->vi); ++plane)
    {
        if (d->process[plane])
        {
            d->vi_pad.width = avs_get_row_size_p(dst, planes[plane]) + 24 * avs_component_size(&fi->vi);
            d->vi_pad.height = avs_get_height_p(dst, planes[plane]) + 8;
            pad[plane] = avs_new_video_frame(fi->env, &d->vi_pad);
        }
    }

    if (d->field < 0)
    {
        int err;
        const int64_t field_based{ avs_prop_get_int(fi->env, avs_get_frame_props_ro(fi->env, src), "_FieldBased", 0, &err) };
        if (err == 0)
        {
            if (field_based == 1)
                field = 0;
            else if (field_based == 2)
                field = 1;

            if (d->field > 1 || field_no_prop > 1)
            {
                if (field_based == 0)
                    field -= 2;

                field = static_cast<int>((n & 1) ? (field == 0) : (field == 1));
            }
        }
        else
        {
            if (field > 1)
            {
                field -= 2;
                field = static_cast<int>((n & 1) ? (field == 0) : (field == 1));
            }
        }
    }
    else
    {
        if (field > 1)
        {
            field -= 2;
            field = static_cast<int>((n & 1) ? (field == 0) : (field == 1));
        }
    }

    try
    {
        d->filter(src, scp, dst, pad, field, d);
    }
    catch (const boost::compute::opencl_error& error)
    {
        const std::string err{ std::string("EEDI3CL: ") + error.error_string() };
        d->err = std::make_unique<char[]>(err.size() + 1);
        strcpy(d->err.get(), err.c_str());
        fi->error = d->err.get();
        avs_release_video_frame(src);
        avs_release_video_frame(scp);
        avs_release_video_frame(dst);

        for (int i{ 0 }; i < 4; ++i)
            avs_release_video_frame(pad[i]);

        return nullptr;
    }

    AVS_Map* props{ avs_get_frame_props_rw(fi->env, dst) };
    avs_prop_set_int(fi->env, props, "_FieldBased", 0, 0);

    if (d->field > 1 || field_no_prop > 1)
    {
        int errNum;
        int errDen;
        unsigned durationNum{ static_cast<unsigned>(avs_prop_get_int(fi->env, props, "_DurationNum", 0, &errNum)) };
        unsigned durationDen{ static_cast<unsigned>(avs_prop_get_int(fi->env, props, "_DurationDen", 0, &errDen)) };
        if (errNum == 0 && errDen == 0)
        {
            muldivRational(&durationNum, &durationDen, 1, 2);
            avs_prop_set_int(fi->env, props, "_DurationNum", durationNum, 0);
            avs_prop_set_int(fi->env, props, "_DurationDen", durationDen, 0);
        }
    }

    avs_release_video_frame(src);
    avs_release_video_frame(scp);

    for (int i{ 0 }; i < 4; ++i)
        avs_release_video_frame(pad[i]);

    return dst;
}

static void AVSC_CC free_EEDI3CL(AVS_FilterInfo* fi)
{
    EEDI3CLData* d{ static_cast<EEDI3CLData*>(fi->user_data) };

    avs_release_clip(d->sclip);
    aligned_free(d->pcosts);
    aligned_free(d->pbackt);
    delete[] d->fpath;
    delete[] d->dmap;

    if (d->tline)
        delete[] d->tline;

    delete d;
}

static int AVSC_CC set_cache_hints_EEDI3CL(AVS_FilterInfo* fi, int cachehints, int frame_range)
{
    return cachehints == AVS_CACHE_GET_MTMODE ? 2 : 0;
}

static AVS_Value AVSC_CC Create_EEDI3CL(AVS_ScriptEnvironment* env, AVS_Value args, void* param)
{
    enum { Clip, Field, Dh, Planes, Alpha, Beta, Gamma, Nrad, Mdis, Hp, Ucubic, Cost3, Vcheck, Vthresh0, Vthresh1, Vthresh2, Sclip, Opt, Device, List_device, Info };

    EEDI3CLData* d{ new EEDI3CLData() };

    AVS_Clip* clip{ avs_new_c_filter(env, &d->fi, avs_array_elt(args, Clip), 1) };
    AVS_Value v{ avs_void };

    d->sclip = (avs_defined(avs_array_elt(args, Sclip))) ? avs_take_clip(avs_array_elt(args, Sclip), env) : nullptr;

    try {
        if (!avs_is_planar(&d->fi->vi))
            throw std::string{ "only planar format is supported" };

        d->field = avs_defined(avs_array_elt(args, Field)) ? avs_as_int(avs_array_elt(args, Field)) : -1;
        d->dh = avs_defined(avs_array_elt(args, Dh)) ? avs_as_bool(avs_array_elt(args, Dh)) : 0;

        const int num_planes{ (avs_defined(avs_array_elt(args, Planes))) ? avs_array_size(avs_array_elt(args, Planes)) : 0 };

        for (int i{ 0 }; i < 4; ++i)
            d->process[i] = (num_planes <= 0);

        for (int i{ 0 }; i < num_planes; ++i)
        {
            const int n{ avs_as_int(*(avs_as_array(avs_array_elt(args, Planes)) + i)) };

            if (n >= avs_num_components(&d->fi->vi))
                throw std::string{ "plane index out of range" };

            if (d->process[n])
                throw std::string{ "plane specified twice" };

            d->process[n] = true;
        }

        float alpha{ avs_defined(avs_array_elt(args, Alpha)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Alpha))) : 0.2f };
        float beta{ avs_defined(avs_array_elt(args, Beta)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Beta))) : 0.25f };
        d->gamma = avs_defined(avs_array_elt(args, Gamma)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Gamma))) : 20.0f;
        const int nrad{ avs_defined(avs_array_elt(args, Nrad)) ? avs_as_int(avs_array_elt(args, Nrad)) : 2 };
        d->mdis = avs_defined(avs_array_elt(args, Mdis)) ? avs_as_int(avs_array_elt(args, Mdis)) : 20;
        d->ucubic = avs_defined(avs_array_elt(args, Ucubic)) ? avs_as_bool(avs_array_elt(args, Ucubic)) : 1;
        const int cost3{ avs_defined(avs_array_elt(args, Cost3)) ? avs_as_bool(avs_array_elt(args, Cost3)) : 1 };
        d->vcheck = avs_defined(avs_array_elt(args, Vcheck)) ? avs_as_int(avs_array_elt(args, Vcheck)) : 2;
        float vthresh0{ avs_defined(avs_array_elt(args, Vthresh0)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Vthresh0))) : 32.0f };
        float vthresh1{ avs_defined(avs_array_elt(args, Vthresh1)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Vthresh1))) : 64.0f };
        d->vthresh2 = avs_defined(avs_array_elt(args, Vthresh2)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Vthresh2))) : 4.0f;
        const int opt{ avs_defined(avs_array_elt(args, Opt)) ? avs_as_int(avs_array_elt(args, Opt)) : -1 };
        const int device_id{ avs_defined(avs_array_elt(args, Device)) ? avs_as_int(avs_array_elt(args, Device)) : -1 };

        if (d->field < -2 || d->field > 3)
            throw std::string{ "field must be -2, -1, 0, 1, 2, or 3" };
        if (!d->dh && (d->fi->vi.height & 1))
            throw std::string{ "height must be mod 2 when dh=False" };
        if (d->dh && d->field > 1)
            throw std::string{ "field must be 0 or 1 when dh=True" };
        if (alpha < 0.0f || alpha > 1.0f)
            throw std::string{ "alpha must be between 0.0 and 1.0 (inclusive)" };
        if (beta < 0.0f || beta > 1.0f)
            throw std::string{ "beta must be between 0.0 and 1.0 (inclusive)" };
        if (alpha + beta > 1.0f)
            throw std::string{ "alpha+beta must be between 0.0 and 1.0 (inclusive)" };
        if (d->gamma < 0.0f)
            throw std::string{ "gamma must be greater than or equal to 0.0" };
        if (nrad < 0 || nrad > 3)
            throw std::string{ "nrad must be between 0 and 3 (inclusive)" };
        if (d->mdis < 1 || d->mdis > 40)
            throw std::string{ "mdis must be between 1 and 40 (inclusive)" };
        if (d->vcheck < 0 || d->vcheck > 3)
            throw std::string{ "vcheck must be 0, 1, 2, or 3" };
        if (d->vcheck && (vthresh0 <= 0.0f || vthresh1 <= 0.0f || d->vthresh2 <= 0.0f))
            throw std::string{ "vthresh0, vthresh1, and vthresh2 must be greater than 0.0" };
        if (opt < -1 || opt > 1)
            throw std::string{ "opt must be -1, 0, or 1" };

        if (device_id >= static_cast<int>(boost::compute::system::device_count()))
            throw std::string{ "device index out of range" };

        if (avs_defined(avs_array_elt(args, List_device)) ? avs_as_bool(avs_array_elt(args, List_device)) : 0)
        {
            avs_release_clip(d->sclip);

            const auto devices{ boost::compute::system::devices() };
            std::string text;

            for (size_t i{ 0 }; i < devices.size(); ++i)
                text += std::to_string(i) + ": " + devices[i].name() + " (" + devices[i].platform().name() + ")" + "\n";

            d->err = std::make_unique<char[]>(text.size() + 1);
            strcpy(d->err.get(), text.c_str());

            AVS_Value cl{ avs_new_value_clip(clip) };
            AVS_Value args_[2]{ cl , avs_new_value_string(d->err.get()) };
            AVS_Value inv{ avs_invoke(d->fi->env, "Text", avs_new_value_array(args_, 2), 0) };
            AVS_Clip* clip1{ avs_take_clip(inv, env) };

            v = avs_new_value_clip(clip1);

            avs_release_clip(clip1);
            avs_release_value(inv);
            avs_release_value(cl);
            avs_release_clip(clip);

            return v;
        }

        boost::compute::device device{ boost::compute::system::default_device() };

        if (device_id > -1)
            device = boost::compute::system::devices().at(device_id);

        boost::compute::context context = boost::compute::context{ device };

        if (avs_defined(avs_array_elt(args, Info)) ? avs_as_bool(avs_array_elt(args, Info)) : 0)
        {
            avs_release_clip(d->sclip);

            std::string text{ "=== Platform Info ==\n" };
            const auto platform{ device.platform() };
            text += "Profile: " + platform.get_info<CL_PLATFORM_PROFILE>() + "\n";
            text += "Version: " + platform.get_info<CL_PLATFORM_VERSION>() + "\n";
            text += "Name: " + platform.get_info<CL_PLATFORM_NAME>() + "\n";
            text += "Vendor: " + platform.get_info<CL_PLATFORM_VENDOR>() + "\n";

            text += "\n";

            text += "=== Device Info ==\n";
            text += "Name: " + device.get_info<CL_DEVICE_NAME>() + "\n";
            text += "Vendor: " + device.get_info<CL_DEVICE_VENDOR>() + "\n";
            text += "Profile: " + device.get_info<CL_DEVICE_PROFILE>() + "\n";
            text += "Version: " + device.get_info<CL_DEVICE_VERSION>() + "\n";
            text += "Max compute units: " + std::to_string(device.get_info<CL_DEVICE_MAX_COMPUTE_UNITS>()) + "\n";
            text += "Max work-group size: " + std::to_string(device.get_info<CL_DEVICE_MAX_WORK_GROUP_SIZE>()) + "\n";
            const auto max_work_item_sizes{ device.get_info<CL_DEVICE_MAX_WORK_ITEM_SIZES>() };
            text += "Max work-item sizes: " + std::to_string(max_work_item_sizes[0]) + ", " + std::to_string(max_work_item_sizes[1]) + ", " + std::to_string(max_work_item_sizes[2]) + "\n";
            text += "2D image max width: " + std::to_string(device.get_info<CL_DEVICE_IMAGE2D_MAX_WIDTH>()) + "\n";
            text += "2D image max height: " + std::to_string(device.get_info<CL_DEVICE_IMAGE2D_MAX_HEIGHT>()) + "\n";
            text += "Image support: " + std::string{ device.get_info<CL_DEVICE_IMAGE_SUPPORT>() ? "CL_TRUE" : "CL_FALSE" } + "\n";
            const auto global_mem_cache_type{ device.get_info<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>() };
            if (global_mem_cache_type == CL_NONE)
                text += "Global memory cache type: CL_NONE\n";
            else if (global_mem_cache_type == CL_READ_ONLY_CACHE)
                text += "Global memory cache type: CL_READ_ONLY_CACHE\n";
            else if (global_mem_cache_type == CL_READ_WRITE_CACHE)
                text += "Global memory cache type: CL_READ_WRITE_CACHE\n";
            text += "Global memory cache size: " + std::to_string(device.get_info<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>() / 1024) + " KB\n";
            text += "Global memory size: " + std::to_string(device.get_info<CL_DEVICE_GLOBAL_MEM_SIZE>() / (1024 * 1024)) + " MB\n";
            text += "Max constant buffer size: " + std::to_string(device.get_info<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>() / 1024) + " KB\n";
            text += "Max constant arguments: " + std::to_string(device.get_info<CL_DEVICE_MAX_CONSTANT_ARGS>()) + "\n";
            text += "Local memory type: " + std::string{ device.get_info<CL_DEVICE_LOCAL_MEM_TYPE>() == CL_LOCAL ? "CL_LOCAL" : "CL_GLOBAL" } + "\n";
            text += "Local memory size: " + std::to_string(device.get_info<CL_DEVICE_LOCAL_MEM_SIZE>() / 1024) + " KB\n";
            text += "Available: " + std::string{ device.get_info<CL_DEVICE_AVAILABLE>() ? "CL_TRUE" : "CL_FALSE" } + "\n";
            text += "Compiler available: " + std::string{ device.get_info<CL_DEVICE_COMPILER_AVAILABLE>() ? "CL_TRUE" : "CL_FALSE" } + "\n";
            text += "OpenCL C version: " + device.get_info<CL_DEVICE_OPENCL_C_VERSION>() + "\n";
            text += "Linker available: " + std::string{ device.get_info<CL_DEVICE_LINKER_AVAILABLE>() ? "CL_TRUE" : "CL_FALSE" } + "\n";
            text += "Image max buffer size: " + std::to_string(device.get_info<size_t>(CL_DEVICE_IMAGE_MAX_BUFFER_SIZE) / 1024) + " KB" + "\n";
            text += "Out of order (on host): " + std::string{ !!(device.get_info<CL_DEVICE_QUEUE_ON_HOST_PROPERTIES>() & 1) ? "CL_TRUE" : "CL_FALSE" } + "\n";
            text += "Out of order (on device): " + std::string{ !!(device.get_info<CL_DEVICE_QUEUE_ON_DEVICE_PROPERTIES>() & 1) ? "CL_TRUE" : "CL_FALSE" };

            d->err = std::make_unique<char[]>(text.size() + 1);
            strcpy(d->err.get(), text.c_str());

            AVS_Value cl{ avs_new_value_clip(clip) };
            AVS_Value args_[2]{ cl, avs_new_value_string(d->err.get()) };
            AVS_Value inv{ avs_invoke(d->fi->env, "Text", avs_new_value_array(args_, 2), 0) };
            AVS_Clip* clip1{ avs_take_clip(inv, env) };

            v = avs_new_value_clip(clip1);

            avs_release_clip(clip1);
            avs_release_value(inv);
            avs_release_value(cl);
            avs_release_clip(clip);

            return v;
        }

        if (d->field == -2 || d->field > 1)
        {
            if (d->fi->vi.num_frames > INT_MAX / 2)
                throw std::string{ "resulting clip is too long" };

            d->fi->vi.num_frames <<= 1;

            unsigned fps_n{ d->fi->vi.fps_numerator };
            unsigned fps_d{ d->fi->vi.fps_denominator };
            muldivRational(&fps_n, &fps_d, 2, 1);
            d->fi->vi.fps_numerator = static_cast<unsigned>(fps_n);
            d->fi->vi.fps_denominator = static_cast<unsigned>(fps_d);
        }

        if (d->dh)
            d->fi->vi.height <<= 1;

        const float remainingWeight{ 1.0f - alpha - beta };

        if (cost3)
            alpha /= 3.0f;

        if (d->vcheck && d->sclip)
        {
            const AVS_VideoInfo* vi1{ avs_get_video_info(d->sclip) };

            if (!avs_is_same_colorspace(&d->fi->vi, vi1))
                throw std::string{ "sclip's format doesn't match" };
            if (vi1->num_frames != d->fi->vi.num_frames)
                throw std::string{ "sclip's number of frames doesn't match" };
            if ((vi1->width != d->fi->vi.width) || (vi1->height != d->fi->vi.height))
                throw std::string{ "sclip's dimension doesn't match" };
        }

        const int comp_size{ avs_component_size(&d->fi->vi) };
        const int iset{ instrset_detect() };
        if (opt == 1 && iset < 2)
            throw std::string{ "opt=1 requires SSE2." };

        if ((opt == -1 && iset >= 2) || opt == 1)
        {
            switch (comp_size)
            {
                case 1: d->filter = filterCL_sse2<uint8_t>; break;
                case 2: d->filter = filterCL_sse2<uint16_t>; break;
                default: d->filter = filterCL_sse2<float>; break;
            }

            d->vectorSize = 4;
        }
        else
        {
            switch (comp_size)
            {
                case 1: d->filter = filterCL_c<uint8_t>; break;
                case 2: d->filter = filterCL_c<uint16_t>; break;
                default: d->filter = filterCL_c<float>; break;
            }

            d->vectorSize = 1;
        }

        if (comp_size < 4)
        {
            d->peak = (1 << avs_bits_per_component(&d->fi->vi)) - 1;
            const float scale{ d->peak / 255.0f };
            beta *= scale;
            d->gamma *= scale;
            vthresh0 *= scale;
            vthresh1 *= scale;
        }
        else
        {
            beta /= 255.0f;
            d->gamma /= 255.0f;
            vthresh0 /= 255.0f;
            vthresh1 /= 255.0f;
        }

        d->tpitch = d->mdis * 2 + 1;
        d->tpitchVector = d->tpitch * d->vectorSize;
        d->mdisVector = d->mdis * d->vectorSize;

        d->rcpVthresh0 = 1.0f / vthresh0;
        d->rcpVthresh1 = 1.0f / vthresh1;
        d->rcpVthresh2 = 1.0f / d->vthresh2;

        switch (comp_size)
        {
            case 1: d->clImageFormat = { CL_R, CL_UNSIGNED_INT8 }; break;
            case 2: d->clImageFormat = { CL_R, CL_UNSIGNED_INT16 }; break;
            default: d->clImageFormat = { CL_R, CL_FLOAT }; break;
        }

        boost::compute::program program;

        try {
            std::ostringstream options;
            options.imbue(std::locale{ "C" });
            options.precision(16);
            options.setf(std::ios::fixed, std::ios::floatfield);
            options << "-cl-denorms-are-zero -cl-fast-relaxed-math -Werror";
            options << " -D ALPHA=" << alpha << "f";
            options << " -D BETA=" << beta << "f";
            options << " -D NRAD=" << nrad;
            options << " -D MDIS=" << d->mdis;
            options << " -D COST3=" << cost3;
            options << " -D REMAINING_WEIGHT=" << remainingWeight << "f";
            options << " -D TPITCH=" << d->tpitch;
            options << " -D VECTOR_SIZE=" << d->vectorSize;
            options << " -D MDIS_VECTOR=" << d->mdisVector;
            options << " -D TPITCH_VECTOR=" << d->tpitchVector;
            options << " -D LOCAL_WORK_SIZE_X=" << (64 / d->vectorSize);
            options << " -D LOCAL_WORK_SIZE_Y=" << d->vectorSize;
            program = boost::compute::program::build_with_source(source, context, options.str());
        }
        catch (const boost::compute::opencl_error& error)
        {
            throw error.error_string() + "\n" + program.build_log();
        }

        d->queue = boost::compute::command_queue{ context, device };

        if (comp_size < 4)
            d->calculateConnectionCosts = program.create_kernel("calculateConnectionCosts_uint");
        else
            d->calculateConnectionCosts = program.create_kernel("calculateConnectionCosts_float");

        d->src = boost::compute::image2d{ context, d->fi->vi.width + 24U, d->fi->vi.height + 8U, boost::compute::image_format{d->clImageFormat}, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY };

        d->ccosts = boost::compute::buffer{ context, d->fi->vi.width * d->tpitchVector * sizeof(cl_float), CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR | CL_MEM_HOST_READ_ONLY };

        d->pcosts = static_cast<float*>(aligned_malloc(d->fi->vi.width * d->tpitchVector * sizeof(float), 16));
        if (!d->pcosts)
            throw std::string{ "malloc failure (pcosts)" };

        d->pbackt = static_cast<int*>(aligned_malloc(d->fi->vi.width * d->tpitchVector * sizeof(int), 16));
        if (!d->pbackt)
            throw std::string{ "malloc failure (pbackt)" };

        d->fpath = new (std::nothrow) int[d->fi->vi.width * d->vectorSize];
        if (!d->fpath)
            throw std::string{ "malloc failure (fpath)" };

        d->dmap = new (std::nothrow) int[d->fi->vi.width * d->fi->vi.height];
        if (!d->dmap)
            throw std::string{ "malloc failure (dmap)" };

        if (d->vcheck)
        {
            d->tline = new (std::nothrow) int[d->fi->vi.width];
            if (!d->tline)
                throw std::string{ "malloc failure (tline)" };
        }
        else
            d->tline = nullptr;

        memcpy(&d->vi_pad, &d->fi->vi, sizeof(AVS_VideoInfo));
        d->vi_pad.pixel_type = AVS_CS_GENERIC_Y;
    }
    catch (const std::string& error)
    {
        const std::string err{ std::string("EEDI3CL: ") + error };
        d->err = std::make_unique<char[]>(err.size() + 1);
        strcpy(d->err.get(), err.c_str());
        v = avs_new_value_error(d->err.get());
    }
    catch (const boost::compute::no_device_found& error)
    {
        const std::string err{ std::string{ "EEDI3CL: " } + error.what() };
        d->err = std::make_unique<char[]>(err.size() + 1);
        strcpy(d->err.get(), err.c_str());
        v = avs_new_value_error(d->err.get());
    }
    catch (const boost::compute::opencl_error& error)
    {
        const std::string err{ std::string{ "EEDI3CL: " } + error.error_string() };
        d->err = std::make_unique<char[]>(err.size() + 1);
        strcpy(d->err.get(), err.c_str());
        v = avs_new_value_error(d->err.get());
    }

    if (!avs_defined(v))
    {
        v = avs_new_value_clip(clip);

        d->fi->user_data = reinterpret_cast<void*>(d);
        d->fi->get_frame = get_frame_EEDI3CL;
        d->fi->set_cache_hints = set_cache_hints_EEDI3CL;
        d->fi->free_filter = free_EEDI3CL;
    }

    avs_release_clip(clip);

    return v;
}

const char* AVSC_CC avisynth_c_plugin_init(AVS_ScriptEnvironment* env)
{
    avs_add_function(env, "EEDI3CL", "c[field]i[dh]b[planes]i*[alpha]f[beta]f[gamma]f[nrad]i[mdis]i[hp]b[ucubic]b[cost3]b[vcheck]i[vthresh0]f[vthresh1]f[vthresh2]f[sclip]c[opt]i[device]i[list_device]b[info]b", Create_EEDI3CL, 0);
    return "EEDI3CL";
}
