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

#include <algorithm>
#include <locale>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

#include "EEDI3CL.cl"
#include "EEDI3CL.h"
#include "EEDI3CL_kernels.cl"

template<typename T>
void copyPlane(void* dstp_, const int dstStride, const void* srcp_, const int srcStride, const int width, const int height) noexcept
{
    const T* srcp{reinterpret_cast<const T*>(srcp_)};
    T* dstp{reinterpret_cast<T*>(dstp_)};

    for (int y{0}; y < height; ++y)
    {
        std::copy_n(srcp, width, dstp);

        srcp += srcStride;
        dstp += dstStride;
    }
}

template<typename T>
void copyPad(const AVS_VideoFrame* src, AVS_VideoFrame* dst, const int plane, const int plane0, const int off, const bool dh,
    AVS_ScriptEnvironment* __restrict env) noexcept
{
    const int srcWidth{static_cast<int>(g_avs_api->avs_get_row_size_p(src, plane) / sizeof(T))};
    const int dstWidth{static_cast<int>(g_avs_api->avs_get_row_size_p(dst, plane0) / sizeof(T))};
    const int srcHeight{g_avs_api->avs_get_height_p(src, plane)};
    const int dstHeight{g_avs_api->avs_get_height_p(dst, plane0)};
    const int srcStride{static_cast<int>(g_avs_api->avs_get_pitch_p(src, plane) / sizeof(T))};
    const int dstStride{static_cast<int>(g_avs_api->avs_get_pitch_p(dst, plane0) / sizeof(T))};
    const T* srcp{reinterpret_cast<const T*>(g_avs_api->avs_get_read_ptr_p(src, plane))};
    T* dstp{reinterpret_cast<T*>(g_avs_api->avs_get_write_ptr_p(dst, plane0))};

    if (dh)
        copyPlane<T>(dstp + dstStride * (4 + off) + 12, dstStride * 2, srcp, srcStride, srcWidth, srcHeight);
    else
        copyPlane<T>(dstp + dstStride * (4 + off) + 12, dstStride * 2, srcp + srcStride * off, srcStride * 2, srcWidth, srcHeight / 2);

    dstp += dstStride * (4 + off);

    for (int y{4 + off}; y < dstHeight - 4; y += 2)
    {
        for (int x{0}; x < 12; ++x)
            dstp[x] = dstp[24 - x];

        for (int x{dstWidth - 12}, c{2}; x < dstWidth; ++x, c += 2)
            dstp[x] = dstp[x - c];

        dstp += dstStride * 2;
    }

    dstp = reinterpret_cast<T*>(g_avs_api->avs_get_write_ptr_p(dst, plane0));

    for (int y{off}; y < 4; y += 2)
        memcpy(dstp + dstStride * y, dstp + dstStride * (8 - y), dstWidth * sizeof(T));

    for (int y{dstHeight - 4 + off}, c{2 + 2 * off}; y < dstHeight; y += 2, c += 4)
        memcpy(dstp + dstStride * y, dstp + dstStride * (y - c), dstWidth * sizeof(T));
}

static void filterCL_c(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, const int field_n, bool use_dh,
    EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi)
{
    constexpr int planes_y[4]{AVS_PLANAR_Y, AVS_PLANAR_U, AVS_PLANAR_V, AVS_PLANAR_A};
    constexpr int planes_r[4]{AVS_PLANAR_R, AVS_PLANAR_G, AVS_PLANAR_B, AVS_PLANAR_A};
    const bool is_rgb{!!avs_is_rgb(&fi->vi)};
    const int* planes{(is_rgb) ? planes_r : planes_y};
    const int comp_size{g_avs_api->avs_component_size(&fi->vi)};

    for (int plane{0}; plane < g_avs_api->avs_num_components(&fi->vi); ++plane)
    {
        if (d->process[plane])
        {
            const int current_plane{planes[plane]};

            AVS_VideoInfo vi_pad;
            memcpy(&vi_pad, &fi->vi, sizeof(AVS_VideoInfo));
            vi_pad.pixel_type = AVS_CS_GENERIC_Y;
            const int src_h_for_pad{
                (use_dh) ? (g_avs_api->avs_get_height_p(src, current_plane) << 1) : g_avs_api->avs_get_height_p(src, current_plane)};
            vi_pad.width = g_avs_api->avs_get_row_size_p(src, current_plane) + 24 * g_avs_api->avs_component_size(&fi->vi);
            vi_pad.height = src_h_for_pad + 8;
            avs_helpers::avs_video_frame_ptr pad(g_avs_api->avs_new_video_frame_a(fi->env, &vi_pad, AVS_FRAME_ALIGN));
            AVS_VideoFrame* pad_raw{pad.get()};

            int peak{d->peak};

            switch (comp_size)
            {
            case 1:
                copyPad<uint8_t>(src, pad_raw, current_plane, AVS_DEFAULT_PLANE, 1 - field_n, use_dh, fi->env);
                break;
            case 2:
                copyPad<uint16_t>(src, pad_raw, current_plane, AVS_DEFAULT_PLANE, 1 - field_n, use_dh, fi->env);
                break;
            case 4:
                copyPad<float>(src, pad_raw, current_plane, AVS_DEFAULT_PLANE, 1 - field_n, use_dh, fi->env);
                peak = (plane == 0 || is_rgb) ? 1 : 0;
                break;
            default:
                break;
            }

            const int paddedWidth{static_cast<int>(g_avs_api->avs_get_row_size_p(pad_raw, AVS_DEFAULT_PLANE) / comp_size)};
            const int dstWidth{static_cast<int>(g_avs_api->avs_get_row_size_p(dst, current_plane) / comp_size)};
            const int paddedHeight{g_avs_api->avs_get_height_p(pad_raw, AVS_DEFAULT_PLANE)};
            const int dstHeight{g_avs_api->avs_get_height_p(dst, current_plane)};
            void* _dstp{g_avs_api->avs_get_write_ptr_p(dst, current_plane)};

            auto& queue{d->queue};
            auto& calculateConnectionCosts{d->calculateConnectionCosts};
            auto& srcImage{d->src};
            auto& dstImage{d->dst};
            auto& vcheckTmpImage{d->vcheck_tmp};
            auto& _ccosts{d->ccosts};
            float* pcosts{d->pcosts.get() + d->mdis};
            int* pbackt{d->pbackt.get() + d->mdis};
            int* fpath_line{d->fpath.get()};
            int* fpath_all_lines{d->dmap.get()};

            const size_t globalWorkSize[]{static_cast<size_t>((dstWidth + 63) & -64), 1};
            constexpr size_t localWorkSize[]{64, 1};
            const int bufferSize{static_cast<int>(dstWidth * d->tpitch * sizeof(cl_float))};

            queue.enqueue_write_image(srcImage, boost::compute::dim(0, 0), boost::compute::dim(paddedWidth, paddedHeight),
                g_avs_api->avs_get_read_ptr_p(pad_raw, AVS_DEFAULT_PLANE), g_avs_api->avs_get_pitch_p(pad_raw, AVS_DEFAULT_PLANE));

            for (int y{4 + field_n}; y < paddedHeight - 4; y += 2)
            {
                calculateConnectionCosts.set_args(srcImage, _ccosts, dstWidth, paddedHeight - 4, y, static_cast<cl_int>(d->ucubic), peak);
                queue.enqueue_nd_range_kernel(calculateConnectionCosts, 2, nullptr, globalWorkSize, localWorkSize);

                float* ccosts{reinterpret_cast<float*>(queue.enqueue_map_buffer(_ccosts, CL_MAP_READ, 0, bufferSize)) + d->mdis};

                *pcosts = *ccosts;
                for (int x{1}; x < dstWidth; ++x)
                {
                    const float* tT{ccosts + d->tpitch * x};
                    const float* ppT{pcosts + d->tpitch * (x - 1)};
                    float* pT{pcosts + d->tpitch * x};
                    int* piT{pbackt + d->tpitch * (x - 1)};

                    const int umax{std::min({x, dstWidth - 1 - x, d->mdis})};
                    const int umax2{std::min({x - 1, dstWidth - x, d->mdis})};
                    for (int u{-umax}; u <= umax; ++u)
                    {
                        int idx{0};
                        float bval{FLT_MAX};

                        for (int v{std::max(-umax2, u - 1)}; v <= std::min(umax2, u + 1); ++v)
                        {
                            const double z{ppT[v] + d->gamma * std::abs(u - v)};
                            const float ccost{static_cast<float>(std::min(z, FLT_MAX * 0.9))};
                            if (ccost < bval)
                            {
                                bval = ccost;
                                idx = v;
                            }
                        }

                        const double z{bval + tT[u]};
                        pT[u] = static_cast<float>(std::min(z, FLT_MAX * 0.9));
                        piT[u] = idx;
                    }
                }

                fpath_line[dstWidth - 1] = 0;

                for (int x{dstWidth - 2}; x >= 0; --x)
                    fpath_line[x] = pbackt[d->tpitch * x + fpath_line[x + 1]];

                const int line_idx{(y - 4 - field_n) / 2};
                memcpy(fpath_all_lines + static_cast<size_t>(line_idx) * dstWidth, fpath_line, dstWidth * sizeof(int));

                queue.enqueue_unmap_buffer(_ccosts, ccosts - d->mdis);
            }

            queue.enqueue_write_buffer(d->fpath_gpu, 0, static_cast<size_t>(dstWidth) * (dstHeight / 2) * sizeof(int), fpath_all_lines);

            const size_t apply_gws[2]{static_cast<size_t>((dstWidth + 15) & -16), static_cast<size_t>((dstHeight / 2 + 15) & -16)};
            const size_t apply_lws[2]{16, 16};

            d->applyInterpolationKernel.set_args(srcImage, dstImage, d->fpath_gpu, d->dmap_gpu, dstWidth, dstHeight, field_n,
                static_cast<cl_int>(d->ucubic), static_cast<cl_int>(use_dh), peak);
            queue.enqueue_nd_range_kernel(d->applyInterpolationKernel, 2, nullptr, apply_gws, apply_lws);

            boost::compute::image2d* final_image{&dstImage};

            if (d->vcheck)
            {
                boost::compute::buffer sclip_buffer;
                int sclip_stride_arg{};
                cl_int use_sclip{};

                if (d->sclip)
                {
                    use_sclip = 1;
                    const void* scpp{g_avs_api->avs_get_read_ptr_p(scp, current_plane)};
                    sclip_stride_arg = g_avs_api->avs_get_pitch_p(scp, current_plane);
                    const size_t sclip_size{static_cast<size_t>(sclip_stride_arg) * g_avs_api->avs_get_height_p(scp, current_plane)};
                    sclip_buffer = boost::compute::buffer(
                        queue.get_context(), sclip_size, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, const_cast<void*>(scpp));
                }
                else
                    sclip_buffer = boost::compute::buffer(queue.get_context(), 1, CL_MEM_READ_ONLY);

                d->vCheckKernel.set_args(srcImage, dstImage, vcheckTmpImage, d->dmap_gpu, sclip_buffer, use_sclip, dstWidth, dstHeight,
                    sclip_stride_arg, field_n, d->vcheck, d->vthresh2, d->rcpVthresh0, d->rcpVthresh1, d->rcpVthresh2, peak);
                queue.enqueue_nd_range_kernel(d->vCheckKernel, 2, nullptr, apply_gws, apply_lws);
                final_image = &vcheckTmpImage;
            }

            queue.enqueue_read_image(*final_image, boost::compute::dim(0, 0), boost::compute::dim(dstWidth, dstHeight), _dstp,
                g_avs_api->avs_get_pitch_p(dst, current_plane));
        }
    }

    d->queue.finish();
}

/* multiplies and divides a rational number, such as a frame duration, in place and reduces the result */
AVS_FORCEINLINE void muldivRational(unsigned* __restrict num, unsigned* __restrict den, int64_t mul, int64_t div) noexcept
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
        int64_t t{a};
        a = b;
        b = t % b;
    }

    if (a < 0)
        a = -a;

    *num /= static_cast<unsigned>(a);
    *den /= static_cast<unsigned>(a);
}

static void transpose_video_frame_cl(
    AVS_VideoFrame* dst, const AVS_VideoFrame* src, EEDI3CLData* __restrict d, const AVS_VideoInfo* __restrict vi)
{
    const int planes_y[4]{AVS_PLANAR_Y, AVS_PLANAR_U, AVS_PLANAR_V, AVS_PLANAR_A};
    const int planes_r[4]{AVS_PLANAR_R, AVS_PLANAR_G, AVS_PLANAR_B, AVS_PLANAR_A};
    const int* planes{(avs_is_rgb(vi) ? planes_r : planes_y)};
    auto& queue = d->queue;

    for (int i{0}; i < g_avs_api->avs_num_components(vi); ++i)
    {
        if (d->process[i])
        {
            const int plane{planes[i]};

            const int comp_size{g_avs_api->avs_component_size(vi)};
            const int src_w{g_avs_api->avs_get_row_size_p(src, plane) / comp_size};
            const int src_h{g_avs_api->avs_get_height_p(src, plane)};
            const int src_pitch{g_avs_api->avs_get_pitch_p(src, plane)};
            const int dst_pitch{g_avs_api->avs_get_pitch_p(dst, plane)};
            const void* srcp{g_avs_api->avs_get_read_ptr_p(src, plane)};
            void* dstp{g_avs_api->avs_get_write_ptr_p(dst, plane)};

            const size_t src_buffer_size{static_cast<size_t>(src_pitch) * src_h};
            const size_t dst_buffer_size{static_cast<size_t>(dst_pitch) * src_w};

            boost::compute::buffer src_buf(
                queue.get_context(), src_buffer_size, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, const_cast<void*>(srcp));
            boost::compute::buffer dst_buf(queue.get_context(), dst_buffer_size, CL_MEM_WRITE_ONLY);

            constexpr size_t TILE_DIM{16};
            const size_t global_work_size[2] = {(static_cast<size_t>(src_w) + TILE_DIM - 1) / TILE_DIM * TILE_DIM,
                (static_cast<size_t>(src_h) + TILE_DIM - 1) / TILE_DIM * TILE_DIM};
            const size_t local_work_size[2] = {TILE_DIM, TILE_DIM};

            auto run_kernel{[&](auto& kernel) {
                kernel.set_args(
                    src_buf, dst_buf, src_w, src_h, static_cast<int>(src_pitch / comp_size), static_cast<int>(dst_pitch / comp_size));
                queue.enqueue_nd_range_kernel(kernel, 2, nullptr, global_work_size, local_work_size);
            }};

            switch (comp_size)
            {
            case 1:
                run_kernel(d->transpose_u8_kernel);
                break;
            case 2:
                run_kernel(d->transpose_u16_kernel);
                break;
            default:
                run_kernel(d->transpose_f32_kernel);
                break;
            }

            queue.enqueue_read_buffer(dst_buf, 0, dst_buffer_size, dstp);
        }
    }

    queue.finish();
}

static AVS_VideoFrame* AVSC_CC get_frame_EEDI3CL(AVS_FilterInfo* __restrict fi, int n)
{
    EEDI3CLData* d{static_cast<EEDI3CLData*>(fi->user_data)};

    const int field_no_prop{[&]() {
        if (d->field == -1)
            return (g_avs_api->avs_get_parity(fi->child, n)) ? 1 : 0;
        else if (d->field == -2)
            return (g_avs_api->avs_get_parity(fi->child, n >> 1)) ? 3 : 2;
        else
            return -1;
    }()};

    int field{(d->field > -1) ? d->field : field_no_prop};

    avs_helpers::avs_video_frame_ptr src{g_avs_api->avs_get_frame(fi->child, (field > 1) ? (n >> 1) : n)};
    AVS_VideoFrame* src_raw{src.get()};

    if (!src_raw)
        return nullptr;

    const AVS_VideoInfo* src_vi{g_avs_api->avs_get_video_info(fi->child)};

    avs_helpers::avs_video_frame_ptr dst{g_avs_api->avs_new_video_frame_p(fi->env, &fi->vi, src_raw)};
    AVS_VideoFrame* dst_raw{dst.get()};
    avs_helpers::avs_video_frame_ptr scp;

    if (d->vcheck && d->sclip)
        scp.reset(g_avs_api->avs_get_frame(d->sclip.get(), n));

    AVS_VideoFrame* scp_raw{scp.get()};

    if (d->field < 0)
    {
        int err;
        const int64_t field_based{
            g_avs_api->avs_prop_get_int(fi->env, g_avs_api->avs_get_frame_props_ro(fi->env, src_raw), "_FieldBased", 0, &err)};
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
        if (!d->dw)
            d->filter(src_raw, scp_raw, dst_raw, field, d->dh, d, fi);
        else
        {
            AVS_VideoInfo vi_inter{fi->vi};
            vi_inter.width = src_vi->width;
            vi_inter.height = (d->dh) ? (src_vi->height << 1) : src_vi->height;
            avs_helpers::avs_video_frame_ptr inter{g_avs_api->avs_new_video_frame_p(fi->env, &vi_inter, nullptr)};
            AVS_VideoFrame* inter_raw{inter.get()};
            d->filter(src_raw, nullptr, inter_raw, field, d->dh, d, fi);

            AVS_VideoInfo vi_t_src{vi_inter};
            std::swap(vi_t_src.width, vi_t_src.height);
            avs_helpers::avs_video_frame_ptr src_t{g_avs_api->avs_new_video_frame_p(fi->env, &vi_t_src, nullptr)};
            AVS_VideoFrame* src_t_raw{src_t.get()};

            AVS_VideoInfo vi_t_dst{vi_t_src};
            vi_t_dst.height <<= 1;
            avs_helpers::avs_video_frame_ptr dst_t{g_avs_api->avs_new_video_frame_p(fi->env, &vi_t_dst, nullptr)};
            AVS_VideoFrame* dst_t_raw{dst_t.get()};

            transpose_video_frame_cl(src_t_raw, inter_raw, d, src_vi);
            d->filter(src_t_raw, nullptr, dst_t_raw, 0, true, d, fi);
            transpose_video_frame_cl(dst_raw, dst_t_raw, d, src_vi);
        }
    }
    catch (const boost::compute::opencl_error& error)
    {
        d->err = "EEDI3CL: " + error.error_string();
        fi->error = d->err.c_str();

        return nullptr;
    }

    AVS_Map* props{g_avs_api->avs_get_frame_props_rw(fi->env, dst_raw)};
    g_avs_api->avs_prop_set_int(fi->env, props, "_FieldBased", 0, 0);

    if (d->dw)
    {
        int errNum;
        int errDen;
        unsigned sarNum{static_cast<unsigned>(g_avs_api->avs_prop_get_int(fi->env, props, "_SARNum", 0, &errNum))};
        unsigned sarDen{static_cast<unsigned>(g_avs_api->avs_prop_get_int(fi->env, props, "_SARDen", 0, &errDen))};
        if (errNum == 0 && errDen == 0)
        {
            muldivRational(&sarNum, &sarDen, 1, 2);
            g_avs_api->avs_prop_set_int(fi->env, props, "_SARNum", sarNum, 0);
            g_avs_api->avs_prop_set_int(fi->env, props, "_SARDen", sarDen, 0);
        }
    }

    if (d->field > 1 || field_no_prop > 1)
    {
        int errNum;
        int errDen;
        unsigned durationNum{static_cast<unsigned>(g_avs_api->avs_prop_get_int(fi->env, props, "_DurationNum", 0, &errNum))};
        unsigned durationDen{static_cast<unsigned>(g_avs_api->avs_prop_get_int(fi->env, props, "_DurationDen", 0, &errDen))};
        if (errNum == 0 && errDen == 0)
        {
            muldivRational(&durationNum, &durationDen, 1, 2);
            g_avs_api->avs_prop_set_int(fi->env, props, "_DurationNum", durationNum, 0);
            g_avs_api->avs_prop_set_int(fi->env, props, "_DurationDen", durationDen, 0);
        }
    }

    return dst.release();
}

static void AVSC_CC free_EEDI3CL(AVS_FilterInfo* __restrict fi) noexcept
{
    delete static_cast<EEDI3CLData*>(fi->user_data);
    fi->user_data = nullptr;
}

static int AVSC_CC set_cache_hints_EEDI3CL(AVS_FilterInfo* __restrict fi, int cachehints, int frame_range) noexcept
{
    return cachehints == AVS_CACHE_GET_MTMODE ? 2 : 0;
}

static AVS_Value AVSC_CC Create_EEDI3CL(AVS_ScriptEnvironment* __restrict env, AVS_Value args, void* param)
{
    enum
    {
        Clip,
        Field,
        Dh,
        Dw,
        Planes,
        Alpha,
        Beta,
        Gamma,
        Nrad,
        Mdis,
        Hp,
        Ucubic,
        Cost3,
        Vcheck,
        Vthresh0,
        Vthresh1,
        Vthresh2,
        Sclip,
        Opt,
        Device,
        List_device,
        Info,
        Luma
    };

    std::unique_ptr<EEDI3CLData> d{std::make_unique<EEDI3CLData>()};

    AVS_FilterInfo* fi{};
    avs_helpers::avs_clip_ptr clip{g_avs_api->avs_new_c_filter(env, &fi, avs_array_elt(args, Clip), 1)};
    AVS_Clip* clip_raw{clip.get()};
    AVS_Value v{avs_void};

    d->sclip = avs_helpers::get_opt_arg<avs_helpers::avs_clip_ptr>(env, args, Sclip).value_or(nullptr);
    AVS_Clip* sclip_raw{d->sclip.get()};

    try
    {
        if (!avs_is_planar(&fi->vi))
            throw std::string{"only planar format is supported"};

        d->field = avs_defined(avs_array_elt(args, Field)) ? avs_as_int(avs_array_elt(args, Field)) : -1;
        d->dh = avs_defined(avs_array_elt(args, Dh)) ? avs_as_bool(avs_array_elt(args, Dh)) : 0;
        d->dw = avs_defined(avs_array_elt(args, Dw)) ? avs_as_bool(avs_array_elt(args, Dw)) : 0;

        avs_helpers::converted_array<int> num_planes_load{avs_helpers::get_opt_array_as_unique_ptr<int>(args, Planes)};
        const int num_planes{num_planes_load.size};
        const int onlyY{(avs_defined(avs_array_elt(args, Luma))) ? avs_as_bool(avs_array_elt(args, Luma)) : 0};

        if (onlyY)
        {
            d->process[0] = true;

            for (int i{1}; i < 4; ++i)
                d->process[i] = false;
        }
        else
        {
            for (int i{0}; i < 4; ++i)
                d->process[i] = (num_planes <= 0);
        }

        for (int i{0}; i < num_planes; ++i)
        {
            const int n{num_planes_load.data[i]};

            if (n >= g_avs_api->avs_num_components(&fi->vi))
                throw std::string{"plane index out of range"};

            if (d->process[n])
                throw std::string{"plane specified twice"};

            d->process[n] = true;
        }

        const int bit_depth{g_avs_api->avs_bits_per_component(&fi->vi)};

        if (onlyY && !avs_is_rgb(&fi->vi))
        {
            if (num_planes > 1)
                throw std::string{"luma cannot be true when processed planes are more than 1"};

            if (!d->process[0])
                throw std::string{"planes=0 must be used for luma=true"};

            switch (bit_depth)
            {
            case 8:
                fi->vi.pixel_type = AVS_CS_GENERIC_Y | AVS_CS_SAMPLE_BITS_8;
                break;
            case 10:
                fi->vi.pixel_type = AVS_CS_GENERIC_Y | AVS_CS_SAMPLE_BITS_10;
                break;
            case 12:
                fi->vi.pixel_type = AVS_CS_GENERIC_Y | AVS_CS_SAMPLE_BITS_12;
                break;
            case 14:
                fi->vi.pixel_type = AVS_CS_GENERIC_Y | AVS_CS_SAMPLE_BITS_14;
                break;
            case 16:
                fi->vi.pixel_type = AVS_CS_GENERIC_Y | AVS_CS_SAMPLE_BITS_16;
                break;
            case 32:
                fi->vi.pixel_type = AVS_CS_GENERIC_Y | AVS_CS_SAMPLE_BITS_32;
                break;
            default:
                break;
            }
        }

        float alpha{avs_defined(avs_array_elt(args, Alpha)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Alpha))) : 0.2f};
        float beta{avs_defined(avs_array_elt(args, Beta)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Beta))) : 0.25f};
        d->gamma = avs_defined(avs_array_elt(args, Gamma)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Gamma))) : 20.0f;
        const int nrad{avs_defined(avs_array_elt(args, Nrad)) ? avs_as_int(avs_array_elt(args, Nrad)) : 2};
        d->mdis = avs_defined(avs_array_elt(args, Mdis)) ? avs_as_int(avs_array_elt(args, Mdis)) : 20;
        d->ucubic = avs_defined(avs_array_elt(args, Ucubic)) ? avs_as_bool(avs_array_elt(args, Ucubic)) : 1;
        const int cost3{avs_defined(avs_array_elt(args, Cost3)) ? avs_as_bool(avs_array_elt(args, Cost3)) : 1};
        d->vcheck = avs_defined(avs_array_elt(args, Vcheck)) ? avs_as_int(avs_array_elt(args, Vcheck)) : 2;
        float vthresh0{
            avs_defined(avs_array_elt(args, Vthresh0)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Vthresh0))) : 32.0f};
        float vthresh1{
            avs_defined(avs_array_elt(args, Vthresh1)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Vthresh1))) : 64.0f};
        d->vthresh2 = avs_defined(avs_array_elt(args, Vthresh2)) ? static_cast<float>(avs_as_float(avs_array_elt(args, Vthresh2))) : 4.0f;
        const int opt{avs_defined(avs_array_elt(args, Opt)) ? avs_as_int(avs_array_elt(args, Opt)) : -1};
        const int device_id{avs_defined(avs_array_elt(args, Device)) ? avs_as_int(avs_array_elt(args, Device)) : -1};

        if (d->field < -2 || d->field > 3)
            throw std::string{"field must be -2, -1, 0, 1, 2, or 3"};

        if (!d->dh && (fi->vi.height & 1))
            throw std::string{"height must be mod 2 when dh=False"};

        if (d->dh && d->field > 1)
            throw std::string{"field must be 0 or 1 when dh=True"};

        if (d->dw && d->field > 1)
            throw std::string{"dw=true is not supported with field=2 or 3"};

        if (d->dw && sclip_raw)
            throw std::string{"sclip is not supported when dw=true"};

        if (alpha < 0.0f || alpha > 1.0f)
            throw std::string{"alpha must be between 0.0 and 1.0 (inclusive)"};

        if (beta < 0.0f || beta > 1.0f)
            throw std::string{"beta must be between 0.0 and 1.0 (inclusive)"};

        if (alpha + beta > 1.0f)
            throw std::string{"alpha+beta must be between 0.0 and 1.0 (inclusive)"};

        if (d->gamma < 0.0f)
            throw std::string{"gamma must be greater than or equal to 0.0"};

        if (nrad < 0 || nrad > 3)
            throw std::string{"nrad must be between 0 and 3 (inclusive)"};

        if (d->mdis < 1 || d->mdis > 40)
            throw std::string{"mdis must be between 1 and 40 (inclusive)"};

        if (d->vcheck < 0 || d->vcheck > 3)
            throw std::string{"vcheck must be 0, 1, 2, or 3"};

        if (d->vcheck && (vthresh0 <= 0.0f || vthresh1 <= 0.0f || d->vthresh2 <= 0.0f))
            throw std::string{"vthresh0, vthresh1, and vthresh2 must be greater than 0.0"};

        if (opt < -1 || opt > 3)
            throw std::string{"opt must be betwen -1 and 3 (inclusive)"};

        const int cpu_instrucs{!(!!(g_avs_api->avs_get_cpu_flags(fi->env) & AVS_CPUF_AVX512F) && (opt < 0 || opt == 3))
                                   ? !(!!(g_avs_api->avs_get_cpu_flags(fi->env) & AVS_CPUF_AVX2) && (opt < 0 || opt == 2))
                                         ? !(!!(g_avs_api->avs_get_cpu_flags(fi->env) & AVS_CPU_SSE2) && (opt < 0 || opt == 1)) ? 0 : 1
                                         : 2
                                   : 3};

        if (opt == 3 && cpu_instrucs != 3)
            throw std::string{"opt=3 requires AVX-512"};
        else if (opt == 2 && cpu_instrucs != 2)
            throw std::string{"opt=2 requires AVX2"};
        else if (opt == 1 && cpu_instrucs != 1)
            throw std::string{"opt=1 requires SSE2"};

        if (device_id >= static_cast<int>(boost::compute::system::device_count()))
            throw std::string{"device index out of range"};

        if (avs_defined(avs_array_elt(args, List_device)) ? avs_as_bool(avs_array_elt(args, List_device)) : 0)
        {
            const auto devices{boost::compute::system::devices()};

            for (size_t i{0}; i < devices.size(); ++i)
                d->err += std::to_string(i) + ": " + devices[i].name() + " (" + devices[i].platform().name() + ")" + "\n";

            AVS_Value cl;
            g_avs_api->avs_set_to_clip(&cl, clip_raw);
            avs_helpers::avs_value_guard cl_guard(cl);
            AVS_Value args_[2]{cl_guard.get(), avs_new_value_string(d->err.c_str())};
            v = g_avs_api->avs_invoke(fi->env, "Text", avs_new_value_array(args_, 2), 0);

            fi->user_data = d.release();
            fi->free_filter = free_EEDI3CL;

            return v;
        }

        boost::compute::device device{boost::compute::system::default_device()};

        if (device_id > -1)
            device = boost::compute::system::devices().at(device_id);

        boost::compute::context context{boost::compute::context{device}};

        if (avs_defined(avs_array_elt(args, Info)) ? avs_as_bool(avs_array_elt(args, Info)) : 0)
        {
            d->err = "=== Platform Info ==\n";
            const auto platform{device.platform()};
            d->err += "Profile: " + platform.get_info<CL_PLATFORM_PROFILE>() + "\n";
            d->err += "Version: " + platform.get_info<CL_PLATFORM_VERSION>() + "\n";
            d->err += "Name: " + platform.get_info<CL_PLATFORM_NAME>() + "\n";
            d->err += "Vendor: " + platform.get_info<CL_PLATFORM_VENDOR>() + "\n";

            d->err += "\n";

            d->err += "=== Device Info ==\n";
            d->err += "Name: " + device.get_info<CL_DEVICE_NAME>() + "\n";
            d->err += "Vendor: " + device.get_info<CL_DEVICE_VENDOR>() + "\n";
            d->err += "Profile: " + device.get_info<CL_DEVICE_PROFILE>() + "\n";
            d->err += "Version: " + device.get_info<CL_DEVICE_VERSION>() + "\n";
            d->err += "Max compute units: " + std::to_string(device.get_info<CL_DEVICE_MAX_COMPUTE_UNITS>()) + "\n";
            d->err += "Max work-group size: " + std::to_string(device.get_info<CL_DEVICE_MAX_WORK_GROUP_SIZE>()) + "\n";
            const auto max_work_item_sizes{device.get_info<CL_DEVICE_MAX_WORK_ITEM_SIZES>()};
            d->err += "Max work-item sizes: " + std::to_string(max_work_item_sizes[0]) + ", " + std::to_string(max_work_item_sizes[1]) +
                      ", " + std::to_string(max_work_item_sizes[2]) + "\n";
            d->err += "2D image max width: " + std::to_string(device.get_info<CL_DEVICE_IMAGE2D_MAX_WIDTH>()) + "\n";
            d->err += "2D image max height: " + std::to_string(device.get_info<CL_DEVICE_IMAGE2D_MAX_HEIGHT>()) + "\n";
            d->err += "Image support: " + std::string{device.get_info<CL_DEVICE_IMAGE_SUPPORT>() ? "CL_TRUE" : "CL_FALSE"} + "\n";
            const auto global_mem_cache_type{device.get_info<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>()};

            if (global_mem_cache_type == CL_NONE)
                d->err += "Global memory cache type: CL_NONE\n";
            else if (global_mem_cache_type == CL_READ_ONLY_CACHE)
                d->err += "Global memory cache type: CL_READ_ONLY_CACHE\n";
            else if (global_mem_cache_type == CL_READ_WRITE_CACHE)
                d->err += "Global memory cache type: CL_READ_WRITE_CACHE\n";

            d->err += "Global memory cache size: " + std::to_string(device.get_info<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>() / 1024) + " KB\n";
            d->err += "Global memory size: " + std::to_string(device.get_info<CL_DEVICE_GLOBAL_MEM_SIZE>() / (1024 * 1024)) + " MB\n";
            d->err += "Max constant buffer size: " + std::to_string(device.get_info<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>() / 1024) + " KB\n";
            d->err += "Max constant arguments: " + std::to_string(device.get_info<CL_DEVICE_MAX_CONSTANT_ARGS>()) + "\n";
            d->err +=
                "Local memory type: " + std::string{device.get_info<CL_DEVICE_LOCAL_MEM_TYPE>() == CL_LOCAL ? "CL_LOCAL" : "CL_GLOBAL"} +
                "\n";
            d->err += "Local memory size: " + std::to_string(device.get_info<CL_DEVICE_LOCAL_MEM_SIZE>() / 1024) + " KB\n";
            d->err += "Available: " + std::string{device.get_info<CL_DEVICE_AVAILABLE>() ? "CL_TRUE" : "CL_FALSE"} + "\n";
            d->err += "Compiler available: " + std::string{device.get_info<CL_DEVICE_COMPILER_AVAILABLE>() ? "CL_TRUE" : "CL_FALSE"} + "\n";
            d->err += "OpenCL C version: " + device.get_info<CL_DEVICE_OPENCL_C_VERSION>() + "\n";
            d->err += "Linker available: " + std::string{device.get_info<CL_DEVICE_LINKER_AVAILABLE>() ? "CL_TRUE" : "CL_FALSE"} + "\n";
            d->err +=
                "Image max buffer size: " + std::to_string(device.get_info<size_t>(CL_DEVICE_IMAGE_MAX_BUFFER_SIZE) / 1024) + " KB" + "\n";
            d->err += "Out of order (on host): " +
                      std::string{!!(device.get_info<CL_DEVICE_QUEUE_ON_HOST_PROPERTIES>() & 1) ? "CL_TRUE" : "CL_FALSE"} + "\n";
            d->err += "Out of order (on device): " +
                      std::string{!!(device.get_info<CL_DEVICE_QUEUE_ON_DEVICE_PROPERTIES>() & 1) ? "CL_TRUE" : "CL_FALSE"};

            AVS_Value cl;
            g_avs_api->avs_set_to_clip(&cl, clip_raw);
            avs_helpers::avs_value_guard cl_guard(cl);
            AVS_Value args_[2]{cl_guard.get(), avs_new_value_string(d->err.c_str())};
            v = g_avs_api->avs_invoke(fi->env, "Text", avs_new_value_array(args_, 2), 0);

            fi->user_data = d.release();
            fi->free_filter = free_EEDI3CL;

            return v;
        }

        if (d->field == -2 || d->field > 1)
        {
            if (fi->vi.num_frames > INT_MAX / 2)
                throw std::string{"resulting clip is too long"};

            fi->vi.num_frames <<= 1;

            unsigned fps_n{fi->vi.fps_numerator};
            unsigned fps_d{fi->vi.fps_denominator};
            muldivRational(&fps_n, &fps_d, 2, 1);
            fi->vi.fps_numerator = static_cast<unsigned>(fps_n);
            fi->vi.fps_denominator = static_cast<unsigned>(fps_d);
        }

        if (d->dw)
        {
            if (fi->vi.width > INT_MAX / 2)
                throw std::string{"resulting clip is too wide"};
            fi->vi.width <<= 1;
        }

        if (d->dh)
            fi->vi.height <<= 1;

        const float remainingWeight{1.0f - alpha - beta};

        if (cost3)
            alpha /= 3.0f;

        if (d->vcheck && sclip_raw)
        {
            const AVS_VideoInfo* vi1{g_avs_api->avs_get_video_info(sclip_raw)};

            if (!((fi->vi.pixel_type == vi1->pixel_type) || (g_avs_api->avs_is_yv12(&fi->vi) && g_avs_api->avs_is_yv12(vi1))))
                throw std::string{"sclip's format doesn't match"};
            if (vi1->num_frames != fi->vi.num_frames)
                throw std::string{"sclip's number of frames doesn't match"};
            if ((vi1->width != fi->vi.width) || (vi1->height != fi->vi.height))
                throw std::string{"sclip's dimension doesn't match"};
        }

        const AVS_VideoInfo* child_vi{g_avs_api->avs_get_video_info(clip_raw)};
        int max_processing_width{child_vi->width};
        int max_processing_height{(d->dh) ? (child_vi->height << 1) : child_vi->height};

        if (d->dw)
        {
            const int pass2_w{(d->dh) ? (child_vi->height << 1) : child_vi->height};
            const int pass2_h{child_vi->width << 1};
            max_processing_width = std::max(max_processing_width, pass2_w);
            max_processing_height = std::max(max_processing_height, pass2_h);
        }

        const int dmap_elements{max_processing_width * max_processing_height};
        const int comp_size{g_avs_api->avs_component_size(&fi->vi)};
        int local_work_size{64};

        switch (cpu_instrucs)
        {
        case 3:
            d->vectorSize = 16;
            local_work_size = 4;

            switch (comp_size)
            {
            case 1:
                d->filter = filterCL_avx512<uint8_t>;
                break;
            case 2:
                d->filter = filterCL_avx512<uint16_t>;
                break;
            default:
                d->filter = filterCL_avx512<float>;
                break;
            }
            break;
        case 2:
            d->vectorSize = 8;
            local_work_size = 8;

            switch (comp_size)
            {
            case 1:
                d->filter = filterCL_avx2<uint8_t>;
                break;
            case 2:
                d->filter = filterCL_avx2<uint16_t>;
                break;
            default:
                d->filter = filterCL_avx2<float>;
                break;
            }
            break;
        case 1:
            d->vectorSize = 4;
            local_work_size = 16;

            switch (comp_size)
            {
            case 1:
                d->filter = filterCL_sse2<uint8_t>;
                break;
            case 2:
                d->filter = filterCL_sse2<uint16_t>;
                break;
            default:
                d->filter = filterCL_sse2<float>;
                break;
            }
            break;
        case 0:
            d->vectorSize = 1;
            d->filter = filterCL_c;
            break;
        default:
            break;
        }

        if (comp_size < 4)
        {
            d->peak = (1 << bit_depth) - 1;
            const float scale{d->peak / 255.0f};
            beta *= scale;
            d->gamma *= scale;
            vthresh0 *= scale;
            vthresh1 *= scale;
        }
        else
        {
            d->peak = 0;
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
        case 1:
            d->clImageFormat = {CL_R, CL_UNSIGNED_INT8};
            break;
        case 2:
            d->clImageFormat = {CL_R, CL_UNSIGNED_INT16};
            break;
        default:
            d->clImageFormat = {CL_R, CL_FLOAT};
            break;
        }

        boost::compute::program program;
        const std::string final_source_code{std::string(source) + std::string(interpolation_kernels_source)};

        try
        {
            std::ostringstream options;
            options.imbue(std::locale{"C"});
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
            options << " -D LOCAL_WORK_SIZE_X=" << local_work_size;
            options << " -D LOCAL_WORK_SIZE_Y=" << d->vectorSize;
            program = boost::compute::program::build_with_source(final_source_code, context, options.str());
        }
        catch (const boost::compute::opencl_error& error)
        {
            throw error.error_string() + "\n" + program.build_log();
        }

        d->queue = boost::compute::command_queue{context, device};

        if (comp_size < 4)
        {
            d->calculateConnectionCosts = program.create_kernel("calculateConnectionCosts_uint");
            d->applyInterpolationKernel = program.create_kernel("applyInterpolation_uint");
            d->vCheckKernel = program.create_kernel("vCheck_uint");
        }
        else
        {
            d->calculateConnectionCosts = program.create_kernel("calculateConnectionCosts_float");
            d->applyInterpolationKernel = program.create_kernel("applyInterpolation_float");
            d->vCheckKernel = program.create_kernel("vCheck_float");
        }

        d->transpose_u8_kernel = program.create_kernel("transpose_u8");
        d->transpose_u16_kernel = program.create_kernel("transpose_u16");
        d->transpose_f32_kernel = program.create_kernel("transpose_f32");

        d->src = boost::compute::image2d{context, static_cast<size_t>(max_processing_width) + 24U,
            static_cast<size_t>(max_processing_height) + 8U, boost::compute::image_format{d->clImageFormat},
            CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY};
        d->dst = boost::compute::image2d{context, static_cast<size_t>(max_processing_width), static_cast<size_t>(max_processing_height),
            boost::compute::image_format{d->clImageFormat}, CL_MEM_READ_WRITE};
        d->vcheck_tmp = boost::compute::image2d{context, static_cast<size_t>(max_processing_width),
            static_cast<size_t>(max_processing_height), boost::compute::image_format{d->clImageFormat}, CL_MEM_WRITE_ONLY};

        const size_t fpath_size{static_cast<size_t>(max_processing_width) * (max_processing_height >> 1) * sizeof(int)};
        d->fpath_gpu = boost::compute::buffer{context, fpath_size, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY};

        const size_t dmap_gpu_size{static_cast<size_t>(dmap_elements) * sizeof(int)};
        d->dmap_gpu = boost::compute::buffer{context, dmap_gpu_size, CL_MEM_READ_WRITE};

        d->ccosts = boost::compute::buffer{context, static_cast<size_t>(max_processing_width) * d->tpitchVector * sizeof(cl_float),
            CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR | CL_MEM_HOST_READ_ONLY};

        d->pcosts = make_unique_aligned_array<float>(static_cast<size_t>(max_processing_width) * d->tpitchVector * sizeof(float), 64);

        if (!d->pcosts)
            throw std::string{"malloc failure (pcosts)"};

        d->pbackt = make_unique_aligned_array<int>(static_cast<size_t>(max_processing_width) * d->tpitchVector * sizeof(int), 64);
        if (!d->pbackt)
            throw std::string{"malloc failure (pbackt)"};

        d->fpath = std::make_unique<int[]>(static_cast<size_t>(max_processing_width) * d->vectorSize);

        if (!d->fpath)
            throw std::string{"malloc failure (fpath)"};

        d->dmap = std::make_unique<int[]>(dmap_elements);

        if (!d->dmap)
            throw std::string{"malloc failure (dmap)"};
    }
    catch (const std::string& error)
    {
        d->err = "EEDI3CL: " + error;
        v = avs_new_value_error(d->err.c_str());
    }
    catch (const boost::compute::no_device_found& error)
    {
        d->err = std::string{"EEDI3CL: "} + error.what();
        v = avs_new_value_error(d->err.c_str());
    }
    catch (const boost::compute::opencl_error& error)
    {
        d->err = "EEDI3CL: " + error.error_string();
        v = avs_new_value_error(d->err.c_str());
    }

    if (!avs_defined(v))
    {
        g_avs_api->avs_set_to_clip(&v, clip_raw);

        fi->user_data = d.release();
        fi->get_frame = get_frame_EEDI3CL;
        fi->set_cache_hints = set_cache_hints_EEDI3CL;
        fi->free_filter = free_EEDI3CL;
    }

    return v;
}

const char* AVSC_CC avisynth_c_plugin_init(AVS_ScriptEnvironment* __restrict env)
{
    static constexpr int REQUIRED_INTERFACE_VERSION{9};
    static constexpr int REQUIRED_BUGFIX_VERSION{2};
    static constexpr std::initializer_list<std::string_view> required_functions{"avs_pool_free", // avs loader helper functions
        "avs_release_clip",                                                                      // avs loader helper functions
        "avs_release_value",                                                                     // avs loader helper functions
        "avs_release_video_frame",                                                               // avs loader helper functions
        "avs_take_clip",                                                                         // avs loader helper functions
        "avs_add_function",
        "avs_new_c_filter",
        "avs_new_video_frame_p",
        "avs_new_video_frame_a",
        "avs_invoke",
        "avs_set_to_clip",
        "avs_get_frame",
        "avs_get_height_p",
        "avs_get_row_size_p",
        "avs_get_pitch_p",
        "avs_get_read_ptr_p",
        "avs_get_video_info",
        "avs_get_write_ptr_p",
        "avs_get_parity",
        "avs_get_cpu_flags",
        "avs_prop_set_int",
        "avs_prop_get_int",
        "avs_bits_per_component",
        "avs_component_size",
        "avs_num_components"};

    if (!avisynth_c_api_loader::get_api(env, REQUIRED_INTERFACE_VERSION, REQUIRED_BUGFIX_VERSION, required_functions))
    {
        std::cerr << avisynth_c_api_loader::get_last_error() << std::endl;
        return avisynth_c_api_loader::get_last_error();
    }

    g_avs_api->avs_add_function(env, "EEDI3CL",
        "c"
        "[field]i"
        "[dh]b"
        "[dw]b"
        "[planes]i*"
        "[alpha]f"
        "[beta]f"
        "[gamma]f"
        "[nrad]i"
        "[mdis]i"
        "[hp]b"
        "[ucubic]b"
        "[cost3]b"
        "[vcheck]i"
        "[vthresh0]f"
        "[vthresh1]f"
        "[vthresh2]f"
        "[sclip]c"
        "[opt]i"
        "[device]i"
        "[list_device]b"
        "[info]b"
        "[luma]b",
        Create_EEDI3CL, 0);
    return "EEDI3CL";
}
