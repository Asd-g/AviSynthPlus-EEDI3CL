#include "EEDI3CL.h"
#include "VCL2/vectorclass.h"

using V_float = Vec16f;
using V_int = std::conditional_t <std::is_same_v<V_float, Vec4f>, Vec4i, std::conditional_t<std::is_same_v<V_float, Vec8f>, Vec8i, Vec16i>>;
using V_ibool = std::conditional_t <std::is_same_v<V_float, Vec4f>, Vec4ib, std::conditional_t<std::is_same_v<V_float, Vec8f>, Vec8ib, Vec16ib>>;

template<typename T>
void filterCL_avx512(const AVS_VideoFrame* __restrict src, const AVS_VideoFrame* __restrict scp, AVS_VideoFrame* __restrict dst, const int field_n, bool use_dh, EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi)
{
    constexpr int planes_y[4]{ AVS_PLANAR_Y, AVS_PLANAR_U, AVS_PLANAR_V, AVS_PLANAR_A };
    constexpr int planes_r[4]{ AVS_PLANAR_R, AVS_PLANAR_G, AVS_PLANAR_B, AVS_PLANAR_A };
    const bool is_rgb{ !!avs_is_rgb(&fi->vi) };
    const int* planes{ (is_rgb) ? planes_r : planes_y };

    for (int plane{ 0 }; plane < avs_num_components(&fi->vi); ++plane)
    {
        if (d->process[plane])
        {
            const int current_plane{ planes[plane] };

            AVS_VideoInfo vi_pad{};
            memcpy(&vi_pad, &fi->vi, sizeof(AVS_VideoInfo));
            vi_pad.pixel_type = AVS_CS_GENERIC_Y;
            const int src_h_for_pad{ (use_dh) ? (avs_get_height_p(src, current_plane) << 1) : avs_get_height_p(src, current_plane) };
            vi_pad.width = avs_get_row_size_p(src, current_plane) + 24 * avs_component_size(&fi->vi);
            vi_pad.height = src_h_for_pad + 8;
            AVS_VideoFrame* pad{ avs_new_video_frame(fi->env, &vi_pad) };

            int peak{ d->peak };

            if constexpr (std::is_same_v<T, uint8_t>)
                copyPad<uint8_t>(src, pad, current_plane, AVS_DEFAULT_PLANE, 1 - field_n, use_dh, fi->env);
            else if constexpr (std::is_same_v<T, uint16_t>)
                copyPad<uint16_t>(src, pad, current_plane, AVS_DEFAULT_PLANE, 1 - field_n, use_dh, fi->env);
            else
            {
                copyPad<float>(src, pad, current_plane, AVS_DEFAULT_PLANE, 1 - field_n, use_dh, fi->env);
                peak = (plane == 0 || is_rgb) ? 1 : 0;
            }

            const int paddedWidth{ avs_get_row_size_p(pad, AVS_DEFAULT_PLANE) / static_cast<int>(sizeof(T)) };
            const int dstWidth{ avs_get_row_size_p(dst, current_plane) / static_cast<int>(sizeof(T)) };
            const int paddedHeight{ avs_get_height_p(pad, AVS_DEFAULT_PLANE) };
            const int dstHeight{ avs_get_height_p(dst, current_plane) };
            void* __restrict _dstp{ avs_get_write_ptr_p(dst, current_plane) };

            auto& queue{ d->queue };
            auto& calculateConnectionCosts{ d->calculateConnectionCosts };
            auto& srcImage{ d->src };
            auto& dstImage{ d->dst };
            auto& vcheckTmpImage{ d->vcheck_tmp };
            auto& _ccosts{ d->ccosts };
            float* pcosts{ d->pcosts + d->mdisVector };
            int* pbackt{ d->pbackt + d->mdisVector };
            int* fpath_line{ d->fpath };
            int* fpath_all_lines{ d->dmap };

            const size_t globalWorkSize[]{ static_cast<size_t>((dstWidth + 3) & -4), static_cast<size_t>(d->vectorSize) };
            constexpr size_t localWorkSize[]{ 4, 16 };
            const int bufferSize{ static_cast<int>(dstWidth * d->tpitchVector * sizeof(cl_float)) };

            queue.enqueue_write_image(srcImage, boost::compute::dim(0, 0), boost::compute::dim(paddedWidth, paddedHeight),
                avs_get_read_ptr_p(pad, AVS_DEFAULT_PLANE), avs_get_pitch_p(pad, AVS_DEFAULT_PLANE));

            for (int y{ 4 + field_n }; y < paddedHeight - 4; y += 2 * d->vectorSize)
            {
                calculateConnectionCosts.set_args(srcImage, _ccosts, dstWidth, paddedHeight - 4, y, static_cast<cl_int>(d->ucubic), peak);
                queue.enqueue_nd_range_kernel(calculateConnectionCosts, 2, nullptr, globalWorkSize, localWorkSize);

                float* ccosts{ reinterpret_cast<float*>(queue.enqueue_map_buffer(_ccosts, CL_MAP_READ, 0, bufferSize)) + d->mdisVector };

                V_float().load(ccosts).store_a(pcosts);
                for (int x{ 1 }; x < dstWidth; ++x)
                {
                    const float* tT{ ccosts + static_cast<size_t>(d->tpitchVector) * x };
                    const float* ppT{ pcosts + static_cast<size_t>(d->tpitchVector) * (x - 1) };
                    float* pT{ pcosts + static_cast<size_t>(d->tpitchVector) * x };
                    int* piT{ pbackt + static_cast<size_t>(d->tpitchVector) * (x - 1) };

                    const int umax{ std::min({ x, dstWidth - 1 - x, d->mdis }) };
                    const int umax2{ std::min({ x - 1, dstWidth - x, d->mdis }) };
                    for (int u{ -umax }; u <= umax; ++u)
                    {
                        V_int idx{ zero_si512() };
                        V_float bval{ FLT_MAX };

                        for (int v{ std::max(-umax2, u - 1) }; v <= std::min(umax2, u + 1); ++v)
                        {
                            const auto z{ V_float().load_a(ppT + static_cast<size_t>(v) * d->vectorSize) + d->gamma * std::abs(u - v) };
                            const auto ccost{ min(z, FLT_MAX * 0.9f) };
                            idx = select(V_ibool(ccost < bval), v, idx);
                            bval = min(ccost, bval);
                        }

                        const auto z{ bval + V_float().load(tT + static_cast<size_t>(u) * d->vectorSize) };
                        min(z, FLT_MAX * 0.9f).store_a(pT + static_cast<size_t>(u) * d->vectorSize);
                        idx.store_nt(piT + static_cast<size_t>(u) * d->vectorSize);
                    }
                }

                for (int vs{ 0 }; vs < d->vectorSize; ++vs)
                {
                    const int realY{ y + 2 * vs };
                    if (realY >= paddedHeight - 4) break;

                    const int* pbackt_line{ pbackt + vs };
                    fpath_line[dstWidth - 1] = 0;
                    for (int x{ dstWidth - 2 }; x >= 0; --x) {
                        fpath_line[x] = pbackt_line[(static_cast<size_t>(d->tpitch) * x + fpath_line[x + 1]) * d->vectorSize];
                    }

                    const int line_idx{ (realY - 4 - field_n) / 2 };
                    memcpy(fpath_all_lines + static_cast<size_t>(line_idx) * dstWidth, fpath_line, static_cast<size_t>(dstWidth) * sizeof(int));
                }

                queue.enqueue_unmap_buffer(_ccosts, ccosts - d->mdisVector);
            }

            queue.enqueue_write_buffer(d->fpath_gpu, 0, static_cast<size_t>(dstWidth) * (dstHeight / 2) * sizeof(int), fpath_all_lines);

            const size_t apply_gws[2] = { static_cast<size_t>((dstWidth + 15) & -16), static_cast<size_t>((dstHeight / 2 + 15) & -16) };
            const size_t apply_lws[2] = { 16, 16 };

            d->applyInterpolationKernel.set_args(srcImage, dstImage, d->fpath_gpu, d->dmap_gpu, dstWidth, dstHeight, field_n, static_cast<cl_int>(d->ucubic), static_cast<cl_int>(use_dh), peak);
            queue.enqueue_nd_range_kernel(d->applyInterpolationKernel, 2, nullptr, apply_gws, apply_lws);

            boost::compute::image2d* final_image{ &dstImage };

            if (d->vcheck)
            {
                boost::compute::buffer sclip_buffer;
                int sclip_stride_arg{};
                cl_int use_sclip{};

                if (d->sclip)
                {
                    use_sclip = 1;
                    const void* scpp{ avs_get_read_ptr_p(scp, current_plane) };
                    sclip_stride_arg = avs_get_pitch_p(scp, current_plane);
                    const size_t sclip_size{ static_cast<size_t>(sclip_stride_arg) * avs_get_height_p(scp, current_plane) };
                    sclip_buffer = boost::compute::buffer(queue.get_context(), sclip_size, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, const_cast<void*>(scpp));
                }
                else
                    sclip_buffer = boost::compute::buffer(queue.get_context(), 1, CL_MEM_READ_ONLY);

                d->vCheckKernel.set_args(srcImage, dstImage, vcheckTmpImage, d->dmap_gpu, sclip_buffer, use_sclip, dstWidth, dstHeight, sclip_stride_arg, field_n, d->vcheck, d->vthresh2, d->rcpVthresh0, d->rcpVthresh1, d->rcpVthresh2, peak);
                queue.enqueue_nd_range_kernel(d->vCheckKernel, 2, nullptr, apply_gws, apply_lws);
                final_image = &vcheckTmpImage;
            }

            queue.enqueue_read_image(*final_image, boost::compute::dim(0, 0), boost::compute::dim(dstWidth, dstHeight), _dstp, avs_get_pitch_p(dst, current_plane));

            avs_release_video_frame(pad);
        }
    }

    d->queue.finish();
}

template void filterCL_avx512<uint8_t>(const AVS_VideoFrame* __restrict src, const AVS_VideoFrame* __restrict scp, AVS_VideoFrame* __restrict dst, const int field_n, bool use_dh, EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
template void filterCL_avx512<uint16_t>(const AVS_VideoFrame* __restrict src, const AVS_VideoFrame* __restrict scp, AVS_VideoFrame* __restrict dst, const int field_n, bool use_dh, EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
template void filterCL_avx512<float>(const AVS_VideoFrame* __restrict src, const AVS_VideoFrame* __restrict scp, AVS_VideoFrame* __restrict dst, const int field_n, bool use_dh, EEDI3CLData* __restrict d, const AVS_FilterInfo* __restrict fi);
