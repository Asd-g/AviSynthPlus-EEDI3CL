#include "EEDI3CL.h"

template<typename T>
void filterCL_sse2(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, AVS_VideoFrame** pad, const int field_n, const EEDI3CLData* const __restrict d)
{
    constexpr int planes_y[4]{ AVS_PLANAR_Y, AVS_PLANAR_U, AVS_PLANAR_V, AVS_PLANAR_A };
    constexpr int planes_r[4]{ AVS_PLANAR_R, AVS_PLANAR_G, AVS_PLANAR_B, AVS_PLANAR_A };
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
            float* pcosts{ d->pcosts + d->mdisVector };
            int* _pbackt{ d->pbackt + d->mdisVector };
            int* fpath{ d->fpath };
            int* _dmap{ d->dmap };
            int* tline{ d->tline };

            const size_t globalWorkSize[] = { static_cast<size_t>((dstWidth + 15) & -16), static_cast<size_t>(d->vectorSize) };
            constexpr size_t localWorkSize[] = { 16, 4 };
            const int bufferSize{ static_cast<int>(dstWidth * d->tpitchVector * sizeof(cl_float)) };

            copyPlane<T>(_dstp + dstStride * (1 - field_n), dstStride * 2, _srcp + srcStride * (4 + 1 - field_n) + 12, srcStride * 2, dstWidth, dstHeight / 2);

            queue.enqueue_write_image(srcImage, boost::compute::dim(0, 0), boost::compute::dim(srcWidth, srcHeight), _srcp, avs_get_pitch_p(pad[plane], planes[0]));

            for (int y{ 4 + field_n }; y < srcHeight - 4; y += 2 * d->vectorSize)
            {
                const int off{ (y - 4 - field_n) >> 1 };

                calculateConnectionCosts.set_args(srcImage, _ccosts, dstWidth, srcHeight - 4, y);
                queue.enqueue_nd_range_kernel(calculateConnectionCosts, 2, nullptr, globalWorkSize, localWorkSize);

                float* ccosts{ reinterpret_cast<float*>(queue.enqueue_map_buffer(_ccosts, CL_MAP_READ, 0, bufferSize)) + d->mdisVector };

                // calculate path costs
                Vec4f().load(ccosts).store_a(pcosts);
                for (int x{ 1 }; x < dstWidth; ++x)
                {
                    const float* tT{ ccosts + d->tpitchVector * x };
                    const float* ppT{ pcosts + d->tpitchVector * (x - 1) };
                    float* pT{ pcosts + d->tpitchVector * x };
                    int* piT{ _pbackt + d->tpitchVector * (x - 1) };

                    const int umax{ std::min({ x, dstWidth - 1 - x, d->mdis }) };
                    const int umax2{ std::min({ x - 1, dstWidth - x, d->mdis }) };
                    for (int u{ -umax }; u <= umax; ++u)
                    {
                        Vec4i idx{ zero_si128() };
                        Vec4f bval{ FLT_MAX };

                        for (int v{ std::max(-umax2, u - 1) }; v <= std::min(umax2, u + 1); ++v)
                        {
                            const Vec4f z{ Vec4f().load_a(ppT + v * d->vectorSize) + d->gamma * std::abs(u - v) };
                            const Vec4f ccost{ min(z, FLT_MAX * 0.9f) };
                            idx = select(Vec4ib(ccost < bval), v, idx);
                            bval = min(ccost, bval);
                        }

                        const Vec4f z{ bval + Vec4f().load(tT + u * d->vectorSize) };
                        min(z, FLT_MAX * 0.9f).store_a(pT + u * d->vectorSize);
                        idx.store_nt(piT + u * d->vectorSize);
                    }
                }

                for (int vs{ 0 }; vs < d->vectorSize; ++vs)
                {
                    const int realY{ 4 + field_n + 2 * (off + vs) };
                    if (realY >= srcHeight - 4)
                        break;

                    const T* srcp{ _srcp + srcStride * realY + 12 };
                    T* dstp{ _dstp + dstStride * (field_n + 2 * (off + vs)) };
                    int* dmap{ _dmap + dstWidth * (off + vs) };

                    const T* src3p{ srcp - srcStride * 3 };
                    const T* src1p{ srcp - srcStride };
                    const T* src1n{ srcp + srcStride };
                    const T* src3n{ srcp + srcStride * 3 };

                    const int* pbackt{ _pbackt + vs };

                    // backtrack
                    fpath[dstWidth - 1] = 0;
                    for (int x{ dstWidth - 2 }; x >= 0; --x)
                        fpath[x] = pbackt[(d->tpitch * x + fpath[x + 1]) * d->vectorSize];

                    interpolate<T>(src3p, src1p, src1n, src3n, nullptr, fpath, dmap, dstp, dstWidth, d->ucubic, d->peak);
                }

                queue.enqueue_unmap_buffer(_ccosts, ccosts - d->mdisVector);
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

template void filterCL_sse2<uint8_t>(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, AVS_VideoFrame** pad, const int field_n, const EEDI3CLData* const __restrict d);
template void filterCL_sse2<uint16_t>(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, AVS_VideoFrame** pad, const int field_n, const EEDI3CLData* const __restrict d);
template void filterCL_sse2<float>(const AVS_VideoFrame* src, const AVS_VideoFrame* scp, AVS_VideoFrame* dst, AVS_VideoFrame** pad, const int field_n, const EEDI3CLData* const __restrict d);
