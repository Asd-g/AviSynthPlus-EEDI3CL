##### 1.1.0:
    Added `dw` parameter (#2).
    Fixed `luma=true` for bit depth >8.
    Changed the cost calculation when `ucubic=true` - now cubic interpolation is used for every stage.
    Fixed boundary handling.
    Reduced CPU-GPU transfers - moved the final interpolation and vCheck steps to GPU.
    Added AVX2 and AVX-512 code.
    Changed to dynamic linking against avs.
    Updated CMake building.
    Fixed memory leak when `info=true` or `list_device=true`.

##### 1.0.2:
    Fixed crashing when unsupported Avs+ used by explicitly throwing error.
    Changed the required Avs+ version.

##### 1.0.1:
    Added parameter `luma`.

##### 1.0.0:
    Initial release. (port of the vs plugin r4)
