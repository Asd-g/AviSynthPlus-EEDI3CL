## Description

EEDI3 works by finding the best non-decreasing (non-crossing) warping between two lines by minimizing a cost functional. The cost is based on neighborhood similarity (favor connecting regions that look similar), the vertical difference created by the interpolated values (favor small differences), the interpolation directions (favor short connections vs long), and the change in interpolation direction from pixel to pixel (favor small changes).

This is [a port of the VapourSynth plugin EEDI3CL](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-EEDI3).

### Requirements:

- AviSynth+ r3688 or later

- Microsoft VisualC++ Redistributable Package 2022 (can be downloaded from [here](https://github.com/abbodi1406/vcredist/releases))

### Usage:

```
EEDI3CL(clip input, int "field", bool "dh", int[] "planes", float "alpha", float "beta", float "gamma", int "nrad", int "mdis", bool "hp", bool "ucubic", bool "cost3", int "vcheck", float "vthresh0", float "vthresh1", float "vthresh2", clip "sclip", int "opt", int "device", bool "list_device", bool "info", bool "luma", bool "dw")
```

### Parameters:

- input<br>
    A clip to process.<br>
    It must be in 8..32-bit planar format.

- field<br>
    Controls the mode of operation (double vs same rate) and which field is kept.<br>
    -2: Double rate (alternates each frame), `_FieldBased` frame property order or if `_FieldBased` is `0`, or missing - AviSynth internal order.<br>
    -1: Same rate, `_FieldBased` frame property order or if `_FieldBased` is `0`, or missing - AviSynth internal order.<br>
    0: Same rate, keep bottom field.<br>
    1: Same rate, keep top field.<br>
    2: Double rate (alternates each frame), starts with bottom.<br>
    3: Double rate (alternates each frame), starts with top.<br>
    Default: -1.

- dh<br>
    Doubles the height of the input.<br>
    Each line of the input is copied to every other line of the output and the missing lines are interpolated.<br>
    If field=0, the input is copied to the odd lines of the output.<br>
    If field=1, the input is copied to the even lines of the output.<br>
    Field must be set to either 0 or 1 when using dh=true.<br>
    Default: False.

- planes<br>
    Sets which planes will be processed.<br>
    Planes that are not processed will contain uninitialized memory.<br>
    Default: [0, 1, 2, 3].

- alpha/beta/gamma<br>
    These trade off line/edge connection vs artifacts created.<br>
    `alpha` and `beta` must be in the range [0,1], and the sum `alpha`+`beta` must be in the range [0,1].<br>
    `alpha` is the weight given to connecting similar neighborhoods. The larger `alpha` is the more lines/edges should be connected.<br>
    `beta` is the weight given to vertical difference created by the interpolation. The larger `beta` is the less edges/lines will be connected (at 1.0 you get no edge directedness at all).<br>
    The remaining weight (1.0-`alpha`-`beta`) is given to interpolation direction (large directions (away from vertical) cost more). So the more weight you have here the more shorter connections will be favored.<br>
    Finally, `gamma` penalizes changes in interpolation direction, the larger `gamma` is the smoother the interpolation field between two lines (range is [0,inf].<br>
    If lines aren't getting connected then increase `alpha` and maybe decrease `beta`/`gamma`. Go the other way if you are getting unwanted artifacts.<br>
    Default: alpha = 0.2, beta = 0.25, gamma = 20.0.

- nrad/mdis<br>
    `nrad` sets the radius used for computing neighborhood similarity. Valid range is [0,3].<br>
    `mdis` sets the maximum connection radius. Valid range is [1,40].<br>
    If `mdis`=20, then when interpolating pixel (50,10) (x,y), the farthest connections allowed would be between (30,9)/(70,11) and (70,9)/(30,11).<br>
    Larger `mdis` will allow connecting lines of smaller slope, but also increases the chance of artifacts.<br>
    Larger `mdis` will be slower.<br>
    Larger `nrad` will be slower.<br>
    Default: nrad = 2, mdis = 20.

- hp/ucubic/cost3<br>
    These are speed vs quality options.<br>
    `hp=True` - use half pel steps, hp=False - use full pel steps. Currently only full pel is implemented and this parameter has no effect.<br>
    `ucubic=True` - use cubic 4 point interpolation, `ucubic=False` - use 2 point linear interpolation.<br>
    `cost3=True` - use 3 neighborhood cost function to define similarity, `cost3=False` - use 1 neighborhood cost function.<br>
    Default: hp = False, ucubic = True, cost3 = True.

- vcheck/vthresh0/vthresh1/vthresh2/sclip<br>
    ```
      vcheck settings:

          0 - no reliability check
          1 - weak reliability check
          2 - med reliability check
          3 - strong reliability check

      If vcheck is greater than 0, then the resulting interpolation is checked for reliability/consistency. Assume
      we interpolated pixel 'fh' below using dir=4 (i.e. averaging pixels bl and cd).

           aa ab ac ad ae af ag ah ai aj ak al am an ao ap
                                eh          el
           ba bb bc bd be bf bg bh bi bj bk bl bm bn bo bp
                    fd          fh          fl
           ca cb cc cd ce cf cg ch ci cj ck cl cm cn co cp
                    gd          gh
           da db dc dd de df dg dh di dj dk dl dm dn do dp

      When checking pixel 'fh' the following is computed:

            d0 = abs((el+fd)/2 - bh)
            d1 = abs((fl+gd)/2 - ch)

            q2 = abs(bh-fh)+abs(ch-fh)
            q3 = abs(el-bl)+abs(fl-bl)
            q4 = abs(fd-cd)+abs(gd-cd)

            d2 = abs(q2-q3)
            d3 = abs(q2-q4)

            mdiff0 = vcheck == 1 ? min(d0,d1) : vcheck == 2 ? ((d0+d1+1)>>1) : max(d0,d1)
            mdiff1 = vcheck == 1 ? min(d2,d3) : vcheck == 2 ? ((d2+d3+1)>>1) : max(d2,d3)

            a0 = mdiff0/vthresh0;
            a1 = mdiff1/vthresh1;
            a2 = max((vthresh2-abs(dir))/vthresh2,0.0f)

            a = min(max(max(a0,a1),a2),1.0f)

            final_value = (1.0-a)*fh + a*cint


        ** If sclip is supplied, cint is the corresponding value from sclip. If sclip isn't supplied,
           then vertical cubic interpolation is used to create it.
    ```

    Default: vcheck = 2, vthresh0 = 32.0, vthresh1 = 64.0, vthresh2 = 4.0, sclip = not specified.

- opt<br>
    Sets which CPU optimizations to use for the pathfinding algorithm.<br>
    -1: Auto-detect.<br>
    0: Use C++ code.<br>
    1: Use SSE2 code.<br>
    2: Use AVX2 code.<br>
    3: Use AVX512 code.<br>
    Default: -1.

- device<br>
    Sets target OpenCL device.<br>
    Use list_device to get the index of the available devices.<br>
    By default the default device is selected.

- list_device<br>
    Whether to draw the devices list on the frame.

- info<br>
    Whether to draw the OpenCL-related info on the frame.

- luma<br>
    Whether the format of the output video is Y when only luma plane is processed.<br>
    It has effect only for YUV clips.<br>
    Default: False.

- dw<br>
    Doubles the width of the input.<br>
    It performs a horizontal interpolation, analogous to how `dh` performs vertical interpolation. This mode is significantly slower as it involves transposing the frame and running the core algorithm twice.<br>
    `field` must not be 2 or 3 when `dw=true`.<br>
    `sclip` is not supported when `dw=true`.<br>
    Default: False.

### Building:

#### Prerequisites
- **Git**
- **CMake** >= 3.25
- A **C++20 capable compiler** (e.g., Visual Studio 2022, GCC 11+, Clang 12+)
- **Boost** libraries (chrono, filesystem, system). Must be findable by CMake.
- An **OpenCL SDK**. Must be findable by CMake.

1.  (Linux) Install prerequisites (example for Debian/Ubuntu):

    ```
    sudo apt-get install build-essential git cmake libboost-chrono-dev libboost-filesystem-dev libboost-system-dev ocl-icd-opencl-dev
    ```

2.  Clone the repository:

    ```
    git clone --depth 1 --shallow-submodules --recursive https://github.com/Asd-g/AviSynthPlus-EEDI3CL
    cd AviSynthPlus-EEDI3CL
    ```

3.  Configure and build the project:

    ```
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j$(nproc)
    ```

4.  (Linux) Install the plugin (optional):

    ```
    sudo cmake --install build
    ```
