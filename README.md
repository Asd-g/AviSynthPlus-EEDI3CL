## Description

EEDI3 works by finding the best non-decreasing (non-crossing) warping between two lines by minimizing a cost functional. The cost is based on neighborhood similarity (favor connecting regions that look similar), the vertical difference created by the interpolated values (favor small differences), the interpolation directions (favor short connections vs long), and the change in interpolation direction from pixel to pixel (favor small changes).

This is [a port of the VapourSynth plugin EEDI3CL](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-EEDI3).

### Requirements:

- AviSynth+ r3682 (can be downloaded from [here](https://gitlab.com/uvz/AviSynthPlus-Builds) until official release is uploaded) (r3689 recommended) or later

- Microsoft VisualC++ Redistributable Package 2022 (can be downloaded from [here](https://github.com/abbodi1406/vcredist/releases))

### Usage:

```
EEDI3CL(clip input, int "field", bool "dh", int[] "planes", float "alpha", float "beta", float "gamma", int "nrad", int "mdis", bool "hp", bool "ucubic", bool "cost3", int "vcheck", float "vthresh0", float "vthresh1", float "vthresh2", clip "sclip", int "opt", int "device", bool "list_device", bool "info", bool "luma")
```

### Parameters:

- input\
    A clip to process.\
    It must be in 8..32-bit planar format.

- field\
    Controls the mode of operation (double vs same rate) and which field is kept.\
    -2: Double rate (alternates each frame), `_FieldBased` frame property order or if `_FieldBased` is `0`, or missing - AviSynth internal order.\
    -1: Same rate, `_FieldBased` frame property order or if `_FieldBased` is `0`, or missing - AviSynth internal order.\
    0: Same rate, keep bottom field.\
    1: Same rate, keep top field.\
    2: Double rate (alternates each frame), starts with bottom.\
    3: Double rate (alternates each frame), starts with top.\
    Default: -1.

- dh\
    Doubles the height of the input.\
    Each line of the input is copied to every other line of the output and the missing lines are interpolated.\
    If field=0, the input is copied to the odd lines of the output.\
    If field=1, the input is copied to the even lines of the output.\
    Field must be set to either 0 or 1 when using dh=true.\
    Default: False.

- planes\
    Sets which planes will be processed.\
    Planes that are not processed will contain uninitialized memory.\
    Default: [0, 1, 2, 3].

- alpha/beta/gamma\
    These trade off line/edge connection vs artifacts created.\
    `alpha` and `beta` must be in the range [0,1], and the sum `alpha`+`beta` must be in the range [0,1].\
    `alpha` is the weight given to connecting similar neighborhoods. The larger `alpha` is the more lines/edges should be connected.\
    `beta` is the weight given to vertical difference created by the interpolation. The larger `beta` is the less edges/lines will be connected (at 1.0 you get no edge directedness at all).\
    The remaining weight (1.0-`alpha`-`beta`) is given to interpolation direction (large directions (away from vertical) cost more). So the more weight you have here the more shorter connections will be favored.\
    Finally, `gamma` penalizes changes in interpolation direction, the larger `gamma` is the smoother the interpolation field between two lines (range is [0,inf].\
    If lines aren't getting connected then increase `alpha` and maybe decrease `beta`/`gamma`. Go the other way if you are getting unwanted artifacts.\
    Default: alpha = 0.2, beta = 0.25, gamma = 20.0.

- nrad/mdis\
    `nrad` sets the radius used for computing neighborhood similarity. Valid range is [0,3].\
    `mdis` sets the maximum connection radius. Valid range is [1,40].\
    If `mdis`=20, then when interpolating pixel (50,10) (x,y), the farthest connections allowed would be between (30,9)/(70,11) and (70,9)/(30,11).\
    Larger `mdis` will allow connecting lines of smaller slope, but also increases the chance of artifacts.\
    Larger `mdis` will be slower.\
    Larger `nrad` will be slower.\
    Default: nrad = 2, mdis = 20.

- hp/ucubic/cost3\
    These are speed vs quality options.\
    `hp=True` - use half pel steps, hp=False - use full pel steps. Currently only full pel is implemented and this parameter has no effect.\
    `ucubic=True` - use cubic 4 point interpolation, `ucubic=False` - use 2 point linear interpolation.\
    `cost3=True` - use 3 neighborhood cost function to define similarity, `cost3=False` - use 1 neighborhood cost function.\
    Default: hp = False, ucubic = True, cost3 = True.

- vcheck/vthresh0/vthresh1/vthresh2/sclip\
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

- opt\
    Sets which cpu optimizations to use.\
    -1: Auto-detect.\
    0: Use C++ code.\
    1: Use SSE2 code.\
    Default: -1.

- device\
    Sets target OpenCL device.\
    Use list_device to get the index of the available devices.\
    By default the default device is selected.

- list_device\
    Whether to draw the devices list on the frame.

- info\
    Whether to draw the OpenCL-related info on the frame.

- luma\
    Whether the format of the output video is Y when only luma plane is processed.\
    It has effect only for YUV clips.\
    Default: False.

### Building:

- Requires `Boost` and `OpenCL`.

- Windows\
    Use solution files.

- Linux
    ```
    Requirements:
        - Git
        - C++17 compiler
        - CMake >= 3.16
    ```
    ```
    git clone https://github.com/Asd-g/AviSynthPlus-EEDI3CL && \
    cd AviSynthPlus-EEDI3CL && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    sudo make install
    ```
