
# BrainSuite Diffusion Pipeline (BDP)

BDP enables fusion of diffusion and structural MRI information for advanced
image and connectivity analysis. It provides various methods for distortion 
correction, co-registration, diffusion modeling (DTI and ODF) and basic ROI-wise
statistic. BDP is flexible and diverse tool which supports wide variety of 
diffusion datasets. For more information and detailed documentation, please 
see BDP website: http://brainsuite.org/processing/diffusion/

This is release 25a (06-Mar-2025) of BDP.

Created by Chitresh Bhushan, Divya Varadarajan, Justin P. Haldar,
           Anand A. Joshi,  David W. Shattuck, and Richard M. Leahy



# Downloading BDP

Starting with version 25a, the latest pre-compiled binaries and source code of BDP can be obtained from http://github.com/BrainSuite/BDP/releases/. Compiled binaries are available for Linux (Intel x64), Mac (Intel x64 and Apple ARM64), and Windows (Intel x64).

Running the compiled binaries requires either the 2024b MATLAB Runtime (free) or an installation of MATLAB2024b (paid).

Source code packages for all released versions of BDP are available from http://github.com/BrainSuite/BDP/.

Older archived pre-compiled binaries and source code of BDP can be obtained from http://brainsuite.org/. 



# Compiling BDP

BDP is written in a combination of MATLAB, C, and C++ code. So, there are two steps
in the compilation process:(1) compiling c/mex files, and (2) compiling MATLAB code
to generate executables and packages. Both steps are performed from the MATLAB
command prompt.

BDP requires working installations of the following for building the binaries:
   - MATLAB 2024b: http://www.mathworks.com/products/matlab/
   - Image Processing Toolbox: http://www.mathworks.com/products/image/
   - MATLAB Compiler: http://www.mathworks.com/products/compiler/
   - MATLAB Parallel Computing Toolbox: https://www.mathworks.com/products/parallel-computing/

   - A supported C/C++ Compiler: https://www.mathworks.com/support/requirements/supported-compilers.html

First open MATLAB and navigate to the `<BDP>/packaging_tools` folder of the 
source code. Then, Compile the c/mex files by running following on the MATLAB
command prompt:

```
>> CompileBDPMexFiles
```

It will compile all required mex functions. Note that you should change the 
generated Makefile in Linux and Mac to use clang instead of gcc. 

Finally, compile the MATLAB code and generate packages by running following:

```
>> CompileBrainSuiteDiffusionPipeline --package 99aRC9 --build 9999
```
 
where, `99aRC9` is the name of release and `9999` is four digit build number.
This will generate relevant tarball/zip for distribution. Note that above step
requires `mcc` command from MATLAB Compiler toolbox.



# BDP License 
The main BDP files are licensed under [GPL-V2-only](https://spdx.org/licenses/GPL-2.0-only.html).
Please see the included license file(s) for more details.


