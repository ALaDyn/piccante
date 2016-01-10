
To build on Mac, you should have brew installed. Then, by running  
`brew install gcc gsl boost`  
the dependencies will be satisfied. Build with `make brew`

To build on Linux, install `boost`, `gsl`, an MPI library (and a C++ compiler if missing) and then adapt the makefile. A plain and simple `make` could also work, if you're lucky.

You can find `make` solutions already existing for many HPC systems available in Italy and Europe.

Debug and profiling `make` solutions, ready or useful to adapt to your system, are also available.

To build on Windows (not Cygwin which relies on the makefile and the dependencies satisfied through its installer), you should have Visual Studio 14 (2015) installed, MS-MPI and GnuWin32.
You can download MS-MPI v7 from [here](https://www.microsoft.com/en-us/download/details.aspx?id=49926), while GnuWin can be installed through [Chocolatey](https://chocolatey.org/)

GnuWin doesn't provide required .lib files, so from an elevated PowerShell you should do

```
PS C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC> .\vcvarsall.bat

PS C:\Program Files (x86)\GnuWin32\lib> &("C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\lib.exe") /def:libgsl.def /machine:x64
   Microsoft (R) Library Manager Version 14.00.23506.0
   Copyright (C) Microsoft Corporation.  All rights reserved.
   Creating library libgsl.lib and object libgsl.exp

PS C:\Program Files (x86)\GnuWin32\lib> &("C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\lib.exe") /def:libgslcblas.def /machine:x64
   Microsoft (R) Library Manager Version 14.00.23506.0
   Copyright (C) Microsoft Corporation.  All rights reserved.
   Creating library libgslcblas.lib and object libgslcblas.exp
 ```
