
To build on Mac, you should have brew installed. Then, by running  
`brew install gcc boost`  
the dependencies will be satisfied. Build with `make brew`

To build on Linux, install `boost`, an MPI library (and a C++ compiler if missing) and then adapt the makefile. A plain and simple `make` could also work, if you're lucky.

You can find `make` solutions already existing for many HPC systems available in Italy and Europe.

Debug and profiling `make` solutions, ready or useful to adapt to your system, are also available.

To build on Windows (not Cygwin which relies on the makefile and the dependencies satisfied through its installer), you should have Visual Studio 14 (2015) installed and MS-MPI.
You can download MS-MPI v7 from [here](https://www.microsoft.com/en-us/download/details.aspx?id=49926)

