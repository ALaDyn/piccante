A useful practical guide to building `piccante` is also available looking at the CI recipes (both travis and appveyor).


## Ubuntu
Open a bash terminal and write  
```
$  sudo apt-get update
$  sudo apt-get install g++ cmake make libboost-all-dev git openmpi-bin openmpi-doc libopenmpi-dev
$  git clone https://github.com/ALaDyn/piccante.git
$  cd piccante
$  mkdir build && cd build
$  cmake .. ; cmake --build . --target install
```


## macOS
1) Install XCode Command Line Tools, if not already installed, with this command in Terminal:
```
$  xcode-select --install
```
2) If not already done, install Homebrew following the [guide](https://brew.sh/index_it.html).  
3) Open a terminal and write
```
$  brew update
$  brew install cmake make boost git open-mpi
$  git clone https://github.com/ALaDyn/piccante.git
$  cd piccante
$  mkdir build && cd build
$  cmake .. ; cmake --build . --target install
```


## Windows (7+) [native Win32]
1) Install or Update Visual Studio. We need at least VS 2015.3 or newer. If not installed, we suggest [Visual Studio 2017 Community](http://visualstudio.com)   
2) If not already done, install [chocolatey](http://chocolatey.org)   
3) Open a Powershell with admin privileges and write
```
PS \>              cinst -y git cmake pscx powershell
```
4) Reboot the PC if requested
5) Let's define a work folder, which will contain also the piccante source tree. We will call it WORKSPACE going on: it could be a "Code" folder in our home, a "codes" folder on our desktop. This path will become a reference for us thanks to an environment variable definition. Open a Powershell without admin privileges and write
```
PS \>              rundll32 sysdm.cpl,EditEnvironmentVariables
```
6) In the window, in the upper part, create a new variable with name WORKSPACE and with value the full path of our work folder defined in the previous step. Add also to the "Path" variable this path (be sure to have a `;` delimiter from other records)
```
%PROGRAMFILES%/CMake/bin;
```
7) If vcpkg is not already installed, please follow the next procedure, otherwise jump to step 9. Close the Powershell and re-open it, again without admin privileges
```
PS \>              cd $env:WORKSPACE
PS Code>           git clone https://github.com/Microsoft/vcpkg.git
PS Code>           cd vcpkg
PS Code\vcpkg>     .\bootstrap-vcpkg.bat 
```
8) Close the Powershell and re-open it with admin privileges
```
PS \>              cd $env:WORKSPACE
PS Code\>          Invoke-WebRequest https://download.microsoft.com/download/B/2/E/B2EB83FE-98C2-4156-834A-E1711E6884FB/MSMpiSetup.exe -OutFile $env:WORKSPACE\msmpi.exe
PS Code>           .\msmpi -unattend
PS Code>           cd vcpkg
PS Code\vcpkg>     .\vcpkg integrate install

```
9) Close the Powershell and re-open it without admin privileges
```
PS \>              cd $env:WORKSPACE
PS Code>           cd vcpkg
PS Code\vcpkg>     .\vcpkg install boost msmpi
```
10) Open a text editor (`notepad.exe`is ok!) and paste the following text, depending on your configuration:

**Visual Studio 2015, Windows 32 bit**
```
Import-Module Pscx
Invoke-BatchFile "${env:PROGRAMFILES}\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86
```
**Visual Studio 2015, Windows 64 bit**
```
Import-Module Pscx
Invoke-BatchFile "${env:PROGRAMFILES(x86)}\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x64
```
**Visual Studio 2017, Windows 32 bit**
```
Import-Module Pscx
Invoke-BatchFile "${env:PROGRAMFILES}\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvarsall.bat" x86
```
**Visual Studio 2017, Windows 64 bit**
```
Import-Module Pscx
Invoke-BatchFile "${env:PROGRAMFILES(x86)}\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
```
Save it in the folder Documents\WindowsPowerShell with the filename `Microsoft.PowerShell_profile.ps1`

11) Close Powershell, reopen it without admin privileges
```
PS \>              cd $env:WORKSPACE
PS Code>           git clone https://github.com/ALaDyn/piccante.git
PS Code>           cd piccante
PS Code\piccante>  mkdir build ; cd build
PS piccante\build> cmake .. "-DCMAKE_TOOLCHAIN_FILE=$env:WORKSPACE\vcpkg\scripts\buildsystems\vcpkg.cmake"
PS piccante\build> cmake --build . --target install --config Release
```


## Windows Subsystem for Linux (WSL)
It's better to use Bash on Ubuntu on Windows only on builds 15063+ (Creators Update and beyond).  
1) If not enabled, activate UoW following the [official guide](https://msdn.microsoft.com/it-it/commandline/wsl/install_guide)   
2) Follow the Ubuntu guide


## Cygwin
1) If not already done, install [chocolatey](http://chocolatey.org)   
2) Open a Powershell with admin privileges and write
```
PS \>              cinst -y cygwin
```
3) Close the Powershell, reopen it without admin privileges
```
PS \>              cd $env:WORKSPACE
PS Code>           Invoke-WebRequest https://cygwin.com/setup-x86_64.exe -OutFile $env:WORKSPACE\cygwin-setup.exe
PS Code>           .\cygwin-setup --quiet-mode --no-shortcuts --no-startmenu --no-desktop --upgrade-also --packages gcc-g++,libopenmpi-devel,cmake,libboost-devel
```
4) Open a Cygwin shell
```
$  cd $WORKSPACE
$  git clone https://github.com/ALaDyn/piccante.git
$  cd piccante
$  mkdir build && cd build
$  cmake .. ; cmake --build . --target install
```



