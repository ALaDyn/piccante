# Replace system-provided `jsoncpp`

In case you have a faulty `jsoncpp` system-wide installation, or in case you have an offline cluster, you need to manually supply another `jsoncpp` installation, otherwise CMake could fail.

If you are offline, please adapt these instructions to your situation. You just have to replace `wget` call with your available way to obtain that file in that folder (`scp`, ...)

```bash
#!/bin/bash

VERSION="1.8.1"

cd ~
mkdir jsoncpp
cd jsoncpp
export JSONCPP_INSTALL_FOLDER=$(pwd)
wget https://github.com/open-source-parsers/jsoncpp/archive/${VERSION}.zip
unzip ${VERSION}.zip
cd jsoncpp-${VERSION}/
mkdir build ; cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=../.. ..
cmake --build . --target install


echo $JSONCPP_INSTALL_FOLDER
# build now piccante passing the JSONCPP_INSTALL_FOLDER
cd ~
cd piccante
mkdir build ; cd build
cmake .. -DCMAKE_PREFIX_PATH=${JSONCPP_INSTALL_FOLDER}
cmake --build . --target install
```
