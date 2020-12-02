# C++ script to get event yield and make pileup histograms

Requires the following header file <https://raw.githubusercontent.com/nlohmann/json/master/single_include/nlohmann/json.hpp>


Compilation flags:

```
g++ -std=c++11 getEventYields.cc `root-config --cflags --glibs`  -o getEventYields
```

Run script:

```
mkdir test
./getEventYields /vols/cms/LLP/files_201117/analysis/2016  test
```

