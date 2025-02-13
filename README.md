# SZp (Also known as fZ-light)

* Major author and developer: Jiajun Huang, Sheng Di
* Supervisor: Sheng Di

This is the official repository of SZp, an ultra-fast error-bounded lossy compressor: CPU version (supporting openMP). 
The design and optimizations of SZp are published under the name fZ-light in SC '24.

## Installation
Configure and build the SZp:
```bash
# Decompress the downloaded files
# Change to the SZp directory and set up the installation directory:
cd SZp
mkdir install

# Run the configuration script:
./configure --prefix=$(pwd)/install --enable-openmp

# Compile the SZp using multiple threads:
make -j

# Install the compiled SZp:
make install

```
