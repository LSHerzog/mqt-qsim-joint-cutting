# Base OS
FROM debian

# Set non-interactive mode to prevent prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required dependencies
RUN apt-get update && apt-get install -y \
    g++ make wget cmake libblas-dev liblapack-dev git && \
    rm -rf /var/lib/apt/lists/*

# Install xtl (dependency for xtensor)
RUN git clone https://github.com/xtensor-stack/xtl.git && \
    cd xtl && \
    cmake . -DCMAKE_INSTALL_PREFIX=/usr/local && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf xtl

# Install xtensor
RUN git clone https://github.com/xtensor-stack/xtensor.git && \
    cd xtensor && \
    cmake . -DCMAKE_INSTALL_PREFIX=/usr/local && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf xtensor

# Install xtensor-blas
RUN git clone https://github.com/xtensor-stack/xtensor-blas.git && \
    cd xtensor-blas && \
    cmake . -DCMAKE_INSTALL_PREFIX=/usr/local && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf xtensor-blas

# Copy relevant files for simulation
COPY ./Makefile /qsim/Makefile
COPY ./apps/ /qsim/apps/
COPY ./circuits/ /qsim/circuits/
COPY ./lib/ /qsim/lib/

# Set working directory
WORKDIR /qsim/

# Compile qsim (should work now with xtensor installed)
RUN make qsim

# Define entry point
ENTRYPOINT ["/qsim/apps/qsim_base.x"]
