# OmpSs-2 environment: ovni + Nanos6 4.3.1 + LLVM/Clang (bsc-pm, github-release-2026.06)
# x86_64 only.
#
# Build:   docker build -t ompss-env:v1 .
#          docker build --build-arg JOBS=8 -t ompss-env:v1 .   # more parallelism if you have the RAM
# Run:     docker run -it --rm ompss-env:v1

FROM --platform=linux/amd64 debian:stable-slim

# Parallel compile jobs. LLVM needs ~2 GB RAM per job; keep this low on small machines.
ARG JOBS=4
ARG NANOS6_VERSION=4.3.1
ARG LLVM_RELEASE=github-release-2026.06
# Pin ovni to a tag or commit for reproducible builds (e.g. --build-arg OVNI_REF=<tag>).
ARG OVNI_REF=master

ENV OVNI_HOME=/root/installation/ovni \
    INSTALLATION_PREFIX=/root/installation/ompss-2
ENV PATH=$INSTALLATION_PREFIX/bin:$OVNI_HOME/bin:$PATH \
    LD_LIBRARY_PATH=$INSTALLATION_PREFIX/lib:$OVNI_HOME/lib

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake ninja-build git wget ca-certificates \
      autoconf automake libtool pkg-config python3 \
      libopenmpi-dev openmpi-bin libnuma-dev libboost-dev \
      gdb vim \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /home

# --- ovni ---------------------------------------------------------------
RUN git clone https://github.com/bsc-pm/ovni.git \
    && cd ovni && git checkout ${OVNI_REF} \
    && cmake -S . -B build -G Ninja -DCMAKE_INSTALL_PREFIX=$OVNI_HOME \
    && ninja -C build -j${JOBS} install \
    && cd /home && rm -rf ovni

# --- Nanos6 (first pass, before the compiler exists) ---------------------
RUN wget -q https://github.com/bsc-pm/nanos6/archive/refs/tags/version-${NANOS6_VERSION}.tar.gz \
    && tar -xzf version-${NANOS6_VERSION}.tar.gz && rm version-${NANOS6_VERSION}.tar.gz \
    && cd nanos6-version-${NANOS6_VERSION} \
    && ./autogen.sh \
    && ./configure --prefix=$INSTALLATION_PREFIX --with-boost \
         --with-ovni=$OVNI_HOME --with-hwloc=embedded --with-jemalloc=embedded \
    && make -j${JOBS} && make install

# --- LLVM / Clang with OmpSs-2 support ----------------------------------
# X86 target only (much less RAM/time); LLVM_PARALLEL_LINK_JOBS needs Ninja.
RUN wget -q https://github.com/bsc-pm/llvm/archive/refs/tags/${LLVM_RELEASE}.tar.gz \
    && tar -xzf ${LLVM_RELEASE}.tar.gz && rm ${LLVM_RELEASE}.tar.gz \
    && cd llvm-${LLVM_RELEASE} \
    && cmake -S llvm -B build -G Ninja \
         -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_INSTALL_PREFIX=$INSTALLATION_PREFIX \
         -DLLVM_ENABLE_PROJECTS=clang \
         -DLLVM_TARGETS_TO_BUILD=X86 \
         -DLLVM_INSTALL_TOOLCHAIN_ONLY=ON \
         -DLLVM_PARALLEL_LINK_JOBS=1 \
         -DCLANG_DEFAULT_OMPSS2_RUNTIME=libnanos6 \
         -DCLANG_DEFAULT_NANOS6_HOME=$INSTALLATION_PREFIX \
         -DCLANG_DEFAULT_NODES_HOME=$INSTALLATION_PREFIX \
    && ninja -C build -j${JOBS} install \
    && cd /home && rm -rf llvm-${LLVM_RELEASE}

# --- Nanos6 (second pass, now aware of the OmpSs-2 clang) ----------------
RUN cd nanos6-version-${NANOS6_VERSION} \
    && ./configure --prefix=$INSTALLATION_PREFIX --with-boost \
         --with-ovni=$OVNI_HOME --with-nanos6-clang=$INSTALLATION_PREFIX \
         --with-hwloc=embedded --with-jemalloc=embedded \
    && make -j${JOBS} && make install \
    && cd /home && rm -rf nanos6-version-${NANOS6_VERSION}

# --- Benchmark sources ----------------------------------------------------
# Compile inside the container on the target machine, since -march=native is CPU-specific:
#   clang++ -O3 -march=native -ffast-math -fompss-2 jacobi_ompss.cpp -lnuma -o jacobi_ompss
COPY src/ /root/test/

WORKDIR /root
CMD ["/bin/bash"]
