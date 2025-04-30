# Dockerfile for simulation pipeline with miniconda, SeqAn, and Sandy

FROM ubuntu:22.04

LABEL maintainer="Cheng-Yuan Liu <danielliu@actgenomics.com>"
ENV RELEASE_DATE=2025.04.28
ENV VERSION=1.0.0
ENV AUTHOR="Cheng-Yuan Liu"

ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/envs/simulation/bin:/opt/conda/bin:$PATH

# Copy necessary files
COPY environment.yml /
COPY Simulator.pm ./

# Update apt source and install base packages
RUN sed -i 's|http://archive.ubuntu.com|http://free.nchc.org.tw/ubuntu|g' /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        wget bzip2 git curl grep sed dpkg make cmake debianutils \
        ca-certificates locales less unzip vim software-properties-common gnupg lsb-release \
        perl zlib1g-dev gcc make cpanminus && \
    locale-gen en_US.UTF-8 && \
    rm -rf /var/lib/apt/lists/*

# Install GCC 12 and G++ 12
RUN apt-get update && \
    add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y gcc-12 g++-12 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100

# Install Clang/LLVM 17
RUN wget https://apt.llvm.org/llvm-snapshot.gpg.key && \
    apt-key add llvm-snapshot.gpg.key && \
    echo "deb http://apt.llvm.org/jammy/ llvm-toolchain-jammy-17 main" > /etc/apt/sources.list.d/llvm.list && \
    apt-get update && \
    apt-get install -y clang-17 lldb-17 lld-17 && \
    update-alternatives --install /usr/bin/clang clang /usr/bin/clang-17 100 && \
    update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-17 100

# Install Miniconda and set up environment
RUN wget --no-check-certificate --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    /opt/conda/bin/conda env update -f /environment.yml && \
    /opt/conda/bin/conda clean --all && \
    rm -f /tmp/miniconda.sh /environment.yml

# Install Sandy using cpanm
RUN cpanm --notest App::Sandy

# Build SeqAn
WORKDIR /tmp
RUN git clone https://github.com/seqan/seqan.git && \
    cd seqan && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    cd / 

RUN cp /tmp/seqan/build/bin/mason_* /opt/conda/envs/simulation/bin && \
    cp /Simulator.pm /usr/local/lib/x86_64-linux-gnu/perl/5.34.0/App/Sandy && \
    rm -rf /tmp/seqan /environment.yml /Simulator.pm

# Default command
CMD ["/bin/bash"]
