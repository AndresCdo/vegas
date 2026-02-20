FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
WORKDIR /src

# Install system dependencies (C++ build tools and HDF5/JSON libraries)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libjsoncpp-dev \
    libhdf5-dev \
    libhdf5-cpp-103-1t64 \
    python3 \
    python3-pip \
    python3-venv \
    git \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Install Python dependencies for analyzers
COPY requirements.txt /src/requirements.txt
RUN pip3 install --no-cache-dir -r /src/requirements.txt

# Copy source code
COPY . /src

# Build with system packages
RUN mkdir -p /src/build \
 && cd /src/build \
 && cmake .. -DCMAKE_BUILD_TYPE=Release \
 && make -j$(nproc)

# Install binary
RUN cp /src/build/vegas /bin/vegas && chmod +x /bin/vegas

# Optional: set entrypoint to vegas
ENTRYPOINT ["/bin/vegas"]
