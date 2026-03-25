# ===========================================================================
# Stage 1 — Build the Rust native viewer binary
# ===========================================================================
FROM rust:latest AS rust-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev \
        libxkbcommon-dev libgtk-3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY viewer/ ./viewer/
RUN cargo build --manifest-path viewer/Cargo.toml --release

# ===========================================================================
# Stage 2 — Final image with pipeline tools + Rust binary
# ===========================================================================
FROM ubuntu:22.04

LABEL maintainer="long_read_visualization"
LABEL description="Long-read pre-processing pipeline for multi-panel IGV.js visualization"

ARG MINIMAP2_VERSION=2.28
ARG SAMTOOLS_VERSION=1.21
ARG HTSLIB_VERSION=1.21

ENV DEBIAN_FRONTEND=noninteractive

# ── System dependencies ─────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        libbz2-dev \
        libcurl4-openssl-dev \
        libdeflate-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        pigz \
        python3 \
        python3-dev \
        python3-pip \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ── minimap2 ────────────────────────────────────────────────────────────────
RUN curl -fsSL \
        "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2" \
    | tar -xjf - -C /opt \
    && ln -s "/opt/minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2" /usr/local/bin/minimap2

# ── htslib (bgzip / tabix) ─────────────────────────────────────────────────
RUN curl -fsSL \
        "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
    | tar -xjf - -C /tmp \
    && cd "/tmp/htslib-${HTSLIB_VERSION}" \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && ldconfig \
    && rm -rf "/tmp/htslib-${HTSLIB_VERSION}"

# ── samtools ────────────────────────────────────────────────────────────────
RUN curl -fsSL \
        "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    | tar -xjf - -C /tmp \
    && cd "/tmp/samtools-${SAMTOOLS_VERSION}" \
    && ./configure --prefix=/usr/local --with-htslib=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && rm -rf "/tmp/samtools-${SAMTOOLS_VERSION}"

# ── Python dependencies ─────────────────────────────────────────────────────
# pysam is used by src/assign_haplotypes.py to read/write BAM/CRAM files
# and inject HP:i: haplotype phase tags.  Build against the htslib installed
# above so that pysam uses the same version and CRAM support is available.
RUN HTSLIB_LIBRARY_DIR=/usr/local/lib \
    HTSLIB_INCLUDE_DIR=/usr/local/include \
    pip3 install --no-cache-dir pysam

# ── Rust binary from builder stage ──────────────────────────────────────────
COPY --from=rust-builder /build/viewer/target/release/long_read_viewer /usr/local/bin/long_read_viewer

# ── Pipeline scripts ────────────────────────────────────────────────────────
COPY src/ /opt/long_read_visualization/src/
COPY scripts/ /opt/long_read_visualization/scripts/
COPY server/ /opt/long_read_visualization/server/
COPY resources/ /opt/long_read_visualization/resources/
RUN chmod +x /opt/long_read_visualization/scripts/*.sh

# Download igv.js for the visualization server
# Note: keep in sync with IGV_JS_VERSION in server/app.py
ARG IGV_JS_VERSION=3.1.3
RUN curl -fsSL \
        "https://cdn.jsdelivr.net/npm/igv@${IGV_JS_VERSION}/dist/igv.min.js" \
        -o /opt/long_read_visualization/server/static/igv.min.js

ENV PATH="/opt/long_read_visualization/scripts:${PATH}"

ENTRYPOINT ["preprocess.sh"]
