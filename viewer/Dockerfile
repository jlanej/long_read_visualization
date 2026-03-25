# viewer/Dockerfile – Multi-stage build for the long_read_viewer Rust binary.
#
# Build context must be the repository root so that both viewer/ source and
# resources/ (toy dataset) are accessible:
#
#   docker build -f viewer/Dockerfile -t long_read_viewer .
#
# Stages:
#   builder  – compiles the release binary
#   tester   – runs `cargo test` with the toy dataset in place
#   runtime  – minimal image containing only the viewer binary

# ─── Stage 1: Compile ────────────────────────────────────────────────────────
FROM rust:1.85-bookworm AS builder

# eframe / egui / winit require these native packages on Linux at link time
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxcb-render0-dev \
        libxcb-shape0-dev \
        libxcb-xfixes0-dev \
        libxkbcommon-dev \
        libwayland-dev \
        libssl-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app/viewer

# Cache dependency compilation before copying real source.
# Build a stub main so Cargo can resolve and compile all deps.
COPY viewer/Cargo.toml viewer/Cargo.lock ./
RUN mkdir src \
    && echo 'fn main() {}' > src/main.rs \
    && cargo build --release \
    && rm -rf src

# Build the real application
COPY viewer/src/ ./src/
RUN touch src/main.rs && cargo build --release

# ─── Stage 2: Test ───────────────────────────────────────────────────────────
FROM builder AS tester

# The test suite references ../resources/toy_dataset/ relative to the viewer/
# working directory (/app/viewer), so we need the resources one level up.
COPY resources/ /app/resources/

RUN cargo test

# ─── Stage 3: Minimal runtime image ──────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

# Runtime shared libraries required by eframe/winit on Linux
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxcb-render0 \
        libxcb-shape0 \
        libxcb-xfixes0 \
        libxkbcommon0 \
        libwayland-client0 \
        libssl3 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/viewer/target/release/long_read_viewer \
                    /usr/local/bin/long_read_viewer

ENTRYPOINT ["long_read_viewer"]
