FROM python:3.11-slim

# hadolint ignore=DL3013
RUN pip install --no-cache-dir --upgrade pip

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        build-essential=12.9 \
        zlib1g-dev=1:1.2.13.dfsg-1 && \
    rm -rf /var/lib/apt/lists/*

# hadolint ignore=DL3042
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,target=/project,rw \
    cd /project && pip install .

RUN oncodrivefml --help

ENTRYPOINT [ "/usr/local/bin/oncodrivefml" ]
