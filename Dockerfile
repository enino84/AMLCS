FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    TZ=America/Bogota \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        make \
        git \
        wget \
        curl \
        unzip \
        ca-certificates \
        pkg-config \
        libnetcdf-dev \
        libnetcdff-dev \
        netcdf-bin \
        libhdf5-dev \
        libopenblas-dev \
        liblapack-dev \
        python3 \
        python3-pip \
        python3-dev \
    && rm -rf /var/lib/apt/lists/*

RUN ln -sf /usr/bin/python3 /usr/local/bin/python

RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel \
 && pip3 install --no-cache-dir \
        "numpy<2.0" \
        scipy \
        pandas \
        scikit-learn \
        matplotlib \
        seaborn \
        netCDF4

ENV PROJECT_ROOT=/opt/Research_SPEEDY
WORKDIR ${PROJECT_ROOT}

RUN printf '#!/usr/bin/env bash\nset -e\ncd "${PROJECT_ROOT}"\nfor d in amlcs models/speedy/t21 to_run; do\n  if [ ! -d "${PROJECT_ROOT}/${d}" ]; then\n    echo "[amlcs] WARNING: ${PROJECT_ROOT}/${d} is missing."\n  fi\ndone\nif [ -d "${PROJECT_ROOT}/to_run" ] && [ -d "${PROJECT_ROOT}/amlcs" ]; then\n  cp -ru "${PROJECT_ROOT}/to_run/." "${PROJECT_ROOT}/amlcs/"\nfi\ncd "${PROJECT_ROOT}/amlcs"\nexec "$@"\n' > /usr/local/bin/amlcs-entrypoint.sh \
 && chmod +x /usr/local/bin/amlcs-entrypoint.sh

ENTRYPOINT ["/usr/local/bin/amlcs-entrypoint.sh"]
CMD ["bash", "-c", "echo '=== Pre-procesado ===' && python3 amlcs_pre.py amlcs_pre_t21.csv && echo '=== LEnKF asimilacion ===' && python3 amlcs_da.py amlcs_da_t21_2.csv && echo '=== Listo ==='"]
