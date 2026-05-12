# syntax=docker/dockerfile:1

FROM mambaorg/micromamba:1.5.10-jammy

LABEL org.opencontainers.image.title="PanResistome"
LABEL org.opencontainers.image.description="Experimental PanResistome all-in-one workflow image"
LABEL org.opencontainers.image.source="https://github.com/Tasnimul-Arabi-Anik/PanResistome"

USER root
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        coreutils \
        curl \
        findutils \
        git \
        grep \
        gzip \
        jq \
        procps \
        sed \
        tar \
        unzip \
        wget && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/panresistome
COPY envs/ /opt/panresistome/envs/

# Keep environments separate to avoid Python/tool conflicts while making their
# command-line tools visible when Nextflow runs with conda.enabled=false.
# GTDB-Tk is intentionally excluded because its database and runtime footprint
# are too large for the first production image target.
RUN micromamba env create -y -f envs/panr2_comprehensive.yaml && \
    micromamba env create -y -f envs/fetchm.yaml && \
    micromamba env create -y -f envs/abricate.yaml && \
    micromamba env create -y -f envs/amrfinderplus.yaml && \
    micromamba env create -y -f envs/checkm2.yaml && \
    micromamba env create -y -f envs/ani.yaml && \
    micromamba env create -y -f envs/quast.yaml && \
    micromamba env create -y -f envs/mash.yaml && \
    micromamba env create -y -f envs/mobsuite.yaml && \
    micromamba env create -y -f envs/genomad.yaml && \
    micromamba env create -y -f envs/organism_typing.yaml && \
    micromamba clean -a -y

ENV PATH="/opt/conda/envs/panr2_comprehensive_env/bin:/opt/conda/envs/fetchm2_env/bin:/opt/conda/envs/abricate_env/bin:/opt/conda/envs/amrfinderplus_env/bin:/opt/conda/envs/checkm2_env/bin:/opt/conda/envs/ani_env/bin:/opt/conda/envs/quast_env/bin:/opt/conda/envs/mash_env/bin:/opt/conda/envs/mobsuite_env/bin:/opt/conda/envs/genomad_env/bin:/opt/conda/envs/organism_typing_env/bin:/opt/conda/bin:${PATH}"
ENV PYTHONNOUSERSITE=1

COPY . /opt/panresistome

RUN python --version && \
    for tool in \
        abricate \
        mlst \
        integron_finder \
        amrfinder \
        checkm2 \
        fastANI \
        skani \
        quast.py \
        mash \
        mob_recon \
        genomad \
        kleborate \
        kaptive \
        ectyper \
        panr; do \
        command -v "${tool}"; \
    done

WORKDIR /work
CMD ["bash"]
