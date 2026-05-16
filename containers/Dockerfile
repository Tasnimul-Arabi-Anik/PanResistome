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

# Keep environments separate to avoid Python/tool conflicts while making their
# command-line tools visible when Nextflow runs with conda.enabled=false.
# GTDB-Tk is intentionally excluded because its database and runtime footprint
# are too large for the first production image target.
COPY envs/panr2_container.yaml /opt/panresistome/envs/panr2_container.yaml
RUN micromamba env create -y -f envs/panr2_container.yaml && micromamba clean -a -y

COPY envs/mobileelementfinder.yaml /opt/panresistome/envs/mobileelementfinder.yaml
RUN micromamba env create -y -f envs/mobileelementfinder.yaml && micromamba clean -a -y

COPY envs/fetchm.yaml /opt/panresistome/envs/fetchm.yaml
RUN micromamba env create -y -f envs/fetchm.yaml && micromamba clean -a -y

COPY envs/abricate.yaml /opt/panresistome/envs/abricate.yaml
RUN micromamba env create -y -f envs/abricate.yaml && micromamba clean -a -y

COPY envs/amrfinderplus.yaml /opt/panresistome/envs/amrfinderplus.yaml
RUN micromamba env create -y -f envs/amrfinderplus.yaml && micromamba clean -a -y

COPY envs/checkm2.yaml /opt/panresistome/envs/checkm2.yaml
RUN micromamba env create -y -f envs/checkm2.yaml && micromamba clean -a -y

COPY envs/ani.yaml /opt/panresistome/envs/ani.yaml
RUN micromamba env create -y -f envs/ani.yaml && micromamba clean -a -y

COPY envs/quast.yaml /opt/panresistome/envs/quast.yaml
RUN micromamba env create -y -f envs/quast.yaml && micromamba clean -a -y

COPY envs/mash.yaml /opt/panresistome/envs/mash.yaml
RUN micromamba env create -y -f envs/mash.yaml && micromamba clean -a -y

COPY envs/mobsuite.yaml /opt/panresistome/envs/mobsuite.yaml
RUN micromamba env create -y -f envs/mobsuite.yaml && micromamba clean -a -y

COPY envs/genomad.yaml /opt/panresistome/envs/genomad.yaml
RUN micromamba env create -y -f envs/genomad.yaml && micromamba clean -a -y

COPY envs/organism_typing.yaml /opt/panresistome/envs/organism_typing.yaml
RUN micromamba env create -y -f envs/organism_typing.yaml && micromamba clean -a -y

ENV PATH="/opt/conda/envs/panr2_container_env/bin:/opt/conda/envs/mobileelementfinder_env/bin:/opt/conda/envs/fetchm2_env/bin:/opt/conda/envs/abricate_env/bin:/opt/conda/envs/amrfinderplus_env/bin:/opt/conda/envs/checkm2_env/bin:/opt/conda/envs/ani_env/bin:/opt/conda/envs/quast_env/bin:/opt/conda/envs/mash_env/bin:/opt/conda/envs/mobsuite_env/bin:/opt/conda/envs/genomad_env/bin:/opt/conda/envs/organism_typing_env/bin:/opt/conda/bin:${PATH}"
ENV PYTHONNOUSERSITE=1

COPY . /opt/panresistome

RUN python --version && \
    for tool in \
        abricate \
        mlst \
        integron_finder \
        mefinder \
        fetchm2 \
        fetchM \
        seqkit \
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
        command -v "${tool}" || exit 1; \
    done && \
    /opt/conda/envs/checkm2_env/bin/python -c "from checkm2 import modelProcessing; modelProcessing.modelProcessor(1); print('checkm2_model_load=PASS')"

WORKDIR /work
CMD ["bash"]
