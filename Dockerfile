# FASTA-Pipeline — RNA-seq 처리 환경
# build:  docker build -t fasta-pipeline .
# run:    docker run --rm -it -v $PWD:/work fasta-pipeline
#
# environment.yml(conda env: rna-align)을 그대로 설치한다.
FROM mambaorg/micromamba:1.5.8

# conda 환경 생성 (environment.yml 기반)
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# 이후 RUN/CMD가 conda 환경 안에서 실행되도록
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# 파이프라인 코드 복사
WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ /work/scripts/
COPY --chown=$MAMBA_USER:$MAMBA_USER config/ /work/config/
COPY --chown=$MAMBA_USER:$MAMBA_USER command.txt /work/command.txt
RUN chmod +x /work/scripts/*.sh

# 주요 도구가 설치됐는지 확인 (빌드 시 검증)
RUN STAR --version && fastp --version && fastqc --version && \
    samtools --version | head -1 && which prefetch fasterq-dump && \
    python -c "import anndata, pandas, numpy; print('py deps OK')"

CMD ["/bin/bash"]
