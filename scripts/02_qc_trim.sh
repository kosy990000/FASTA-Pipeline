#!/bin/bash
# ============================================
# QC 및 Trimming 스크립트
# FastQC로 품질 검사, fastp로 어댑터 제거 및 품질 필터링
#
# 사용법: ./02_qc_trim.sh <SRR_ID> [SRR_ID2] ...
# 예시: ./02_qc_trim.sh SRR27543183
# ============================================

set -euo pipefail

# 스크립트 디렉토리 기준으로 설정 파일 로드
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

# 사용법 출력 함수
usage() {
    echo "사용법: $0 <SRR_ID> [SRR_ID2] ..."
    echo "예시: $0 SRR27543183"
    exit 1
}

# QC 및 Trimming 실행 함수
run_qc_trim() {
    local SRR="$1"

    # 입력 파일 경로 (raw FASTQ)
    local R1="${RAW_DIR}/${SRR}/${SRR}_1.fastq.gz"
    local R2="${RAW_DIR}/${SRR}/${SRR}_2.fastq.gz"
    local SE="${RAW_DIR}/${SRR}/${SRR}.fastq.gz"  # single-end 파일

    log_info "${SRR} QC 및 Trimming 시작..."

    # ----------------------------------------
    # Single-end / Paired-end 자동 감지
    # ----------------------------------------
    local IS_PAIRED=false
    if [[ -f "${R1}" ]] && [[ -f "${R2}" ]]; then
        IS_PAIRED=true
        log_info "  - Paired-end 데이터 감지됨"
    elif [[ -f "${R1}" ]]; then
        # R1만 있는 경우 (single-end)
        SE="${R1}"
        log_info "  - Single-end 데이터 감지됨 (${SRR}_1.fastq.gz)"
    elif [[ -f "${SE}" ]]; then
        log_info "  - Single-end 데이터 감지됨 (${SRR}.fastq.gz)"
    else
        log_error "입력 파일을 찾을 수 없음: ${SRR}"
        return 1
    fi

    # ----------------------------------------
    # 1단계: FastQC 품질 검사
    # raw 데이터의 시퀀싱 품질 확인
    # - per base quality score
    # - GC content
    # - adapter contamination 등
    # ----------------------------------------
    log_info "[1/3] FastQC 품질 검사 실행 중..."
    mkdir -p "${QC_DIR}/fastqc"

    if [[ "${IS_PAIRED}" == true ]]; then
        fastqc "${R1}" "${R2}" \
            -o "${QC_DIR}/fastqc" \
            -t "${THREADS}" \
            2>&1 | tee "${LOGS_DIR}/${SRR}_fastqc.log"
    else
        fastqc "${SE}" \
            -o "${QC_DIR}/fastqc" \
            -t "${THREADS}" \
            2>&1 | tee "${LOGS_DIR}/${SRR}_fastqc.log"
    fi

    # ----------------------------------------
    # 2단계: fastp Trimming
    # - 어댑터 자동 감지 및 제거
    # - 저품질 염기 트리밍
    # - 너무 짧은 리드 필터링
    # ----------------------------------------
    log_info "[2/3] fastp trimming 실행 중..."
    mkdir -p "${TRIMMED_DIR}/${SRR}"

    if [[ "${IS_PAIRED}" == true ]]; then
        fastp \
            -i "${R1}" \
            -I "${R2}" \
            -o "${TRIMMED_DIR}/${SRR}/${SRR}_1.trim.fastq.gz" \
            -O "${TRIMMED_DIR}/${SRR}/${SRR}_2.trim.fastq.gz" \
            -w "${THREADS}" \
            -h "${QC_DIR}/${SRR}_fastp.html" \
            -j "${QC_DIR}/${SRR}_fastp.json" \
            2>&1 | tee "${LOGS_DIR}/${SRR}_fastp.log"
    else
        fastp \
            -i "${SE}" \
            -o "${TRIMMED_DIR}/${SRR}/${SRR}.trim.fastq.gz" \
            -w "${THREADS}" \
            -h "${QC_DIR}/${SRR}_fastp.html" \
            -j "${QC_DIR}/${SRR}_fastp.json" \
            2>&1 | tee "${LOGS_DIR}/${SRR}_fastp.log"
    fi

    # ----------------------------------------
    # 3단계: 원본 raw FASTQ 파일 삭제
    # trimmed 파일이 생성되었으므로 용량 절약
    # ----------------------------------------
    log_info "[3/3] raw FASTQ 파일 정리 중..."
    cleanup_after_trim "${SRR}"

    log_info "${SRR} QC 및 Trimming 완료!"
}

# Trimming 후 정리 함수
cleanup_after_trim() {
    local SRR="$1"
    local SRR_RAW_DIR="${RAW_DIR}/${SRR}"

    # trimmed 파일 경로
    local TRIM_R1="${TRIMMED_DIR}/${SRR}/${SRR}_1.trim.fastq.gz"
    local TRIM_R2="${TRIMMED_DIR}/${SRR}/${SRR}_2.trim.fastq.gz"
    local TRIM_SE="${TRIMMED_DIR}/${SRR}/${SRR}.trim.fastq.gz"

    local TRIM_OK=false

    # Paired-end trimmed 파일 확인
    if [[ -f "${TRIM_R1}" ]] && [[ -f "${TRIM_R2}" ]]; then
        if [[ -s "${TRIM_R1}" ]] && [[ -s "${TRIM_R2}" ]]; then
            TRIM_OK=true
            rm -f "${SRR_RAW_DIR}/${SRR}_1.fastq.gz"
            rm -f "${SRR_RAW_DIR}/${SRR}_2.fastq.gz"
            log_info "  - raw FASTQ 파일 삭제됨 (paired-end)"
        fi
    # Single-end trimmed 파일 확인
    elif [[ -f "${TRIM_SE}" ]] && [[ -s "${TRIM_SE}" ]]; then
        TRIM_OK=true
        rm -f "${SRR_RAW_DIR}/${SRR}.fastq.gz"
        rm -f "${SRR_RAW_DIR}/${SRR}_1.fastq.gz"  # single-end가 _1로 저장된 경우
        log_info "  - raw FASTQ 파일 삭제됨 (single-end)"
    fi

    if [[ "${TRIM_OK}" == true ]]; then
        # 빈 디렉토리 삭제
        rmdir "${SRR_RAW_DIR}" 2>/dev/null || true
    else
        log_error "  - trimmed 파일이 없거나 비어있음. raw 파일 유지."
    fi
}

# ============================================
# 메인 실행부
# ============================================
if [[ $# -lt 1 ]]; then
    usage
fi

# 필요한 디렉토리 초기화
init_dirs

# 각 SRR ID에 대해 QC 및 Trimming 실행
for SRR in "$@"; do
    run_qc_trim "${SRR}"
done

log_info "모든 QC 및 Trimming 완료!"
