#!/bin/bash
# ============================================
# SRR 데이터 다운로드 스크립트
# NCBI SRA에서 FASTQ 파일을 다운로드하고 압축
#
# 사용법: ./01_download.sh <SRR_ID> [SRR_ID2] ...
# 예시: ./01_download.sh SRR27543183 SRR27543184
# ============================================

set -euo pipefail

# 스크립트 디렉토리 기준으로 설정 파일 로드
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

# 사용법 출력 함수
usage() {
    echo "사용법: $0 <SRR_ID> [SRR_ID2] ..."
    echo "예시: $0 SRR27543183 SRR27543184"
    exit 1
}

# SRR 다운로드 및 처리 함수
download_srr() {
    local SRR="$1"
    local SRR_DIR="${RAW_DIR}/${SRR}"

    log_info "${SRR} 처리 시작..."

    # 디렉토리 생성 및 이동
    mkdir -p "${SRR_DIR}"
    cd "${SRR_DIR}"

    # ----------------------------------------
    # 1단계: SRA 파일 다운로드 (prefetch)
    # ----------------------------------------
    log_info "[1/4] ${SRR} SRA 파일 다운로드 중..."
    if prefetch "${SRR}" 2>&1 | tee "${LOGS_DIR}/${SRR}_prefetch.log"; then
        log_info "  - ${SRR} 다운로드 성공"
    else
        log_error "  - ${SRR} 다운로드 실패! 네트워크를 확인하세요."
        return 1  # 실패 시 다음 단계로 가지 않고 함수 종료
    fi
    # ----------------------------------------
    # 2단계: FASTQ 변환 (fasterq-dump)
    # paired-end 데이터는 _1, _2로 분리됨
    # ----------------------------------------
    log_info "[2/4] FASTQ 형식으로 변환 중..."
    fasterq-dump "${SRR}" --split-files --threads "${THREADS}" \
        2>&1 | tee "${LOGS_DIR}/${SRR}_fasterq.log"

    # ----------------------------------------
    # 3단계: FASTQ 파일 압축 (gzip)
    # 용량 절약을 위해 압축
    # single-end와 paired-end 모두 처리
    # ----------------------------------------
    log_info "[3/4] FASTQ 파일 압축 중..."
    for fq in "${SRR}"*.fastq; do
        if [[ -f "$fq" ]]; then
            gzip -f "$fq"
            log_info "  - $(basename "$fq") 압축 완료"
        fi
    done

    # ----------------------------------------
    # 4단계: 불필요한 SRA 캐시 파일 삭제
    # prefetch가 다운로드한 .sra 파일 정리
    # ----------------------------------------
    log_info "[4/4] SRA 캐시 파일 정리 중..."
    cleanup_download "${SRR}"

    log_info "${SRR} 다운로드 완료!"
}

# 다운로드 후 정리 함수
cleanup_download() {
    local SRR="$1"

    # SRA 캐시 디렉토리 삭제 (보통 ~/ncbi/public/sra/ 또는 현재 디렉토리)
    if [[ -d "${SRR}" ]]; then
        rm -rf "${SRR}"
        log_info "  - ${SRR}/ 디렉토리 삭제됨"
    fi

    # .sra 파일 삭제
    if [[ -f "${SRR}.sra" ]]; then
        rm -f "${SRR}.sra"
        log_info "  - ${SRR}.sra 파일 삭제됨"
    fi

    # 홈 디렉토리의 ncbi 캐시도 정리
    local NCBI_CACHE="${HOME}/ncbi/public/sra/${SRR}.sra"
    if [[ -f "${NCBI_CACHE}" ]]; then
        rm -f "${NCBI_CACHE}"
        log_info "  - NCBI 캐시 삭제됨"
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

# 각 SRR ID에 대해 다운로드 실행
for SRR in "$@"; do
    download_srr "${SRR}"
done

log_info "모든 다운로드 완료!"
