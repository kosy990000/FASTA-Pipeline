#!/bin/bash
# ============================================
# RNA-seq 전체 파이프라인 실행 스크립트
# 다운로드 -> QC/Trimming -> Alignment -> Count 추출
#
# 사용법: ./run_pipeline.sh <SRR_ID> [SRR_ID2] ...
# 예시: ./run_pipeline.sh SRR27543183 SRR27543184
# ============================================

set -euo pipefail

# 스크립트 디렉토리
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

# 사용법 출력
usage() {
    echo "============================================"
    echo "RNA-seq 분석 파이프라인"
    echo "============================================"
    echo ""
    echo "사용법: $0 <SRR_ID> [SRR_ID2] ..."
    echo ""
    echo "옵션:"
    echo "  --skip-download    다운로드 단계 건너뛰기"
    echo "  --skip-qc          QC/Trimming 단계 건너뛰기"
    echo "  --skip-align       Alignment 단계 건너뛰기"
    echo "  --merge            마지막에 count 파일 병합"
    echo ""
    echo "예시:"
    echo "  $0 SRR27543183"
    echo "  $0 SRR27543183 SRR27543184 --merge"
    echo "  $0 SRR27543183 --skip-download"
    exit 1
}

# ============================================
# 인자 파싱
# ============================================
SKIP_DOWNLOAD=false
SKIP_QC=false
SKIP_ALIGN=false
DO_MERGE=false
SRR_LIST=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
            ;;
        --skip-qc)
            SKIP_QC=true
            shift
            ;;
        --skip-align)
            SKIP_ALIGN=true
            shift
            ;;
        --merge)
            DO_MERGE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            SRR_LIST+=("$1")
            shift
            ;;
    esac
done

if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
    usage
fi

# ============================================
# 파이프라인 실행
# ============================================
log_info "============================================"
log_info "RNA-seq 파이프라인 시작"
log_info "샘플 수: ${#SRR_LIST[@]}"
log_info "샘플 목록: ${SRR_LIST[*]}"
log_info "============================================"

# 필요한 디렉토리 초기화
init_dirs

# ----------------------------------------
# 1단계: 다운로드
# ----------------------------------------
if [[ "${SKIP_DOWNLOAD}" == "false" ]]; then
    log_info ""
    log_info "========== [1/4] 다운로드 단계 =========="
    for SRR in "${SRR_LIST[@]}"; do
        bash "${SCRIPT_DIR}/01_download.sh" "${SRR}"
    done
else
    log_info "[건너뜀] 다운로드 단계"
fi

# ----------------------------------------
# 2단계: QC 및 Trimming
# ----------------------------------------
if [[ "${SKIP_QC}" == "false" ]]; then
    log_info ""
    log_info "========== [2/4] QC/Trimming 단계 =========="
    for SRR in "${SRR_LIST[@]}"; do
        bash "${SCRIPT_DIR}/02_qc_trim.sh" "${SRR}"
    done
else
    log_info "[건너뜀] QC/Trimming 단계"
fi

# ----------------------------------------
# 3단계: STAR Alignment
# ----------------------------------------
if [[ "${SKIP_ALIGN}" == "false" ]]; then
    log_info ""
    log_info "========== [3/4] Alignment 단계 =========="
    for SRR in "${SRR_LIST[@]}"; do
        bash "${SCRIPT_DIR}/03_align.sh" "${SRR}"
    done
else
    log_info "[건너뜀] Alignment 단계"
fi

# ----------------------------------------
# 4단계: Count 추출
# ----------------------------------------
log_info ""
log_info "========== [4/4] Count 추출 단계 =========="

# conda 환경 활성화 (Python 스크립트용)
# mamba activate rna-align 또는 conda activate rna-align

for SRR in "${SRR_LIST[@]}"; do
    python3 "${SCRIPT_DIR}/04_extract_counts.py" --strand "${STRAND_TYPE}" "${SRR}"
done

# ----------------------------------------
# (선택) Count 파일 병합
# ----------------------------------------
if [[ "${DO_MERGE}" == "true" ]]; then
    log_info ""
    log_info "========== Count 파일 병합 =========="
    python3 "${SCRIPT_DIR}/04_extract_counts.py" --merge
fi

# ============================================
# 완료 메시지
# ============================================
log_info ""
log_info "============================================"
log_info "파이프라인 완료!"
log_info "============================================"
log_info ""
log_info "출력 파일 위치:"
log_info "  - BAM 파일: ${ALIGNED_DIR}/"
log_info "  - Count 파일: ${COUNTS_DIR}/"
log_info "  - QC 리포트: ${QC_DIR}/"
log_info "  - 로그 파일: ${LOGS_DIR}/"
log_info ""
