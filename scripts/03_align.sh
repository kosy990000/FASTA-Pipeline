#!/bin/bash
# ============================================
# STAR Alignment 스크립트
# RNA-seq 리드를 reference genome에 정렬
#
# 사용법: ./03_align.sh <SRR_ID> [SRR_ID2] ...
# 예시: ./03_align.sh SRR27543183
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

# STAR Alignment 실행 함수
run_alignment() {
    local SRR="$1"

    # 입력 파일 경로 (trimmed FASTQ)
    local R1="${TRIMMED_DIR}/${SRR}/${SRR}_1.trim.fastq.gz"
    local R2="${TRIMMED_DIR}/${SRR}/${SRR}_2.trim.fastq.gz"
    local SE="${TRIMMED_DIR}/${SRR}/${SRR}.trim.fastq.gz"  # single-end 파일

    # 출력 디렉토리
    local OUT_DIR="${ALIGNED_DIR}/${SRR}"

    log_info "${SRR} STAR alignment 시작..."

    # ----------------------------------------
    # Single-end / Paired-end 자동 감지
    # ----------------------------------------
    local IS_PAIRED=false
    if [[ -f "${R1}" ]] && [[ -f "${R2}" ]]; then
        IS_PAIRED=true
        log_info "  - Paired-end 데이터 감지됨"
    elif [[ -f "${R1}" ]]; then
        SE="${R1}"
        log_info "  - Single-end 데이터 감지됨 (${SRR}_1.trim.fastq.gz)"
    elif [[ -f "${SE}" ]]; then
        log_info "  - Single-end 데이터 감지됨 (${SRR}.trim.fastq.gz)"
    else
        log_error "입력 파일을 찾을 수 없음: ${SRR}"
        return 1
    fi

    if [[ ! -d "${STAR_INDEX}" ]]; then
        log_error "STAR 인덱스를 찾을 수 없음: ${STAR_INDEX}"
        return 1
    fi

    mkdir -p "${OUT_DIR}"

    # ----------------------------------------
    # 1단계: STAR Alignment 실행
    # 주요 옵션 설명:
    # --readFilesCommand zcat : gzip 파일 직접 읽기
    # --outSAMtype BAM SortedByCoordinate : 정렬된 BAM 출력
    # --quantMode GeneCounts : gene count 동시 생성
    # ----------------------------------------
    log_info "[1/2] STAR alignment 실행 중..."

    if [[ "${IS_PAIRED}" == true ]]; then
        STAR \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${R1}" "${R2}" \
            --readFilesCommand zcat \
            --runThreadN "${STAR_THREADS}" \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outFileNamePrefix "${OUT_DIR}/${SRR}_" \
            2>&1 | tee "${LOGS_DIR}/${SRR}_star.log"
    else
        STAR \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${SE}" \
            --readFilesCommand zcat \
            --runThreadN "${STAR_THREADS}" \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outFileNamePrefix "${OUT_DIR}/${SRR}_" \
            2>&1 | tee "${LOGS_DIR}/${SRR}_star.log"
    fi

    # ----------------------------------------
    # 2단계: 불필요한 임시 파일 정리
    # trimmed FASTQ 및 STAR 임시 파일 삭제
    # ----------------------------------------
    log_info "[2/2] 임시 파일 정리 중..."
    cleanup_after_align "${SRR}"

    log_info "${SRR} STAR alignment 완료!"
}

# Alignment 후 정리 함수
cleanup_after_align() {
    local SRR="$1"
    local OUT_DIR="${ALIGNED_DIR}/${SRR}"

    # BAM 파일이 정상적으로 생성되었는지 확인
    local BAM_FILE="${OUT_DIR}/${SRR}_Aligned.sortedByCoord.out.bam"
    local COUNTS_FILE="${OUT_DIR}/${SRR}_ReadsPerGene.out.tab"

    if [[ -f "${BAM_FILE}" ]] && [[ -f "${COUNTS_FILE}" ]]; then
        # BAM 파일 크기 확인 (최소 1MB 이상)
        local BAM_SIZE=$(stat -f%z "${BAM_FILE}" 2>/dev/null || stat -c%s "${BAM_FILE}" 2>/dev/null)

        if [[ "${BAM_SIZE}" -gt 1000000 ]]; then
            # trimmed FASTQ 파일 삭제 (용량 절약)
            rm -f "${TRIMMED_DIR}/${SRR}/${SRR}_1.trim.fastq.gz"
            rm -f "${TRIMMED_DIR}/${SRR}/${SRR}_2.trim.fastq.gz"
            rm -f "${TRIMMED_DIR}/${SRR}/${SRR}.trim.fastq.gz"
            log_info "  - trimmed FASTQ 파일 삭제됨 (${SRR})"

            # 빈 trimmed 디렉토리 삭제
            rmdir "${TRIMMED_DIR}/${SRR}" 2>/dev/null || true
        else
            log_error "  - BAM 파일이 너무 작음. trimmed 파일 유지."
        fi
    else
        log_error "  - alignment 출력 파일이 없음. trimmed 파일 유지."
    fi

    # STAR 임시 디렉토리 삭제
    if [[ -d "${OUT_DIR}/${SRR}__STARtmp" ]]; then
        rm -rf "${OUT_DIR}/${SRR}__STARtmp"
        log_info "  - STAR 임시 디렉토리 삭제됨"
    fi

    # 불필요한 STAR 로그 파일 정리 (선택적)
    # Log.out, Log.progress.out 등은 유지하되 SJ.out.tab은 삭제
    # (SJ.out.tab은 splice junction 정보로 보통 필요 없음)
    if [[ -f "${OUT_DIR}/${SRR}_SJ.out.tab" ]]; then
        rm -f "${OUT_DIR}/${SRR}_SJ.out.tab"
        log_info "  - SJ.out.tab 삭제됨"
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

# 각 SRR ID에 대해 alignment 실행
for SRR in "$@"; do
    run_alignment "${SRR}"
done

log_info "모든 alignment 완료!"
