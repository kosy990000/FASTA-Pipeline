#!/bin/bash
# ============================================
# GSE txt 파일 기반 RNA-seq 파이프라인 실행 스크립트
# GSE txt 파일에서 SRR 목록과 메타정보를 읽어
# 샘플별로 파이프라인을 실행
#
# 사용법: ./run_from_gse.sh <GSE_FILE.txt> [옵션]
# 예시:
#   ./run_from_gse.sh ../config/GSE119340.txt
#   ./run_from_gse.sh ../config/GSE119340.txt --skip-download
#   ./run_from_gse.sh ../config/GSE119340.txt --dry-run
# ============================================

set -euo pipefail

# 스크립트 디렉토리 기준으로 설정 파일 로드
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

# 사용법 출력
usage() {
    echo "============================================"
    echo "GSE 기반 RNA-seq 파이프라인"
    echo "============================================"
    echo ""
    echo "사용법: $0 <GSE_FILE.txt> [옵션]"
    echo ""
    echo "GSE 파일 형식 (탭 구분):"
    echo "  SRR_ID  SAMPLE_NAME  LAYOUT  STRAND  SPECIES"
    echo "  (템플릿: config/GSE_TEMPLATE.txt 참고)"
    echo ""
    echo "옵션:"
    echo "  --skip-download    다운로드 단계 건너뛰기"
    echo "  --skip-qc          QC/Trimming 단계 건너뛰기"
    echo "  --skip-align       Alignment 단계 건너뛰기"
    echo "  --skip-count       Count 추출 단계 건너뛰기"
    echo "  --merge            마지막에 count 파일 병합"
    echo "  --dry-run          실행하지 않고 파싱 결과만 출력"
    echo ""
    echo "예시:"
    echo "  $0 ../config/GSE119340.txt"
    echo "  $0 ../config/GSE119340.txt --skip-download --merge"
    echo "  $0 ../config/GSE119340.txt --dry-run"
    exit 1
}

# ============================================
# SPECIES별 Reference 경로 설정
# ============================================
get_ref_paths() {
    local SPECIES="$1"
    case "${SPECIES}" in
        mouse)
            REF_DIR="${PROJECT_DIR}/ref/mm10_ensembl95"
            STAR_INDEX="${REF_DIR}/STARindex"
            GENOME_FA="${REF_DIR}/Mus_musculus.GRCm38.dna.primary_assembly.fa"
            GTF_FILE="${REF_DIR}/Mus_musculus.GRCm38.95.gtf"
            ;;
        human)
            REF_DIR="${PROJECT_DIR}/ref/hg38_ensembl"
            STAR_INDEX="${REF_DIR}/STARindex"
            GENOME_FA="${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
            GTF_FILE="${REF_DIR}/Homo_sapiens.GRCh38.gtf"
            ;;
        *)
            log_error "지원하지 않는 SPECIES: ${SPECIES}"
            return 1
            ;;
    esac
}

# ============================================
# 인자 파싱
# ============================================
GSE_FILE=""
SKIP_DOWNLOAD=false
SKIP_QC=false
SKIP_ALIGN=false
SKIP_COUNT=false
DO_MERGE=false
DRY_RUN=false

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
        --skip-count)
            SKIP_COUNT=true
            shift
            ;;
        --merge)
            DO_MERGE=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            if [[ -z "${GSE_FILE}" ]]; then
                GSE_FILE="$1"
            else
                log_error "알 수 없는 인자: $1"
                usage
            fi
            shift
            ;;
    esac
done

if [[ -z "${GSE_FILE}" ]]; then
    log_error "GSE 파일을 지정해주세요."
    usage
fi

if [[ ! -f "${GSE_FILE}" ]]; then
    log_error "GSE 파일을 찾을 수 없음: ${GSE_FILE}"
    exit 1
fi

# ============================================
# GSE 파일 파싱
# ============================================
GSE_NAME=$(basename "${GSE_FILE}" .txt)

SRR_IDS=()
SAMPLE_NAMES=()
LAYOUTS=()
STRANDS=()
SPECIES_LIST=()

LINE_NUM=0
while IFS=$'\t' read -r SRR_ID SAMPLE_NAME LAYOUT STRAND SPECIES || [[ -n "${SRR_ID}" ]]; do
    LINE_NUM=$((LINE_NUM + 1))

    # 주석/빈 줄 건너뛰기
    [[ -z "${SRR_ID}" ]] && continue
    [[ "${SRR_ID}" =~ ^# ]] && continue

    # 필수 컬럼 검증
    if [[ -z "${SAMPLE_NAME}" || -z "${LAYOUT}" || -z "${STRAND}" || -z "${SPECIES}" ]]; then
        log_error "Line ${LINE_NUM}: 컬럼이 부족합니다. 탭으로 5개 컬럼을 구분해주세요."
        log_error "  -> ${SRR_ID}"
        exit 1
    fi

    # SRR ID 형식 검증
    if [[ ! "${SRR_ID}" =~ ^[SDE]RR[0-9]+$ ]]; then
        log_error "Line ${LINE_NUM}: 잘못된 SRR ID 형식: ${SRR_ID}"
        exit 1
    fi

    # LAYOUT 검증
    if [[ "${LAYOUT}" != "paired" && "${LAYOUT}" != "single" ]]; then
        log_error "Line ${LINE_NUM}: LAYOUT은 'paired' 또는 'single'이어야 합니다: ${LAYOUT}"
        exit 1
    fi

    # STRAND 검증
    if [[ "${STRAND}" != "reverse" && "${STRAND}" != "forward" && "${STRAND}" != "unstranded" ]]; then
        log_error "Line ${LINE_NUM}: STRAND는 'reverse', 'forward', 'unstranded' 중 하나여야 합니다: ${STRAND}"
        exit 1
    fi

    # SPECIES 검증
    if [[ "${SPECIES}" != "mouse" && "${SPECIES}" != "human" ]]; then
        log_error "Line ${LINE_NUM}: SPECIES는 'mouse' 또는 'human'이어야 합니다: ${SPECIES}"
        exit 1
    fi

    SRR_IDS+=("${SRR_ID}")
    SAMPLE_NAMES+=("${SAMPLE_NAME}")
    LAYOUTS+=("${LAYOUT}")
    STRANDS+=("${STRAND}")
    SPECIES_LIST+=("${SPECIES}")

done < "${GSE_FILE}"

SAMPLE_COUNT=${#SRR_IDS[@]}

if [[ ${SAMPLE_COUNT} -eq 0 ]]; then
    log_error "GSE 파일에서 샘플을 찾을 수 없습니다."
    exit 1
fi

# ============================================
# 파싱 결과 출력
# ============================================
log_info "============================================"
log_info "GSE 파일: ${GSE_FILE}"
log_info "GSE 이름: ${GSE_NAME}"
log_info "샘플 수: ${SAMPLE_COUNT}"
log_info "============================================"

printf "%-4s  %-15s  %-20s  %-8s  %-12s  %-6s\n" \
    "#" "SRR_ID" "SAMPLE_NAME" "LAYOUT" "STRAND" "SPECIES"
printf "%s\n" "----  ---------------  --------------------  --------  ------------  ------"

for i in $(seq 0 $((SAMPLE_COUNT - 1))); do
    printf "%-4d  %-15s  %-20s  %-8s  %-12s  %-6s\n" \
        $((i + 1)) "${SRR_IDS[$i]}" "${SAMPLE_NAMES[$i]}" \
        "${LAYOUTS[$i]}" "${STRANDS[$i]}" "${SPECIES_LIST[$i]}"
done
echo ""

# dry-run이면 여기서 종료
if [[ "${DRY_RUN}" == "true" ]]; then
    log_info "[dry-run] 실제 실행 없이 종료합니다."
    exit 0
fi

# ============================================
# GSE별 하위 디렉토리 설정
# 각 폴더 밑에 GSE 번호 폴더를 만들어 정리
# 예: qc/GSE119340/, logs/GSE119340/ 등
# ============================================
export RAW_DIR="${PROJECT_DIR}/raw/srr/${GSE_NAME}"
export TRIMMED_DIR="${PROJECT_DIR}/trimmed/${GSE_NAME}"
export QC_DIR="${PROJECT_DIR}/qc/${GSE_NAME}"
export ALIGNED_DIR="${PROJECT_DIR}/aligned/${GSE_NAME}"
export COUNTS_DIR="${PROJECT_DIR}/counts/${GSE_NAME}"
export LOGS_DIR="${PROJECT_DIR}/logs/${GSE_NAME}"

GSE_COUNTS_DIR="${PROJECT_DIR}/${GSE_NAME}-count"
mkdir -p "${GSE_COUNTS_DIR}"

# 필요한 디렉토리 초기화 (GSE 하위 폴더 포함)
init_dirs

# ============================================
# 파이프라인 실행 (샘플별)
# ============================================
FAILED_SAMPLES=()

for i in $(seq 0 $((SAMPLE_COUNT - 1))); do
    SRR="${SRR_IDS[$i]}"
    SAMPLE="${SAMPLE_NAMES[$i]}"
    LAYOUT="${LAYOUTS[$i]}"
    STRAND="${STRANDS[$i]}"
    SPECIES="${SPECIES_LIST[$i]}"

    log_info ""
    log_info "=========================================="
    log_info "샘플 [$((i + 1))/${SAMPLE_COUNT}]: ${SRR} (${SAMPLE})"
    log_info "  LAYOUT=${LAYOUT}, STRAND=${STRAND}, SPECIES=${SPECIES}"
    log_info "=========================================="

    # SPECIES에 맞는 reference 경로 설정
    get_ref_paths "${SPECIES}"
    export STAR_INDEX GTF_FILE GENOME_FA

    # --- 1단계: 다운로드 ---
    if [[ "${SKIP_DOWNLOAD}" == "false" ]]; then
        log_info "[1/4] 다운로드..."
        if ! bash "${SCRIPT_DIR}/01_download.sh" "${SRR}"; then
            log_error "${SRR} 다운로드 실패. 건너뜁니다."
            FAILED_SAMPLES+=("${SRR}:download")
            continue
        fi
    else
        log_info "[건너뜀] 다운로드"
    fi

    # --- 2단계: QC/Trimming ---
    if [[ "${SKIP_QC}" == "false" ]]; then
        log_info "[2/4] QC/Trimming..."
        if ! bash "${SCRIPT_DIR}/02_qc_trim.sh" "${SRR}"; then
            log_error "${SRR} QC/Trimming 실패. 건너뜁니다."
            FAILED_SAMPLES+=("${SRR}:qc_trim")
            continue
        fi
    else
        log_info "[건너뜀] QC/Trimming"
    fi

    # --- 3단계: Alignment ---
    if [[ "${SKIP_ALIGN}" == "false" ]]; then
        log_info "[3/4] Alignment..."
        export STRAND_TYPE="${STRAND}"
        if ! bash "${SCRIPT_DIR}/03_align.sh" "${SRR}"; then
            log_error "${SRR} Alignment 실패. 건너뜁니다."
            FAILED_SAMPLES+=("${SRR}:align")
            continue
        fi
    else
        log_info "[건너뜀] Alignment"
    fi

    # --- 4단계: Count 추출 ---
    if [[ "${SKIP_COUNT}" == "false" ]]; then
        log_info "[4/4] Count 추출..."
        if ! python3 "${SCRIPT_DIR}/04_extract_counts.py" --strand "${STRAND}" "${SRR}"; then
            log_error "${SRR} Count 추출 실패."
            FAILED_SAMPLES+=("${SRR}:count")
            continue
        fi
    else
        log_info "[건너뜀] Count 추출"
    fi

    log_info "${SRR} (${SAMPLE}) 완료!"
done

# ============================================
# (선택) Count 파일 병합
# ============================================
if [[ "${DO_MERGE}" == "true" && "${SKIP_COUNT}" == "false" ]]; then
    log_info ""
    log_info "========== Count 파일 병합 =========="
    python3 "${SCRIPT_DIR}/04_extract_counts.py" --merge
fi

# ============================================
# 최종 결과 요약
# ============================================
log_info ""
log_info "============================================"
log_info "파이프라인 완료: ${GSE_NAME}"
log_info "============================================"
log_info "  총 샘플: ${SAMPLE_COUNT}"
log_info "  성공: $((SAMPLE_COUNT - ${#FAILED_SAMPLES[@]}))"
log_info "  실패: ${#FAILED_SAMPLES[@]}"

if [[ ${#FAILED_SAMPLES[@]} -gt 0 ]]; then
    log_info ""
    log_info "실패한 샘플:"
    for fail in "${FAILED_SAMPLES[@]}"; do
        log_info "  - ${fail}"
    done
fi

log_info ""
log_info "출력 파일 위치 (GSE: ${GSE_NAME}):"
log_info "  - Raw FASTQ:  raw/srr/${GSE_NAME}/"
log_info "  - Trimmed:    trimmed/${GSE_NAME}/"
log_info "  - QC 리포트:  qc/${GSE_NAME}/"
log_info "  - BAM 파일:   aligned/${GSE_NAME}/"
log_info "  - Count 파일: counts/${GSE_NAME}/"
log_info "  - 로그 파일:  logs/${GSE_NAME}/"
log_info "  - GSE Count:  ${GSE_NAME}-count/"
