# FASTA-Pipeline 🧬

NCBI SRA의 RNA-seq raw read(FASTQ)를 받아 **다운로드 → QC/Trimming → STAR 정렬 → gene count → h5ad(AnnData)** 까지 자동 처리하는 bulk RNA-seq 파이프라인.

- **실행 환경**: conda `rna-align` (`environment.yml`)
- **Reference**: mouse mm10 / Ensembl release-95 (GRCm38)
- **입력**: SRR ID 또는 GSE 샘플시트(`config/*.txt`)
- **최종 산출물**: `final_data/GSE*.h5ad` (count matrix + 메타데이터)

---

## 파이프라인 단계

| 단계 | 스크립트 | 도구 | 입력 → 출력 |
|------|----------|------|-------------|
| 1. 다운로드 | `scripts/01_download.sh` | sra-tools (prefetch, fasterq-dump) | SRR ID → `raw/srr/*.fastq.gz` |
| 2. QC/Trim | `scripts/02_qc_trim.sh` | FastQC, fastp | raw FASTQ → `trimmed/*.trim.fastq`, `qc/` 리포트 |
| 3. 정렬 | `scripts/03_align.sh` | STAR | trimmed FASTQ → `aligned/*.bam`, `ReadsPerGene.out.tab` |
| 4. count 추출 | `scripts/04_extract_counts.py` | pandas | STAR count → `counts/` (`--merge`로 병합) |
| 5. h5ad 변환 | `scripts/05_convert_to_h5ad.py` | anndata | counts + GEO/SRA 메타데이터 → `final_data/GSE*.h5ad` |

**오케스트레이터**

| 스크립트 | 역할 |
|----------|------|
| `scripts/run_pipeline.sh <SRR_ID> ...` | SRR 단위로 1→4단계 일괄 실행 |
| `scripts/run_from_gse.sh <GSE.txt>` | GSE 샘플시트 기반으로 샘플별 전체 실행 (권장) |

---

## 설치

### conda
```bash
conda env create -f environment.yml
conda activate rna-align
```

### Docker
```bash
docker build -t fasta-pipeline .
docker run --rm -it -v $PWD:/work fasta-pipeline
```

---

## 사용법

```bash
conda activate rna-align

# (최초 1회) reference 준비 — 게놈/GTF 다운로드 + STAR 인덱스 생성
#   command.txt 의 명령 참고

# A. GSE 샘플시트 기반 실행 (권장)
./scripts/run_from_gse.sh config/GSE196908.txt
./scripts/run_from_gse.sh config/GSE196908.txt --skip-download   # 다운로드 건너뛰기
./scripts/run_from_gse.sh config/GSE196908.txt --dry-run         # 실행 계획만 출력

# B. SRR 단위 실행
./scripts/run_pipeline.sh SRR27543183 SRR27543184

# C. h5ad 변환 (메타데이터 포함)
python scripts/05_convert_to_h5ad.py GSE196908
```

---

## GSE 샘플시트 형식

탭(`\t`) 구분, `#`은 주석. 템플릿: `config/GSE_TEMPLATE.txt`

```
# SRR_ID      SAMPLE_NAME   LAYOUT   STRAND    SPECIES
SRR7722937    Control_1     paired   reverse   mouse
SRR7722938    Control_2     paired   reverse   mouse
SRR7722939    Treatment_1   paired   reverse   mouse
```

| 컬럼 | 설명 |
|------|------|
| SRR_ID | SRA Run Accession (필수) |
| SAMPLE_NAME | 샘플 이름 (공백 불가) |
| LAYOUT | paired / single |
| STRAND | reverse / forward / unstranded |
| SPECIES | mouse / human |

---

## 디렉터리 구조

```
Fastaq/
├── config/
│   ├── config.sh          # 공통 설정 (경로·스레드·strand·유틸 함수)
│   ├── GSE196908.txt       # GSE 샘플시트 예시
│   └── GSE_TEMPLATE.txt    # 샘플시트 템플릿
├── scripts/                # 01~05 단계 + run_pipeline / run_from_gse
├── ref/mm10_ensembl95/     # STAR 인덱스 + 게놈 FASTA + GTF (git 미포함)
├── raw/, trimmed/, qc/, aligned/, counts/, final_data/, logs/   # 산출물 (git 미포함)
├── command.txt             # reference 준비 명령 메모
├── environment.yml         # conda 환경 정의
└── Dockerfile              # 컨테이너 빌드 정의
```

---

## 설정 (`config/config.sh`)

| 항목 | 기본값 |
|------|--------|
| CONDA_ENV | `rna-align` |
| Reference | `ref/mm10_ensembl95/` (GRCm38, Ensembl 95) |
| STRAND_TYPE | `reverse` (Illumina TruSeq Stranded dUTP) |
| THREADS / STAR_THREADS | 8 / 8 |

> 모든 스크립트가 `config.sh`를 `source`해 경로·파라미터를 공유한다. 환경변수로 override 가능(GSE별 하위 폴더 지원). 다른 종/버전을 쓰려면 `REF_DIR`/`GENOME_FA`/`GTF_FILE`를 바꾸고 STAR 인덱스를 재생성한다.

---

## 사용 도구

STAR 2.7.10b · fastp 0.22.0 · FastQC 0.12.1 · samtools 1.6 · subread 2.0.6 · sra-tools 3.4.1 · anndata · pandas · numpy (python 3.11) — 전체 목록은 `environment.yml` 참고.
