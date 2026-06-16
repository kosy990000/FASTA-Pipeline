# 🧬 RNA-seq FASTQ 처리 파이프라인 (Fastaq)

NCBI SRA의 RNA-seq raw read(FASTQ)를 받아
다운로드 → QC/Trimming → STAR 정렬 → gene count → h5ad(AnnData)까지
자동 처리하는 파이프라인.

> 코드 위치: `/data/project/shinyoung/Fastaq` · 실행 환경: conda `rna-align`
> Reference: mouse mm10 / Ensembl release-95 (GRCm38)

---

## 1. 파이프라인 단계

| 단계 | 스크립트 | 도구 | 입력 → 출력 |
|------|----------|------|-------------|
| 1. 다운로드 | `01_download.sh` | sra-tools (prefetch, fasterq-dump) | SRR ID → `raw/srr/*.fastq.gz` |
| 2. QC/Trim | `02_qc_trim.sh` | FastQC, fastp | raw FASTQ → `trimmed/*.trim.fastq`, `qc/` 리포트 |
| 3. 정렬 | `03_align.sh` | STAR | trimmed FASTQ → `aligned/*.bam` + `ReadsPerGene.out.tab` |
| 4. count 추출 | `04_extract_counts.py` | pandas | STAR count → `counts/` (병합 가능 `--merge`) |
| 5. h5ad 변환 | `05_convert_to_h5ad.py` | anndata | counts + GEO/SRA 메타데이터 → `final_data/GSE*.h5ad` |

**오케스트레이터**
| 스크립트 | 역할 |
|----------|------|
| `run_pipeline.sh <SRR_ID> ...` | SRR 단위로 1→4단계 일괄 실행 |
| `run_from_gse.sh <GSE.txt>` | GSE 샘플시트 기반으로 샘플별 전체 실행 |

---

## 2. 디렉터리 구조

```
Fastaq/
├── config/
│   ├── config.sh          # 공통 설정 (경로·스레드·strand·유틸 함수)
│   ├── GSE196908.txt       # GSE 샘플시트 예시 (테스트 1샘플)
│   └── GSE_TEMPLATE.txt    # 샘플시트 템플릿
├── scripts/                # 01~05 단계 + run_pipeline / run_from_gse
├── ref/mm10_ensembl95/     # STAR 인덱스 + 게놈 FASTA + GTF
├── raw/srr/                # 원본 FASTQ
├── trimmed/                # trimming된 FASTQ
├── qc/                     # FastQC / fastp 리포트
├── aligned/                # BAM + STAR 출력
├── counts/                 # gene count
├── final_data/             # 최종 h5ad
├── logs/                   # 단계별 로그
├── command.txt             # reference 준비(게놈 다운로드 + STAR index) 명령 메모
└── environment.yml         # conda 환경 정의 (rna-align)
```

---

## 3. 핵심 설정 (`config/config.sh`)

| 항목 | 값 |
|------|-----|
| CONDA_ENV | `rna-align` |
| Reference | `ref/mm10_ensembl95/` (GRCm38, Ensembl 95) |
| STRAND_TYPE | `reverse` (Illumina TruSeq Stranded dUTP, 가장 흔함) |
| THREADS / STAR_THREADS | 8 / 8 |

> 모든 스크립트가 `config.sh`를 `source`해 경로·파라미터를 공유. 환경변수로 override 가능 (GSE별 하위 폴더 지원).

---

## 4. GSE 샘플시트 형식 (탭 구분, `#`은 주석)

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

## 5. 실행 방법

```bash
conda activate rna-align

# (최초 1회) reference 준비 — 게놈/GTF 다운로드 + STAR 인덱스 생성
#   command.txt 참고

# A. SRR 단위 실행
./scripts/run_pipeline.sh SRR27543183 SRR27543184

# B. GSE 샘플시트 기반 실행 (권장)
./scripts/run_from_gse.sh config/GSE196908.txt
./scripts/run_from_gse.sh config/GSE196908.txt --skip-download
./scripts/run_from_gse.sh config/GSE196908.txt --dry-run

# C. h5ad 변환 (메타데이터 포함)
python scripts/05_convert_to_h5ad.py GSE196908
```

---

## 6. 사용 도구 / 버전 (environment.yml: `rna-align`)

| 도구 | 버전 | 용도 |
|------|------|------|
| python | 3.11 | 04·05 스크립트 |
| STAR | 2.7.10b | 정렬 |
| fastp | 0.22.0 | 어댑터 제거·품질 필터 |
| fastqc | 0.12.1 | 품질 검사 |
| samtools | 1.6 | BAM 처리 |
| subread | 2.0.6 | featureCounts(보조) |
| pandas / numpy | 2.3.3 / 2.4.1 | count 처리 |
| openjdk | 11 | FastQC 실행 |
| sra-tools | 3.4.1 | prefetch / fasterq-dump (다운로드) |
| anndata | 0.12.x | h5ad 저장 |

---

## 7. 주의사항

- ✅ **environment.yml 보강 완료**: 누락돼 있던 `sra-tools`(1단계)·`anndata`(5단계)를 추가함
  → `conda env create -f environment.yml` 한 번으로 전 단계 재현 가능
- reference는 mouse mm10(Ensembl 95)에 고정 — 다른 종/버전은 `config.sh`의 `REF_DIR`/`GENOME_FA`/`GTF_FILE`와 STAR 인덱스 재생성 필요
- strand 기본값 `reverse` — 라이브러리 종류에 맞게 `STRAND_TYPE` 조정
