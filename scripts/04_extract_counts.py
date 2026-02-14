#!/usr/bin/env python3
"""
Gene Count 추출 스크립트
STAR의 ReadsPerGene.out.tab 파일에서 gene count를 추출

사용법:
    python 04_extract_counts.py <SRR_ID> [SRR_ID2] ...
    python 04_extract_counts.py --merge  # 모든 count 파일 병합

예시:
    python 04_extract_counts.py SRR27543183
    python 04_extract_counts.py --merge
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

# ============================================
# 프로젝트 경로 설정
# 환경변수가 있으면 사용 (GSE별 하위 폴더 지원)
# ============================================
PROJECT_DIR = Path(__file__).parent.parent
ALIGNED_DIR = Path(os.environ.get("ALIGNED_DIR", PROJECT_DIR / "aligned"))
COUNTS_DIR = Path(os.environ.get("COUNTS_DIR", PROJECT_DIR / "counts"))


def extract_counts(srr_id: str, strand: str = "unstranded") -> pd.Series:
    """
    STAR ReadsPerGene 파일에서 gene count 추출

    STAR ReadsPerGene.out.tab 파일 형식:
        - 첫 4줄: 메타데이터 (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
        - 이후: gene_id, unstranded, strand1(forward), strand2(reverse)
    """
    input_file = ALIGNED_DIR / srr_id / f"{srr_id}_ReadsPerGene.out.tab"

    if not input_file.exists():
        print(f"[에러] 파일을 찾을 수 없음: {input_file}", file=sys.stderr)
        return None

    # STAR 출력 파일 읽기
    df = pd.read_csv(input_file, sep="\t", header=None)

    # 첫 4줄 메타데이터 제거
    df = df.iloc[4:]
    df.columns = ["gene", "unstranded", "strand1", "strand2"]

    # 스트랜드 선택
    strand_map = {"unstranded": "unstranded", "forward": "strand1", "reverse": "strand2"}
    col = strand_map.get(strand, "unstranded")

    # count 추출
    counts = pd.to_numeric(df[col], errors="coerce")
    counts.index = df["gene"].values
    counts.name = srr_id

    return counts


def save_counts(srr_id: str, counts: pd.Series) -> Path:
    """개별 샘플의 count를 파일로 저장"""
    COUNTS_DIR.mkdir(parents=True, exist_ok=True)
    output_file = COUNTS_DIR / f"{srr_id}.counts.txt"
    counts.to_csv(output_file, sep="\t", header=True)
    print(f"[저장됨] {output_file}")
    return output_file


def merge_all_counts() -> pd.DataFrame:
    """모든 count 파일을 하나의 matrix로 병합 (DESeq2용)"""
    count_files = [f for f in COUNTS_DIR.glob("*.counts.txt") if "merged" not in f.name]

    if not count_files:
        print("[에러] count 파일을 찾을 수 없음", file=sys.stderr)
        return None

    print(f"[정보] {len(count_files)}개 파일 병합 중...")

    dfs = [pd.read_csv(f, sep="\t", index_col=0) for f in sorted(count_files)]
    merged = pd.concat(dfs, axis=1)

    output_file = COUNTS_DIR / "merged_counts.txt"
    merged.to_csv(output_file, sep="\t")
    print(f"[저장됨] {output_file} ({merged.shape[0]} genes x {merged.shape[1]} samples)")

    return merged


def main():
    parser = argparse.ArgumentParser(description="STAR 출력에서 gene count 추출")
    parser.add_argument("srr_ids", nargs="*", help="SRR ID 목록")
    parser.add_argument("--merge", action="store_true", help="count 파일 병합")
    parser.add_argument("--strand", default="unstranded",
                       choices=["unstranded", "forward", "reverse"])

    args = parser.parse_args()

    if args.merge:
        merge_all_counts()
        return

    if not args.srr_ids:
        parser.print_help()
        sys.exit(1)

    for srr_id in args.srr_ids:
        print(f"[처리 중] {srr_id}...")
        counts = extract_counts(srr_id, args.strand)
        if counts is not None:
            save_counts(srr_id, counts)


if __name__ == "__main__":
    main()
