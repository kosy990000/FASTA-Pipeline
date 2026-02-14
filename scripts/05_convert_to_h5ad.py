#!/usr/bin/env python3
"""
counts/GSE{ID}/ 안의 raw count 파일들을 메타데이터와 함께 h5ad로 변환

메타데이터 소스:
  1) SraRunTable.csv  -> SRR-GSM 매핑 + SRA 메타데이터
  2) GSE*_family.soft -> GEO Sample 메타데이터 (title, characteristics 등)

두 파일이 없으면 GEO/SRA에서 자동 다운로드 시도

사용법:
  python 05_convert_to_h5ad.py GSE196908
  python 05_convert_to_h5ad.py GSE196908 --no-download
"""

import argparse
import os
import glob
import gzip
import sys
import urllib.request
import subprocess

import numpy as np
import pandas as pd
import anndata as ad

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ============================================
# 다운로드 함수
# ============================================

def download_soft(gse_id, out_dir):
    """GEO에서 SOFT 파일 다운로드"""
    soft_path = os.path.join(out_dir, f"{gse_id}_family.soft")
    if os.path.exists(soft_path):
        print(f"  SOFT 파일 이미 존재: {soft_path}")
        return soft_path

    url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{gse_id[:-3]}nnn/{gse_id}/soft/{gse_id}_family.soft.gz"
    )
    gz_path = soft_path + ".gz"
    print(f"  SOFT 다운로드: {url}")
    try:
        urllib.request.urlretrieve(url, gz_path)
        with gzip.open(gz_path, 'rb') as f_in, open(soft_path, 'wb') as f_out:
            f_out.write(f_in.read())
        os.remove(gz_path)
        print(f"  저장: {soft_path}")
        return soft_path
    except Exception as e:
        print(f"  SOFT 다운로드 실패: {e}")
        if os.path.exists(gz_path):
            os.remove(gz_path)
        return None


def download_sra_table(gse_id, out_dir):
    """Entrez eutils API로 SraRunTable(runinfo) 다운로드"""
    csv_path = os.path.join(out_dir, "SraRunTable.csv")
    if os.path.exists(csv_path) and os.path.getsize(csv_path) > 0:
        print(f"  SraRunTable 이미 존재: {csv_path}")
        return csv_path

    print(f"  SraRunTable 다운로드 시도 (Entrez API)...")
    try:
        import xml.etree.ElementTree as ET

        # Step 1: GEO -> BioProject ID
        geo_url = (
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
            f"?acc={gse_id}&targ=self&form=text&view=brief"
        )
        with urllib.request.urlopen(geo_url, timeout=30) as resp:
            geo_text = resp.read().decode()

        bioproject = None
        for line in geo_text.split('\n'):
            if 'BioProject' in line and 'PRJNA' in line:
                bioproject = line.split('PRJNA')[-1].split('/')[0].strip()
                bioproject = f"PRJNA{bioproject}"
                break

        if not bioproject:
            raise ValueError("BioProject ID를 찾을 수 없음")
        print(f"  BioProject: {bioproject}")

        # Step 2: esearch로 SRA ID 목록
        search_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=sra&term={bioproject}&retmax=500&usehistory=y"
        )
        with urllib.request.urlopen(search_url, timeout=30) as resp:
            search_xml = resp.read().decode()

        root = ET.fromstring(search_xml)
        webenv = root.findtext('WebEnv')
        query_key = root.findtext('QueryKey')
        count = int(root.findtext('Count', '0'))
        print(f"  SRA 레코드 {count}개 발견")

        if count == 0:
            raise ValueError("SRA 레코드 없음")

        # Step 3: efetch로 runinfo CSV 다운로드
        fetch_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=sra&query_key={query_key}&WebEnv={webenv}"
            f"&rettype=runinfo&retmode=csv&retmax={count}"
        )
        with urllib.request.urlopen(fetch_url, timeout=60) as resp:
            runinfo = resp.read().decode()

        # 헤더 추가 (efetch runinfo에 헤더가 없는 경우)
        lines = [l for l in runinfo.strip().split('\n') if l.strip()]
        if lines and not lines[0].startswith('Run'):
            header = (
                "Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,"
                "avgLength,size_MB,AssemblyName,download_path,Experiment,"
                "SampleName,LibraryStrategy,LibrarySelection,LibrarySource,"
                "LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,"
                "BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,"
                "SampleType,TaxID,ScientificName,SampleName2,g1k_pop_code,"
                "source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,"
                "Affection_Status,Analyte_Type,Histological_Type,Body_Site,"
                "CenterName,Submission,dbgap_study_accession,Consent,"
                "RunHash,ReadHash"
            )
            lines.insert(0, header)

        with open(csv_path, 'w') as f:
            f.write('\n'.join(lines) + '\n')

        print(f"  저장: {csv_path}")
        return csv_path

    except Exception as e:
        print(f"  다운로드 실패: {e}")
        print(f"  수동으로 다운로드하세요:")
        print(f"  https://www.ncbi.nlm.nih.gov/Traces/study/?acc={gse_id}")
        print(f"  -> 'Metadata' 버튼 -> {csv_path} 에 저장")
        return None


# ============================================
# SOFT 파일 파싱
# ============================================

def parse_soft(soft_path):
    """SOFT 파일에서 샘플별 메타데이터 추출 (모든 characteristics 자동 파싱)"""
    samples = {}
    current_gsm = None

    with open(soft_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('^SAMPLE'):
                current_gsm = line.split('=')[1].strip()
                samples[current_gsm] = {'GSM': current_gsm}
            elif current_gsm:
                if line.startswith('!Sample_title'):
                    samples[current_gsm]['title'] = line.split('=', 1)[1].strip()
                elif line.startswith('!Sample_source_name_ch1'):
                    samples[current_gsm]['tissue'] = line.split('=', 1)[1].strip()
                elif line.startswith('!Sample_characteristics_ch1'):
                    char = line.split('=', 1)[1].strip()
                    if ':' in char:
                        key, val = char.split(':', 1)
                        key = key.strip().replace(' ', '_').lower()
                        samples[current_gsm][key] = val.strip()

    return samples


# ============================================
# SraRunTable 파싱
# ============================================

def parse_sra_table(csv_path):
    """SraRunTable.csv -> SRR-GSM 매핑 + 메타데이터"""
    df = pd.read_csv(csv_path)
    print(f"  SraRunTable: {len(df)} rows, 컬럼: {list(df.columns)}")

    # Run 컬럼 찾기
    run_col = None
    for c in df.columns:
        if c.lower().strip() in ('run', 'run_accession'):
            run_col = c
            break
    if run_col is None:
        run_col = df.columns[0]

    # GSM(Sample Name) 컬럼 찾기
    gsm_col = None
    for c in df.columns:
        if c.lower().strip() in ('sample name', 'sample_name', 'samplename',
                                  'samplename2', 'geo_accession'):
            # GSM으로 시작하는 값이 있는지 확인
            if df[c].astype(str).str.startswith('GSM').any():
                gsm_col = c
                break
    if gsm_col is None:
        for c in df.columns:
            if df[c].astype(str).str.startswith('GSM').any():
                gsm_col = c
                break

    if gsm_col is None:
        print(f"  [경고] GSM 컬럼을 찾을 수 없음 -> SRR을 샘플명으로 사용")
        return {}, {}

    srr_to_gsm = dict(zip(df[run_col], df[gsm_col]))

    # 메타데이터 수집 (run, gsm 외 모든 컬럼)
    sra_meta = {}
    skip_cols = {run_col, gsm_col}
    meta_cols = [c for c in df.columns if c not in skip_cols]

    for _, row in df.iterrows():
        gsm = row[gsm_col]
        meta = {'SRR': row[run_col]}
        for c in meta_cols:
            val = row[c]
            if pd.notna(val):
                key = c.replace(' ', '_')
                meta[key] = val
        sra_meta[gsm] = meta

    return srr_to_gsm, sra_meta


# ============================================
# 메인 로직
# ============================================

def main():
    parser = argparse.ArgumentParser(
        description="counts/GSE{ID}/ -> h5ad 변환 (SraRunTable + SOFT 메타데이터)"
    )
    parser.add_argument("gse_id", help="GSE accession (예: GSE196908)")
    parser.add_argument("--counts-dir", default=None, help="counts 기본 디렉토리")
    parser.add_argument("--no-download", action="store_true", help="자동 다운로드 비활성화")
    args = parser.parse_args()

    gse_id = args.gse_id
    counts_base = args.counts_dir or os.path.join(PROJECT_DIR, "counts")
    gse_dir = os.path.join(counts_base, gse_id)

    print(f"{'='*50}")
    print(f" {gse_id} -> h5ad 변환")
    print(f"{'='*50}")

    # ── 1. Count 파일 확인 ──
    if not os.path.isdir(gse_dir):
        print(f"에러: 디렉토리 없음: {gse_dir}")
        sys.exit(1)

    count_files = sorted(glob.glob(os.path.join(gse_dir, "*.counts.txt")))
    if not count_files:
        print(f"에러: count 파일 없음: {gse_dir}/*.counts.txt")
        sys.exit(1)
    print(f"\n[1/5] Count 파일 {len(count_files)}개 발견")

    # ── 2. 메타데이터 파일 준비 (없으면 다운로드) ──
    print(f"\n[2/5] 메타데이터 파일 준비")

    # SOFT 파일
    soft_path = None
    soft_candidates = glob.glob(os.path.join(gse_dir, "*.soft"))
    if soft_candidates:
        soft_path = soft_candidates[0]
        print(f"  SOFT 파일 발견: {soft_path}")
    elif not args.no_download:
        soft_path = download_soft(gse_id, gse_dir)

    # SraRunTable
    sra_path = None
    sra_candidate = os.path.join(gse_dir, "SraRunTable.csv")
    if os.path.exists(sra_candidate) and os.path.getsize(sra_candidate) > 0:
        sra_path = sra_candidate
        print(f"  SraRunTable 발견: {sra_path}")
    elif not args.no_download:
        # 빈 파일이면 삭제 후 다시 다운로드
        if os.path.exists(sra_candidate) and os.path.getsize(sra_candidate) == 0:
            os.remove(sra_candidate)
        sra_path = download_sra_table(gse_id, gse_dir)

    # ── 3. 메타데이터 파싱 ──
    print(f"\n[3/5] 메타데이터 파싱")

    soft_meta = {}
    if soft_path and os.path.exists(soft_path):
        soft_meta = parse_soft(soft_path)
        print(f"  SOFT: {len(soft_meta)}개 샘플")
    else:
        print(f"  SOFT 파일 없음 -> 건너뜀")

    srr_to_gsm = {}
    sra_meta = {}
    if sra_path and os.path.exists(sra_path):
        srr_to_gsm, sra_meta = parse_sra_table(sra_path)
        print(f"  SRR-GSM 매핑: {len(srr_to_gsm)}개")
    else:
        print(f"  SraRunTable 없음 -> SRR을 샘플명으로 사용")

    # ── 4. Count 데이터 로드 ──
    print(f"\n[4/5] Count 데이터 로드")
    count_data = {}
    gene_ids = None

    for f in count_files:
        srr = os.path.basename(f).replace(".counts.txt", "")
        sample_id = srr_to_gsm.get(srr, srr)

        df = pd.read_csv(f, sep="\t", index_col=0)
        if gene_ids is None:
            gene_ids = df.index.tolist()
        count_data[sample_id] = df.iloc[:, 0].values
        print(f"  {srr} -> {sample_id}")

    count_matrix = pd.DataFrame(count_data).T
    count_matrix.columns = gene_ids
    print(f"  Matrix: {count_matrix.shape[0]} samples x {count_matrix.shape[1]} genes")

    # ── 5. 메타데이터 병합 (SOFT + SraRunTable) ──
    meta_list = []
    for sample_id in count_matrix.index:
        row = {'sample_id': sample_id}

        # SOFT 메타데이터 (title, tissue, characteristics)
        s = soft_meta.get(sample_id, {})
        for k, v in s.items():
            row[k] = v

        # SraRunTable 메타데이터 (SOFT에 없는 것만 추가)
        sr = sra_meta.get(sample_id, {})
        for k, v in sr.items():
            if k not in row:
                row[k] = v

        meta_list.append(row)

    meta_df = pd.DataFrame(meta_list)
    meta_df.index = count_matrix.index

    # 필수 메타데이터 컬럼만 유지
    # SOFT 유래: 실험 조건 관련 (동적으로 감지)
    # SraRunTable 유래: 핵심 시퀀싱 정보만
    sra_keep = {'SRR', 'Platform', 'Model', 'LibraryLayout', 'LibraryStrategy',
                'Organism', 'ScientificName', 'BioProject', 'BioSample'}
    # SOFT에서 온 컬럼 (sample_id, GSM, title, tissue + characteristics)
    soft_keys = set()
    for info in soft_meta.values():
        soft_keys.update(info.keys())

    keep_cols = []
    for c in meta_df.columns:
        if c in soft_keys or c in sra_keep or c == 'sample_id':
            keep_cols.append(c)
    meta_df = meta_df[keep_cols]

    # ── 6. h5ad 저장 ──
    print(f"\n[5/5] h5ad 저장")
    adata = ad.AnnData(
        X=count_matrix.values.astype(np.float32),
        obs=meta_df,
        var=pd.DataFrame({"gene_id": gene_ids}, index=gene_ids)
    )

    final_dir = os.path.join(PROJECT_DIR, "final_data")
    os.makedirs(final_dir, exist_ok=True)
    output_file = os.path.join(final_dir, f"{gse_id}.h5ad")
    adata.write(output_file)

    print(f"  저장 완료: {output_file}")
    print(f"  Shape: {adata.shape}")
    print(f"  obs 컬럼: {list(adata.obs.columns)}")
    print(f"\n샘플 메타데이터:")
    print(adata.obs.to_string())


if __name__ == "__main__":
    main()
