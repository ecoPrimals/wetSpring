#!/usr/bin/env python3
"""
genome_search.py - Human genome search for ecoPrimals/wetSpring.

Search GRCh38 by gene name, gene ID, chromosome region, or sequence.

Usage:
  python3 genome_search.py gene TP53
  python3 genome_search.py gene BRCA1
  python3 genome_search.py gene 7157
  python3 genome_search.py region chr17:7668421-7687490
  python3 genome_search.py sequence ATCGATCGATCG
  python3 genome_search.py stats
"""
import gzip, sys, os, json, hashlib
from pathlib import Path

REF_DIR = Path(__file__).parent.parent / "data" / "reference" / "GRCh38"
FASTA_GZ = REF_DIR / "GRCh38.p14_genomic.fna.gz"
GENE_INFO = REF_DIR / "Homo_sapiens.gene_info.gz"

CHR_MAP = {
    "chr1":"NC_000001.11","chr2":"NC_000002.12","chr3":"NC_000003.12",
    "chr4":"NC_000004.12","chr5":"NC_000005.10","chr6":"NC_000006.12",
    "chr7":"NC_000007.14","chr8":"NC_000008.11","chr9":"NC_000009.12",
    "chr10":"NC_000010.11","chr11":"NC_000011.10","chr12":"NC_000012.12",
    "chr13":"NC_000013.11","chr14":"NC_000014.9","chr15":"NC_000015.10",
    "chr16":"NC_000016.10","chr17":"NC_000017.11","chr18":"NC_000018.10",
    "chr19":"NC_000019.10","chr20":"NC_000020.11","chr21":"NC_000021.9",
    "chr22":"NC_000022.11","chrX":"NC_000023.11","chrY":"NC_000024.10",
    "chrM":"NC_012920.1",
}
ACC_TO_CHR = {v:k for k,v in CHR_MAP.items()}

def search_gene(query):
    if not GENE_INFO.exists(): return []
    results = []
    ql = query.lower()
    is_num = query.isdigit()
    with gzip.open(GENE_INFO,"rt") as f:
        f.readline()
        for line in f:
            fs = line.strip().split("\t")
            if len(fs)<15: continue
            gid,sym,_,_,syns,_,chrom,loc,desc,gtype = fs[1],fs[2],fs[3],fs[3],fs[4],fs[5],fs[6],fs[7],fs[8],fs[9]
            match = False
            if is_num and gid==query: match=True
            elif not is_num:
                if sym.lower()==ql: match=True
                elif ql in syns.lower().split("|"): match=True
                elif len(query)>=3 and ql in sym.lower(): match=True
            if match:
                results.append({"gene_id":gid,"symbol":sym,"synonyms":syns,"chromosome":chrom,"map_location":loc,"description":desc,"type":gtype})
    return results

def extract_region(chrom, start, end):
    if not FASTA_GZ.exists(): return None
    acc = CHR_MAP.get(chrom, chrom)
    seq = []; in_target = False; pos = 0
    with gzip.open(FASTA_GZ,"rt") as f:
        for line in f:
            if line.startswith(">"):
                if in_target: break
                if line.split()[0][1:]==acc: in_target=True
                continue
            if in_target:
                ls = line.strip()
                le = pos+len(ls)
                if le>start and pos<end:
                    s=max(0,start-pos); e=min(len(ls),end-pos)
                    seq.append(ls[s:e])
                pos=le
                if pos>=end: break
    return "".join(seq) if seq else None

def search_sequence(query, max_results=5):
    if not FASTA_GZ.exists(): return []
    qu = query.upper(); k=len(qu); results=[]; cur_chr=""; pos=0; buf=""
    with gzip.open(FASTA_GZ,"rt") as f:
        for line in f:
            if line.startswith(">"):
                acc=line.split()[0][1:]; cur_chr=ACC_TO_CHR.get(acc,acc); pos=0; buf=""; continue
            seq=line.strip().upper()
            buf = buf[-(k-1):]+seq if buf else seq
            ss = max(0,len(buf)-len(seq)-k+1)
            idx=buf.find(qu,ss)
            while idx!=-1:
                ap=pos-(len(buf)-len(seq))+idx
                results.append({"chromosome":cur_chr,"position":ap,"strand":"+"})
                if len(results)>=max_results: return results
                idx=buf.find(qu,idx+1)
            pos+=len(seq)
    return results

def genome_stats():
    if not FASTA_GZ.exists(): print("No reference found"); return
    print("=== GRCh38.p14 Human Reference Genome ===")
    print(f"  File: {FASTA_GZ}")
    print(f"  Size: {FASTA_GZ.stat().st_size/1024/1024:.0f} MB (compressed)")
    nc=0; tb=0; chroms=[]
    with gzip.open(FASTA_GZ,"rt") as f:
        cl=0; cn=""
        for line in f:
            if line.startswith(">"):
                if cn: chroms.append((cn,cl)); tb+=cl
                cn=line.strip()[1:].split()[0]; cl=0; nc+=1
            else: cl+=len(line.strip())
        if cn: chroms.append((cn,cl)); tb+=cl
    print(f"  Sequences: {nc}")
    print(f"  Total: {tb:,} bp (~{tb/1e9:.2f} Gbp)")
    print("\n  Primary chromosomes:")
    for acc,l in chroms[:25]:
        n=ACC_TO_CHR.get(acc,acc)
        print(f"    {n:6s} ({acc:15s}): {l:>13,} bp")
    if GENE_INFO.exists():
        gc=sum(1 for _ in gzip.open(GENE_INFO,"rt"))-1
        print(f"\n  Gene annotations: {gc:,}")
    print(f"\n  Provenance: NCBI RefSeq GCF_000001405.40")
    print(f"  License: Public Domain")

def main():
    if len(sys.argv)<2: print(__doc__); return
    cmd=sys.argv[1].lower()
    if cmd=="gene" and len(sys.argv)>=3:
        q=" ".join(sys.argv[2:])
        rs=search_gene(q)
        if not rs: print(f"No genes found for '{q}'")
        else:
            print(f"Found {len(rs)} gene(s) for '{q}':\n")
            for g in rs[:20]:
                print(f"  {g['symbol']} (GeneID: {g['gene_id']})")
                print(f"    {g['description']}")
                print(f"    Chromosome: {g['chromosome']} ({g['map_location']})")
                print(f"    Type: {g['type']}")
                if g["synonyms"]!="-": print(f"    Synonyms: {g['synonyms']}")
                print()
    elif cmd=="region" and len(sys.argv)>=3:
        parts=sys.argv[2].replace("-",":").split(":")
        if len(parts)!=3: print("Format: chr17:7668421-7687490"); return
        ch,s,e=parts[0],int(parts[1]),int(parts[2])
        print(f"Extracting {ch}:{s}-{e}...")
        sq=extract_region(ch,s,e)
        if sq:
            print(f"  Length: {len(sq)} bp")
            print(f"  First 80: {sq[:80]}")
            gc=(sq.upper().count('G')+sq.upper().count('C'))/max(1,len(sq))*100
            print(f"  GC: {gc:.1f}%")
            print(f"  SHA-256: {hashlib.sha256(sq.encode()).hexdigest()[:32]}")
        else: print("Region not found")
    elif cmd=="sequence" and len(sys.argv)>=3:
        q=sys.argv[2].upper()
        print(f"Searching {q} ({len(q)}bp) in GRCh38...")
        rs=search_sequence(q)
        if not rs: print("No matches")
        else:
            print(f"Found {len(rs)} match(es):")
            for r in rs: print(f"  {r['chromosome']}:{r['position']} ({r['strand']})")
    elif cmd=="stats": genome_stats()
    elif cmd=="json" and len(sys.argv)>=4:
        if sys.argv[2]=="gene": print(json.dumps(search_gene(" ".join(sys.argv[3:])),indent=2))
    else: print(__doc__)

if __name__=="__main__": main()
