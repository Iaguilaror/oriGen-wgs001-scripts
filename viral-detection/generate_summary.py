# Author: Eugenio Guzmán <eugenio.guzman@tec.mx>

import datetime
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs, flush=True)

eprint(datetime.datetime.now(), "importing libraries")

import os
import numpy as np
import pandas as pd
import threading
import queue
import time
from itertools import combinations

# INPUT: Sample list (oriGen-1427.txt)
samps = set()
with open("oriGen-1427.txt") as f:
    line_ = f.readline()
    while line_:
        line = line_[:-1]
        samps.add(line)
        line_ = f.readline()

# INPUT: Virus database index (virusdb.fa.fai)
vdbi = pd.read_csv("virusdb.fa.fai", sep="\t", names=["chrom", "len", "offset", "linesize", "fulllinesize"]).set_index("chrom")

dt_fdf = np.dtype([
    ('sample', np.uint16),
    ('qseqid', np.uint64),
    ('flag', np.uint16),
    ('mapq', np.uint8),
    ('hgchrom', np.uint16),
    ('hgpos', np.uint32),
    ('hgmatechrom', np.uint16),
    ('hgmatepos', np.uint32),
    ('qlen', np.uint8),
    ('sseqid', np.uint16),
    ('pident', np.uint8),
    ('length', np.uint8),
    ('mismatch', np.uint8),
    ('gapopen', np.uint8),
    ('qstart', np.uint8),
    ('qend', np.uint8),
    ('sstart', np.uint32),
    ('send', np.uint32),
    ('evalue', np.int16),
    ('bitscore', np.uint16),
])


# INPUT: Reference genome database index (Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai)
chroms = pd.read_csv("Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai", sep="\t", names=["chrom", "len", "offset", "linesize", "fulllinesize"])

# INPUT: Virus database index (virusdb.fa.fai)
genomes = pd.read_csv("virusdb.fa.fai", sep="\t", names=["id", "len", "start", "linesize", "linesize2"])

version = "000"
filename_base = f"viralsummary_{version}.{{}}.parquet"
idx_genomes = genomes["id"].to_dict()

# INPUT: Sample list (oriGen-1427.txt)
idx_samples = pd.read_csv("oriGen-1427.txt", sep="\t", names=["id"])["id"].to_dict()

idx_chroms = chroms["chrom"].to_dict()
idx_chroms[np.iinfo(np.uint16).max] = "*"

chroms = chroms.set_index("chrom")
chromid_to_chrom = {i: k for i, k in enumerate(chroms.index)}

n_worker_threads = 10
n_sectors = 10
rcounts_filter_thresh = 10
scounts_filter_thresh = 3
chunk_size = int(1e6)
notif_interval = int(5e7)
limit = None

os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

work_queue = queue.Queue()
lock = threading.Lock()
thread_local = threading.local()
M = np.zeros((len(idx_genomes),len(idx_genomes)), np.uint64)
B = np.zeros((len(idx_genomes)), np.float64)

def worker(wq):
    global lock, M, B
    thread_local.M = np.zeros((len(idx_genomes),len(idx_genomes)), np.uint64)
    thread_local.B = np.zeros((len(idx_genomes)), np.float64)
    while True:
        task = wq.get()
        if task is None:
            break
        combos = task
        if len(combos) >= 2:
            for (i,j) in combinations(combos, 2):
                thread_local.M[i,j] += 1
        
        wq.task_done()
    with lock:
        M += thread_local.M


worker_threads = []
for _ in range(n_worker_threads):
    thread = threading.Thread(target=worker, args=(work_queue,))
    thread.start()
    worker_threads.append(thread)

rcounts = pd.Series(name="read_count", index=pd.MultiIndex.from_arrays([np.array([], dtype=np.uint16), np.array([], dtype=np.uint16)], names=["sample", "sseqid"]))
scounts = pd.DataFrame(index=pd.MultiIndex.from_arrays([np.array([], dtype=np.uint16), np.array([], dtype=np.uint16)], names=["sample", "sseqid"]), columns=range(n_sectors))

readn = None
read_chunks = 0
acc_lines = 0
next_notif = notif_interval
M_diag = np.zeros(len(idx_genomes), np.uint64)

eprint(datetime.datetime.now(), "about to read data")

# INPUT: Binary outfmt6 data (out.data)
with open('out.data', 'rb') as f:
    carry = pd.DataFrame()
    while True:
        chunk = np.fromfile(f, dtype=dt_fdf, count=chunk_size)[["sample", "qseqid", "sseqid", "flag", "bitscore", "mapq", "sstart", "send"]]
        acc_lines += chunk.shape[0]
        read_chunks += 1
        if acc_lines >= next_notif:
            next_notif += notif_interval
            eprint(datetime.datetime.now(), f"lines so far {acc_lines}, chunks so far: {read_chunks}")
        no_more_chunks = True if not chunk.size else False
        if not no_more_chunks:
            chunk_df = pd.DataFrame(chunk)
            curr_mateid = (
                chunk_df["qseqid"].iloc[-1] + ((chunk_df["flag"].iloc[-1] & 128) // 128),
                int(chunk_df["sample"].iloc[-1])
            )
            chunk_df = chunk_df[(chunk_df["flag"] & 0x04) != 0].copy()
            chunk_df["mateid"] = chunk_df["qseqid"] + ((chunk_df["flag"] & 128) // 128)
            islast = (chunk_df["mateid"] == curr_mateid[0]) & (chunk_df["sample"] == curr_mateid[1])
            new_carry = chunk_df[islast]
            chunk_df = pd.concat([carry, chunk_df[~islast]], ignore_index=True)
            carry = new_carry
        else:
            chunk_df = carry

        b = ((chunk_df["sample"].values[:-1] == chunk_df["sample"].values[1:]) & (chunk_df["mateid"].values[:-1] == chunk_df["mateid"].values[1:]))
        r = np.concatenate([[0],(np.arange(chunk_df.shape[0]-1) + 1)[~b]])
        d = np.concatenate([np.diff(r), [chunk_df.shape[0] - r[-1]]])
        s = chunk_df["sseqid"].values
        bs = chunk_df["bitscore"].values
        rde1 = r[d==1]
        dgt1 = d>1
        
        for i, bi in zip(s[rde1], bs[rde1]):
            M_diag[i] += 1
            B[i] += bi/10

        for ri, di in zip(r[dgt1], d[dgt1]):
            while work_queue.qsize() >= 1000:
                while work_queue.qsize() >= 100:
                    time.sleep(0.001)
            ss = s[ri:ri+di]
            if np.all(ss[:-1] == ss[1:]):
                M_diag[ss[0]] += 1
                B[ss[0]] += np.mean(b[ri:ri+di]/10)
            else:
                u = np.unique(ss)
                work_queue.put(u)

        chunk_df["sector"] = ((chunk_df["sstart"]+chunk_df["send"])/(chunk_df["sseqid"].map(genomes["len"]))*(10//2)).astype(int)
        svcs = chunk_df[["sample", "sseqid", "sector"]].value_counts()
        #A[svcs.index.get_level_values(0), svcs.index.get_level_values(1), svcs.index.get_level_values(2)] += svcs
        raw_sector_pairs = svcs.reset_index().pivot(index=["sample", "sseqid"], columns="sector", values="count")
        scounts = scounts.add(raw_sector_pairs, fill_value=0).astype(float)

        raw_pairs = chunk_df.groupby(["sample", "sseqid"], sort=False)["mateid"].nunique().rename("read_count")
        rcounts = rcounts.add(raw_pairs, fill_value=0)

        if no_more_chunks:
            break

        if limit and read_chunks >= limit:
            break
work_queue.join()

for _ in range(n_worker_threads):
    work_queue.put(None)

for thread in worker_threads:
    thread.join()

os.environ.pop('OMP_NUM_THREADS', None)
os.environ.pop('OPENBLAS_NUM_THREADS', None)
os.environ.pop('MKL_NUM_THREADS', None)
os.environ.pop('VECLIB_MAXIMUM_THREADS', None)
os.environ.pop('NUMEXPR_NUM_THREADS', None)
eprint(datetime.datetime.now(), "finishing step 1")

np.fill_diagonal(M, np.diagonal(M) + M_diag)
#B += B_local
B /= np.diagonal(M)
unique_bitscores = pd.Series(B, index=list(idx_genomes.values()), name="unique_bitscore").rename_axis(index="sseqid").to_frame().reset_index()
unique_bitscores.to_parquet(filename_base.format("unique_bitscores"))

cm = pd.DataFrame(M+M.T-np.triu(np.tril(M)),index=list(idx_genomes.values()),columns=list(idx_genomes.values()))
cm.to_parquet(filename_base.format("comappings"))
rcounts_raw = rcounts.astype(int)
rcounts = rcounts_raw.to_frame().reset_index()
rcounts["sseqid"] = rcounts["sseqid"].astype(pd.CategoricalDtype(idx_genomes.keys())).cat.rename_categories(idx_genomes)
rcounts["sample"] = rcounts["sample"].astype(pd.CategoricalDtype(idx_samples.keys())).cat.rename_categories(idx_samples)
rcounts = rcounts.convert_dtypes()
rcounts.to_parquet(filename_base.format("raw_counts"))

scounts_raw = scounts.fillna(0).astype(int)
scounts_raw = scounts_raw.loc[np.sum(scounts_raw, axis=1) > 0]
scounts = scounts_raw.fillna(0).astype(int).reset_index()
scounts["sseqid"] = scounts["sseqid"].astype(pd.CategoricalDtype(idx_genomes.keys())).cat.rename_categories(idx_genomes)
scounts["sample"] = scounts["sample"].astype(pd.CategoricalDtype(idx_samples.keys())).cat.rename_categories(idx_samples)
scounts.columns = scounts.columns.astype(str)
scounts.to_parquet(filename_base.format("f4_sector_counts"))

f1 = (np.sum(scounts_raw>0,axis=1) >= scounts_filter_thresh)
f2 = (rcounts_raw >= rcounts_filter_thresh)
f12 = f1 & f2

combos_pass_filter = rcounts_raw.loc[f12]
#adf = pd.DataFrame(A.reshape(len(idx_samples),len(idx_genomes)*n_sectors), columns=pd.MultiIndex.from_product([genomes["id"], map(str, range(n_sectors))], names=["sseqid", "sector"]), index=idx_samples.values())
#adf = adf.loc[:, (adf>0).any()]
#adf.to_parquet(filename_base.format("f4_sector_counts"))
eprint(datetime.datetime.now(), f"{combos_pass_filter.shape[0]} combinations pass filter (f1: {f1.sum()}, f2: {f2.sum()})")
eprint(datetime.datetime.now(), "ready for step 2")

chunk_size = int(1e6)
read_chunks = 0
acc_lines = 0
next_notif = notif_interval
with open('out.data', 'rb') as f:
    outs = []
    carry = pd.DataFrame()
    while True:
        chunk = np.fromfile(f, dtype=dt_fdf, count=chunk_size)[[
            "sample", "qseqid", "sseqid", "flag", "hgmatechrom",
            "bitscore", "pident", "hgmatepos", "mapq", "length",
            "qlen", "qstart", "qend"
        ]]
        acc_lines += chunk.shape[0]
        read_chunks += 1
        if acc_lines >= next_notif:
            next_notif += notif_interval
            eprint(datetime.datetime.now(), f"lines so far {acc_lines}, chunks so far: {read_chunks}")
        no_more_chunks = True if not chunk.size else False
        if not no_more_chunks:
            chunk_df = pd.DataFrame(chunk)
            curr_sample = int(chunk_df.iloc[-1]["sample"])
            chunk_df = chunk_df[(chunk_df["flag"] & (0x04|0x08)) != 0]
            chunk_df = chunk_df.merge(combos_pass_filter, on=["sample", "sseqid"]).drop(columns="read_count").copy()
            chunk_df["mateid"] = chunk_df["qseqid"] + ((chunk_df["flag"] & 128) // 128)
            carry = pd.concat([carry, chunk_df], ignore_index=True)
            islast = (chunk_df["sample"] == curr_sample)
            carry = chunk_df[islast]
            chunk_df = chunk_df[~islast]
        else:
            chunk_df = carry
            carry = pd.DataFrame()

        chunk_df = chunk_df.merge(combos_pass_filter, on=["sample", "sseqid"], how="left")
        chunk_df = chunk_df.groupby(["sample", "sseqid", "qseqid"]).filter(lambda x: ((x["flag"] & 4) != 0).any())
        if chunk_df.shape[0] > 0:
            chunk_df["bitscore"] /= 10
            chunk_df["pident"] /= np.iinfo(np.uint8).max
            chunk_df["qcov"] = ((chunk_df["qend"]-chunk_df["qstart"]).abs()+1)/chunk_df["qlen"]
            chunk_df["qcovident"] = chunk_df["pident"]*(chunk_df["length"]/chunk_df["qlen"])
            chunk_df = chunk_df.merge((chunk_df.groupby(["sample", "sseqid", "qseqid"], sort=False)["mateid"].nunique() > 1).rename("mate_maps_to_same_virus"), on=["sample", "sseqid", "qseqid"], how="left")
            chromposses = chunk_df[~chunk_df["mate_maps_to_same_virus"] & ((chunk_df["flag"] & 8) == 0)].groupby(["sample", "sseqid", "qseqid"], sort=False).first()[["hgmatechrom", "hgmatepos"]].reset_index()
            best_chroms = chromposses[["sample", "sseqid", "hgmatechrom"]].value_counts(sort=False).reset_index().set_index("hgmatechrom").groupby(["sample", "sseqid"], sort=False)["count"].agg(chrom="idxmax", mates_at_chrom="max")
            chromposses = chromposses.merge(best_chroms, on=["sample", "sseqid"], how="left")
            pos_dists = chromposses[chromposses["hgmatechrom"] == chromposses["chrom"]].groupby(["sample", "sseqid"], sort=False)["hgmatepos"].agg(pos_mean="mean",pos_stdev="std")

            outs.append(best_chroms.join(pos_dists, how="outer").join(chunk_df.groupby(["sample", "sseqid"], sort=False).agg(**{
                "high_hg_mapq_reads": ("mapq", lambda x: (x > 10).sum()),
                "mean_hg_mapq": ("mapq", "mean"),
                "fraction_hg_unmapped": ("flag", lambda x: ((x & 4) > 0).mean()),
                "mean_bitscore": ("bitscore", "mean"),
                "mean_identity": ("pident", "mean"),
                "mean_coverage": ("qcov", "mean"),
                "mean_identical_coverage": ("qcovident", "mean"),
                "mean_alnlen": ("length", "mean"),
                "complete_pairs": ("mate_maps_to_same_virus", "sum"),
                "read_count": ("mateid", "nunique"),
                "mapped_mate_count": ("hgmatechrom", lambda x: (x<np.iinfo(np.uint16).max).sum()),
            }), how="outer"))

        if no_more_chunks:
            break

        if limit and read_chunks >= limit:
            break

eprint(datetime.datetime.now(), "finishing step 2")
df = pd.concat(outs).reset_index()

df["sseqid"] = df["sseqid"].astype(pd.CategoricalDtype(idx_genomes.keys())).cat.rename_categories(idx_genomes)
df["sample"] = df["sample"].astype(pd.CategoricalDtype(idx_samples.keys())).cat.rename_categories(idx_samples)
df["chrom"] = df["chrom"].astype(pd.CategoricalDtype(chromid_to_chrom.keys())).cat.rename_categories(chromid_to_chrom.values())
df = df.convert_dtypes().merge(rcounts.rename(columns={"read_count": "flag4_readcount"}), how="left", on=["sample", "sseqid"])
df.to_parquet(filename_base.format("sample_vs_virus"))