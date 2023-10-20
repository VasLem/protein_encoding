# %%
from init import SAVE_DIR
import requests
import bs4
import re

parser = bs4.BeautifulSoup()
data = requests.get("https://awi.cuhk.edu.cn/dbPTM/download.php").text
import time
import os
from tqdm import tqdm

PTM_DIR = os.path.join(SAVE_DIR, "ptm")

PTM_SRC_DIR = os.path.join(SAVE_DIR, "ptm", "src")
PTM_EX_DIR = PTM_SRC_DIR.replace("src", "raw")
if __name__ == "__main__":
    # %% From experimental sites
    os.makedirs(PTM_SRC_DIR, exist_ok=True)
    # experimental
    for ptm in tqdm(re.findall("download\/experiment\/([\w-]+)\.gz", data)):
        outPath = os.path.join(PTM_SRC_DIR, f"{ptm}.gz")
        if os.path.exists(outPath):
            continue
        with open(outPath, "wb") as out:
            out.write(
                requests.get(
                    f"https://awi.cuhk.edu.cn/dbPTM/download/experiment/{ptm}.gz",
                    allow_redirects=True,
                ).content
            )
        time.sleep(0.1)
    # %%
    import shutil
    import tarfile

    shutil.rmtree(PTM_EX_DIR, ignore_errors=True)
    os.makedirs(PTM_EX_DIR, exist_ok=True)
    for ptm in os.listdir(PTM_SRC_DIR):
        if not ptm.endswith(".gz"):
            continue
        ptm = ptm.replace(".gz", "")
        inPath = os.path.join(PTM_SRC_DIR, ptm + ".gz")
        outPath = os.path.join(PTM_EX_DIR, ptm)
        with tarfile.open(inPath, "r") as inp:
            inp.extract(ptm, outPath)

    # %%
    import pandas as pd

    ts = []
    failed = []
    for ptm in os.listdir(PTM_EX_DIR):
        inPath = os.path.join(PTM_EX_DIR, ptm, ptm)
        ts.append(pd.read_csv(inPath, sep="\t", header=None))

    df = pd.concat(ts, axis=0)
    # %%
    df.head()
    # %%
    df = df[[0, 1, 2, 3]]
    df.columns = ["protein", "accession", "residue", "ptm"]

    # %%
    df.head()
    # %%
    df.duplicated().any()
    # %%
    outPath = os.path.join(SAVE_DIR, "ptm", "data.csv")
    df.to_csv(outPath, index=False)
    # %%
    with open(outPath, "r") as inp:
        for c, line in enumerate(inp.readlines()):
            print(line)
            if c == 9:
                break
    # %%
