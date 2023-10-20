# %%
from init import SAVE_DIR
from download_proteins import PROT_DIR
from download_ptms import PTM_DIR
import json
import pandas as pd
import os

from tqdm import tqdm

tqdm.pandas()

TRANS_DIR = os.path.join(SAVE_DIR, "input_to_prot_encoder")
PTM_ENC_PATH = os.path.join(TRANS_DIR, "ptms_encoder.pt")
MERGED_PATH = os.path.join(TRANS_DIR, "merged.json")
# %%
if __name__ == "__main__":
    os.makedirs(TRANS_DIR, exist_ok=True)
    # %%
    print("Loading data..")
    with open(os.path.join(PROT_DIR, "data.json"), "r") as inp:
        proteins = json.load(inp)["data"]

    ptms = pd.read_csv(os.path.join(PTM_DIR, "data.csv"), index_col=None)
    # %%
    ptms.shape, ptms.head()
    # %%
    ptms.duplicated().any(axis=0)
    # %%
    from sklearn.base import BaseEstimator, TransformerMixin
    import numpy as np

    class PTMEncoder(TransformerMixin, BaseEstimator):
        def __init__(self):
            self.ptms = {}

        @property
        def inv_ptms(self):
            return {x: v for v, x in self.ptms.items()}

        def fit(self, ptms):
            assert all(isinstance(x, str) for x in ptms)
            self.ptms = {
                p: chr(ord("Z") + c + 1) for c, p in enumerate(sorted(set(ptms)))
            }
            return self

        def partial_fit(self, ptms):
            assert all(isinstance(x, str) for x in ptms)
            self.ptms.update(
                {
                    p: chr(ord("Z") + c + 1 + len(self.ptms))
                    for c, p in enumerate(sorted(set(ptms)))
                }
            )
            return self

        def transform(self, ptms):
            is_str = isinstance(ptms, str)
            if is_str:
                ptms = [ptms]
            ret = [self.ptms[x] for x in ptms]
            if is_str:
                return ret[0]
            if isinstance(ptms, (np.ndarray, pd.Series)):
                r = ptms[:]
                r[:] = ret
                ret = r
            return ret

        def inverse_transform(self, ptms):
            is_str = isinstance(ptms, str)
            if is_str:
                ptms = [ptms]
            ret = [self.inv_ptms[x] for x in ptms]
            if is_str:
                return ret[0]
            if isinstance(ptms, (np.ndarray, pd.Series)):
                r = ptms[:]
                r[:] = ret
                ret = r
            return ret

    print("Encoding PTMs..")
    import joblib

    ptmsEncoder = PTMEncoder().fit(ptms["ptm"])

    joblib.dump(ptmsEncoder, PTM_ENC_PATH)

    ptms["ptm_t"] = ptmsEncoder.transform(ptms["ptm"])

    # %% merge pmt into proteins
    print("Checking ptms and proteins matches..")

    # %% there are duplicates
    proteins = pd.DataFrame(proteins).head()
    proteins.drop_duplicates(inplace=True)
    # %%
    proteins.duplicated("protein_name").any()
    assert len(proteins["protein_name"]) == len(proteins)
    # %%
    pnames = proteins["protein_name"].values.tolist()
    accnames = proteins["acc"].values.tolist()
    sequences = {}

    ptmspnames = set(ptms["protein"].unique())
    found = list(ptmspnames.intersection(pnames))
    notFoundPtmProteins = list(ptmspnames.difference(found))
    if notFoundPtmProteins:
        print("Proteins in ptm db not found in uniprot:", len(notFoundPtmProteins))
        print("Out of", len(ptmspnames))
    # %% get missing proteins
    notFoundPtmProteins = [x for x in notFoundPtmProteins if isinstance(x, str)]
    meanSize = np.mean([len(x) for x in notFoundPtmProteins if x]) + 1  # the + operator
    numQueries = int(min(2000 / meanSize, 500))
    from download_proteins import parseUniProtRet

    HEADER = {"Accept": "application/json"}

    def urlGen(x):
        return (
            "https://rest.uniprot.org/uniprotkb/search?size=500&query="
            f"{' || '.join(x)}*"
            "&fields=id,protein_name,sequence,cc_function"
        )

    import requests
    import time

    add = []
    print("Downloading missing proteins..")
    iterator = tqdm(list(range(len(notFoundPtmProteins) // numQueries + 1)))
    for i in iterator:
        subset = notFoundPtmProteins[numQueries * i : numQueries * (i + 1)]
        url = urlGen(subset)
        while True:
            try:
                ret = requests.get(
                    url,
                    headers=HEADER,
                )
                break
            except BaseException:
                time.sleep(0.5)
        x = parseUniProtRet(ret)
        print(x)
        break
        add.extend(x)
        iterator.set_description(
            f"Found {len(x)} out of {len(subset)}. Total: {len(add)}"
        )
        time.sleep(0.2)
    if add:
        proteins.extend(add)
        found.extend([x["protein_name"] for x in add])
        with open(os.path.join(PROT_DIR, "data.json"), "w") as out:
            json.dump({"data": proteins}, out)
    notFoundPtmProteins = list(ptmspnames.difference(found))
    if add:
        print(
            "Updated number of proteins not found in UniProt after requerying:",
            len(notFoundPtmProteins),
        )

    # %%
    print("Merging proteins with ptms..")

    def makeSequence(subset):
        sequence = [x for x in proteins[pnames.index(protein_name)]["seq"]]
        for i, s in subset.iterrows():
            try:
                sequence[s["residue"] - 1] += s["ptm_t"]
            except IndexError:
                print(
                    f'Invalid residue {s["residue"]} for protein {protein_name} with length {len(sequence)}'
                )
                pass
        return sequence

    ptms_seqs = (
        ptms[ptms["protein"].isin(found)]
        .groupby("protein")
        .progress_apply(makeSequence)
    )
    # %%
    notFound = [x for x in pnames if x not in found]
    sequences = ptms_seqs.copy()
    for c, protein_name in tqdm(enumerate(notFound)):
        if protein_name not in found:
            sequences[protein_name] = " ".join([x for x in proteins[c]["seq"]])
            continue

    # %%
    import json

    with open(MERGED_PATH, "w") as out:
        json.dump(sequences, out)
