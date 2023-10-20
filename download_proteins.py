# %%
from init import SAVE_DIR
import os
import requests
from traceback import print_stack


def parseUniProtRet(r):
    j = r.json()["results"]
    results = []
    for it in j:
        comm = None
        if "sequence" not in it:
            continue
        if "comments" not in it:
            pass
        elif not it["comments"]:
            pass
        else:
            comm = "\n".join(x["texts"][0]["value"] for x in it["comments"])
        results.append(
            dict(
                protein_name=it["uniProtkbId"],
                acc=it["primaryAccession"],
                desc=comm,
                seq=it["sequence"]["value"],
            )
        )
    return results


PROT_DIR = os.path.join(SAVE_DIR, "proteins")
if __name__ == "__main__":
    os.makedirs(PROT_DIR, exist_ok=True)
    # %%
    HEADER = {"Accept": "application/json"}
    URL = "https://rest.uniprot.org/uniprotkb/search?size=500&query=*&fields=id,protein_name,sequence,cc_function"
    import time
    from tqdm import tqdm

    results = []
    pcnt = 0
    parsedCnt = 0
    unparsedCnt = 0
    with requests.Session() as s:
        with tqdm() as pbar:
            while True:
                while True:
                    try:
                        r = s.get(URL, headers=HEADER)
                        break
                    except Exception as err:
                        print("Failed getting URL. Retrying, error was: \n")
                        print(err)

                if "Link" not in r.headers:
                    break
                URL = r.headers["Link"]
                URL = URL[1 : URL.find("; rel") - 1]
                d = parseUniProtRet(r)
                pparsedCnt = sum([1 for x in d if x["desc"]])
                punparsedCnt = len(d) - pparsedCnt
                results.extend(d)
                parsedCnt += pparsedCnt
                unparsedCnt += punparsedCnt
                pbar.update(pcnt)
                pbar.set_description(
                    str(parsedCnt + unparsedCnt)
                    + " parsed, "
                    + str(unparsedCnt)
                    + " have no content"
                )
                time.sleep(0.1)
                pcnt += 1

    # %%
    import json

    with open(os.path.join(PROT_DIR, "data.json"), "w") as out:
        json.dump({"data": results}, out)
