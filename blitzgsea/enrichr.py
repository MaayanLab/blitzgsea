from typing import List
import itertools
import urllib.request
import shutil
import json
import os
import re
import ssl

def get_library(library: str):
    """
    Load gene set library from Enrichr as dictionary.

    Parameters:
    signature (string): Name of Enrichr gene set library.

    Returns:
    Gene set library as dictionary.
    """
    return read_gmt(load_library(library))

def list_libraries():
    print(get_config())
    return(load_json(get_config()["LIBRARY_LIST_URL"])["library"])

def load_library(library: str, overwrite: bool = False, verbose: bool = False) -> str:
    if not os.path.exists(get_data_path()+library or overwrite):
        if verbose:
            print("Download Enrichr geneset library")
        urlretrieve(get_config()["LIBRARY_DOWNLOAD_URL"]+library, get_data_path()+library)
    else:
        if verbose:
            print("File cached. To reload use load_library(\""+library+"\", overwrite=True) instead.")
    lib = read_gmt(get_data_path()+library)
    if verbose:
        print("# genesets: "+str(len(lib)))
    return(get_data_path()+library)

def print_libraries():
    """
    Retrieve all names of gene set libraries from Enrichr.

    Returns:
    Array of gene set library names.
    """
    libs = list_libraries()
    for i in range(0, len(libs)):
        print(str(i)+" - "+libs[i])

def read_gmt(gmt_file: str, background_genes: List[str]=[], verbose=False):
    file = open(gmt_file, 'r')
    lines = file.readlines()
    library = {}
    background_set = {}
    if len(background_genes) > 1:
        background_genes = [x.upper() for x in background_genes]
        background_set = set(background_genes)
    for line in lines:
        sp = line.strip().split("\t")
        sp2 = [re.sub(",.*", "",value) for value in sp[2:]]
        sp2 = [x.upper() for x in sp2 if x]
        if len(background_genes) > 2:
            geneset = list(set(sp2).intersection(background_set))
            if len(geneset) > 0:
                library[sp[0]] = geneset
        else:
            if len(sp2) > 0:
                library[sp[0]] = sp2
    ugenes = list(set(list(itertools.chain.from_iterable(library.values()))))
    if verbose:
        print("Library loaded. Library contains "+str(len(library))+" gene sets. "+str(len(ugenes))+" unique genes found.")
    return library

def load_json(url):
    context = ssl._create_unverified_context()
    req = urllib.request.Request(url)
    r = urllib.request.urlopen(req, context=context).read()
    return(json.loads(r.decode('utf-8')))

def get_config():
    config_url = os.path.join(
        os.path.dirname(__file__),
        'data/config.json')
    with open(config_url) as json_file:
        data = json.load(json_file)
    return(data)

def get_data_path() -> str:
    path = os.path.join(
        os.path.dirname(__file__),
        'data/'
    )
    return(path)

def urlretrieve(req, filename):
    context = ssl._create_unverified_context()
    with urllib.request.urlopen(req, context=context) as fr:
        with open(filename, 'wb') as fw:
            shutil.copyfileobj(fr, fw)