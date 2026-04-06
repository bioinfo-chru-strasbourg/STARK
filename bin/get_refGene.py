import os
import re
import gzip
import shutil
import argparse
import requests
from typing import List, Optional, Tuple

# Usage:
# python get_refGene.py -a hg19 -o refGene.gtf -p p13 #--keep-gz

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/"

ASSEMBLY_ALIASES = {
    "hg19": "GRCh37",
    "hg38": "GRCh38",
}


def list_assembly_directories() -> List[str]:
    """
    Fetch and parse the FTP directory listing to retrieve all assembly folders.
    """
    response = requests.get(BASE_URL)
    response.raise_for_status()

    pattern = r'href="(GCF_[^"/]+)/"'
    return re.findall(pattern, response.text)


def filter_assemblies(directories: List[str], assembly: str) -> List[str]:
    """
    Filter directories matching the requested assembly.
    """
    if assembly not in ASSEMBLY_ALIASES:
        raise ValueError(f"Unknown assembly: {assembly}")

    grc_name = ASSEMBLY_ALIASES.get(assembly, assembly)
    return [d for d in directories if grc_name in d]


def extract_patch_and_version(directory: str) -> Tuple[int, int]:
    """
    Extract (patch_number, accession_version) from directory name.

    Example:
    GCF_000001405.25_GRCh37.p13 --> (13, 25)
    """
    patch_match = re.search(r"\.p(\d+)", directory)
    version_match = re.search(r"GCF_\d+\.(\d+)", directory)

    patch = int(patch_match.group(1)) if patch_match else -1
    version = int(version_match.group(1)) if version_match else -1

    return patch, version


def select_best_match(directories: List[str], patch: Optional[str]) -> str:
    """
    Select the best matching directory:
    - If patch specified --> exact match
    - Otherwise --> highest patch, then highest accession version
    """
    if not directories:
        raise ValueError("No matching assemblies found")

    if patch:
        for d in directories:
            if d.endswith(f".{patch}"):
                return d
        raise ValueError(f"Patch {patch} not found")

    # Sort by (patch, version)
    sorted_dirs = sorted(
        directories, key=lambda d: extract_patch_and_version(d), reverse=True
    )

    return sorted_dirs[0]


def build_gtf_url(directory: str) -> str:
    """
    Construct GTF file URL.
    """
    return f"{BASE_URL}{directory}/{directory}_genomic.gtf.gz"


def download_file(url: str, output_path: str) -> str:
    """
    Download file from URL.
    """
    response = requests.get(url, stream=True)
    response.raise_for_status()

    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)

    return output_path


def decompress_gzip(input_path: str, output_path: str) -> str:
    """
    Decompress a .gz file.
    """
    with gzip.open(input_path, "rb") as f_in:
        with open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return output_path


def resolve_output_name(downloaded_file: str, user_output: Optional[str]) -> str:
    """
    Determine final output filename.
    """
    if user_output:
        return user_output

    # Default: remove .gz
    if downloaded_file.endswith(".gz"):
        return downloaded_file[:-3]

    return downloaded_file


def main():
    parser = argparse.ArgumentParser(
        description="Download GTF gene annotation from NCBI RefSeq"
    )

    parser.add_argument(
        "-a",
        "--assembly",
        default="hg19",
        choices=["hg19", "hg38"],
        help="Genome assembly (default: hg19)",
    )

    parser.add_argument(
        "-p",
        "--patch",
        help="Patch version (e.g., p13). If omitted, latest is selected.",
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output filename (default: same as downloaded, without .gz)",
    )

    parser.add_argument(
        "-k", "--keep-gz", action="store_true", help="Keep the downloaded .gz file"
    )

    args = parser.parse_args()

    print(f"[INFO] Assembly: {args.assembly}")
    print(f"[INFO] Patch: {args.patch or 'latest'}")

    # Step 1: discover directories
    directories = list_assembly_directories()

    # Step 2: filter
    filtered = filter_assemblies(directories, args.assembly)

    # Step 3: select best match
    selected = select_best_match(filtered, args.patch)
    print(f"[INFO] Selected: {selected}")

    # Step 4: build URL
    url = build_gtf_url(selected)
    print(f"[INFO] Downloading: {url}")

    gz_filename = url.split("/")[-1]

    # Step 5: download
    download_file(url, gz_filename)

    # Step 6: resolve output name
    output_file = resolve_output_name(gz_filename, args.output)

    print(f"[INFO] Decompressing --> {output_file}")

    # Step 7: decompress
    decompress_gzip(gz_filename, output_file)

    # Step 8: cleanup
    if not args.keep_gz:
        os.remove(gz_filename)

    print(f"[SUCCESS] File ready: {output_file}")


if __name__ == "__main__":
    main()
