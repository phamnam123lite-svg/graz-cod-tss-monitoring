"""Download and verify the original Zenodo archive, then extract locally."""
from pathlib import Path
import hashlib
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "graz/data-set.zip"
URL = "https://zenodo.org/records/15783832/files/data-set.zip?download=1"
EXPECTED_MD5 = "307525c96633e977e144357e014acfc8"

def main():
    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    if not ARCHIVE.exists():
        partial = ARCHIVE.with_suffix(".zip.part")
        request = urllib.request.Request(URL, headers={"User-Agent": "Graz-COD-TSS-reproducibility/1.0"})
        print("Downloading the original data archive from Zenodo...", flush=True)
        with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as dest:
            while chunk := response.read(1024 * 1024):
                dest.write(chunk)
        if hashlib.md5(partial.read_bytes()).hexdigest() != EXPECTED_MD5:
            raise ValueError("Downloaded archive checksum differs from the study source")
        partial.replace(ARCHIVE)
    if hashlib.md5(ARCHIVE.read_bytes()).hexdigest() != EXPECTED_MD5:
        raise ValueError("Archive checksum differs from the study source")
    target = (ROOT / "graz/raw").resolve()
    with zipfile.ZipFile(ARCHIVE) as source:
        for entry in source.infolist():
            if not (target / entry.filename).resolve().is_relative_to(target):
                raise ValueError("Unsafe archive path")
        source.extractall(target)
    print("Verified original archive and extracted raw data")

if __name__ == "__main__":
    main()
