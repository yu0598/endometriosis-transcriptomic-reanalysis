"""Retrieve one official raw archive and validate/extract its 36 FE gzip files."""
import csv
import gzip
import hashlib
from pathlib import Path
import tarfile
import urllib.request

root = Path.cwd()
out = root / 'data_raw/GSE120103_rawFE_v13'
out.mkdir(parents=True, exist_ok=True)
source = root / 'outputs/biomarker_transfer_v13/phase2_investigation/GSE120103_processing_sources'
filelist = list(csv.DictReader(open(source / 'GSE120103_filelist.txt'), delimiter='\t'))
archive_record = next(r for r in filelist if r['#Archive/File'] == 'Archive')
expected = {r['Name']: int(r['Size']) for r in filelist if r['#Archive/File'] == 'File'}
assert len(expected) == 36
url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120103/suppl/' + archive_record['Name']
expected_size = int(archive_record['Size'])
with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=30) as response:
    server_size = int(response.headers['Content-Length'])
assert server_size == expected_size, (server_size, expected_size)
archive = out / archive_record['Name']
partial = archive.with_suffix('.tar.part')
if not archive.exists():
    offset = partial.stat().st_size if partial.exists() else 0
    assert offset <= expected_size
    request = urllib.request.Request(url, headers={'Range': f'bytes={offset}-'} if offset else {})
    with urllib.request.urlopen(request, timeout=45) as response:
        if offset:
            assert response.status == 206
            assert response.headers['Content-Range'].startswith(f'bytes {offset}-')
        with open(partial, 'ab' if offset else 'wb') as handle:
            checkpoint = offset // (10 * 1024 * 1024)
            while block := response.read(1024 * 1024):
                handle.write(block)
                current = handle.tell()
                if current // (10 * 1024 * 1024) > checkpoint:
                    checkpoint = current // (10 * 1024 * 1024)
                    print(f'Downloaded {current}/{expected_size} bytes', flush=True)
    assert partial.stat().st_size == expected_size
    partial.rename(archive)
assert archive.stat().st_size == expected_size
archive_sha = hashlib.file_digest(open(archive, 'rb'), 'sha256').hexdigest()
rows = []
with tarfile.open(archive, 'r:') as tar:
    members = tar.getmembers()
    assert len(members) == 36
    assert {m.name for m in members} == set(expected)
    for member in members:
        assert member.isfile() and Path(member.name).name == member.name
        assert member.size == expected[member.name]
        path = out / member.name
        if not path.exists():
            with tar.extractfile(member) as source_handle, open(path, 'wb') as dest:
                while block := source_handle.read(1024 * 1024):
                    dest.write(block)
        assert path.stat().st_size == member.size
        compressed_sha = hashlib.file_digest(open(path, 'rb'), 'sha256').hexdigest()
        digest = hashlib.sha256()
        nbytes = 0
        with gzip.open(path, 'rb') as handle:
            while block := handle.read(1024 * 1024):
                digest.update(block)
                nbytes += len(block)
        rows.append(dict(gsm=member.name.split('_')[0], file=member.name,
                         filelist_bytes=expected[member.name], extracted_bytes=path.stat().st_size,
                         compressed_sha256=compressed_sha, uncompressed_bytes=nbytes,
                         uncompressed_sha256=digest.hexdigest(), gzip_crc_eof_verified=True,
                         archive_url=url, archive_bytes=expected_size,
                         server_content_length=server_size, archive_sha256=archive_sha,
                         safe_regular_basename_member=True, tar_header_checksum_verified=True))
        print('Validated', member.name.split('_')[0], flush=True)
with open(out / 'rawFE_retrieval_manifest.csv', 'w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print('COMPLETE', len(rows), 'archive_sha256', archive_sha, flush=True)
