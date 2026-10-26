"""Retrieve whole public matrices and platform annotations; no group effects."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import datetime,gzip,hashlib,json,time,urllib.request

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'sources/biomarker_transfer_v13/clinical_inputs';OUT.mkdir(exist_ok=True)
resources=[]
for accession in ['GSE51981','GSE120103']:
    name=accession+'_series_matrix.txt.gz'
    resources.append((accession,name,f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/matrix/{name}'))
for accession in ['GPL570','GPL6480']:
    name=accession+'.annot.gz'
    resources.append((accession,name,f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{accession[:-3]}nnn/{accession}/annot/{name}'))

def retrieve(item):
    accession,name,url=item;record={'accession':accession,'filename':name,'url':url}
    try:
        req=urllib.request.Request(url,headers={'Range':'bytes=0-0','Accept-Encoding':'identity'})
        with urllib.request.urlopen(req,timeout=35) as response:
            assert response.status==206
            content_range=response.headers['Content-Range'];assert content_range.startswith('bytes 0-0/')
            expected=int(content_range.split('/')[1]);assert response.read()==b'\x1f'
            record.update(expected_server_bytes=expected,last_modified=response.headers.get('Last-Modified'),
                          initial_response_headers=dict(response.headers))
        path=OUT/name;partial=OUT/(name+'.partial')
        if not path.exists():
            offset=partial.stat().st_size if partial.exists() else 0;assert offset<=expected
            with partial.open('ab') as dest:
                while offset<expected:
                    end=min(offset+2*1024*1024,expected)-1
                    for attempt in range(3):
                        try:
                            request=urllib.request.Request(url,headers={'Range':f'bytes={offset}-{end}','Accept-Encoding':'identity'})
                            with urllib.request.urlopen(request,timeout=35) as response:
                                assert response.status==206 and response.headers['Content-Range']==f'bytes {offset}-{end}/{expected}'
                                block=response.read()
                            assert len(block)==end-offset+1;break
                        except Exception:
                            if attempt==2:raise
                            time.sleep(1)
                    dest.write(block);dest.flush();offset=end+1
            assert partial.stat().st_size==expected
            gzip.decompress(partial.read_bytes());partial.rename(path)
        compressed=path.read_bytes();assert len(compressed)==expected
        payload=gzip.decompress(compressed)
        record.update(status='complete_bytes_and_gzip_CRC',bytes=len(compressed),
                      sha256=hashlib.sha256(compressed).hexdigest(),uncompressed_bytes=len(payload),
                      uncompressed_sha256=hashlib.sha256(payload).hexdigest())
    except Exception as exc:
        record.update(status='retrieval_failed_not_no_data',error_type=type(exc).__name__,error=str(exc))
    print(accession,record['status'],flush=True)
    return record

with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(retrieve,resources))
(OUT/'retrieval.json').write_text(json.dumps({'completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'records':records,'scope':'Technical input retrieval only; no group expression contrasts or prediction models.'},indent=2)+'\n')
assert all(x['status']=='complete_bytes_and_gzip_CRC' for x in records)
