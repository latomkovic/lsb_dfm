import os
import set_bathy
import requests

target=set_bathy.merged_2m_path

url="https://tds.sfei.org/thredds/fileServer/lsb_master_dem/lsb_master_dem-20171116-test.tif"

if os.path.exists(target):
    raise Exception('%s already exists - will not download'%target)


r=requests.get(url,stream=True)
byte_sum=0
thresh=102400
bucket=0

with open(target,'wb') as fp:
    for chunk in r.iter_content(chunk_size=1024):
        if chunk:
            fp.write(chunk)
            byte_sum+=len(chunk)
            bucket+=len(chunk)
            if bucket>thresh:
                print("%6.3f Mbytes"%(byte_sum/1.e6))
                bucket=0



