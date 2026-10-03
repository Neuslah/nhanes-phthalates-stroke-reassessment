"""Acquire verified official sources outside the release candidate (Python stdlib)."""
from pathlib import Path
import argparse,csv,hashlib,json,os,subprocess,urllib.request,shutil

def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig')))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest().upper()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--package-root',type=Path,default=Path(__file__).resolve().parents[1])
 p.add_argument('--input-root',type=Path,required=True)
 p.add_argument('--source-cache',type=Path,required=True)
 p.add_argument('--translated-cache',type=Path)
 p.add_argument('--rscript',default='Rscript')
 a=p.parse_args();package=a.package_root.resolve();inputs=a.input_root.resolve();cache=a.source_cache.resolve()
 if package==inputs or package in inputs.parents:raise RuntimeError('Input directory must be outside the public package')
 if inputs.exists() and any(inputs.iterdir()):raise RuntimeError('Refuse to overwrite existing input directory')
 inputs.mkdir(parents=True,exist_ok=True);cache.mkdir(parents=True,exist_ok=True)
 env=os.environ.copy()
 if os.name=='nt':
  for k in list(env):
   if k.startswith('LC_') or k in ['LANG','LANGUAGE']:env.pop(k)
  env['LANG']='English_United States.utf8'
 cmd=[a.rscript,'--vanilla',str(package/'R/00_acquire_translated_inputs.R'),str(package),str(cache),str(a.translated_cache.resolve()) if a.translated_cache else '-',str(inputs)]
 with (inputs/'source-acquisition.log').open('wb') as f:subprocess.run(cmd,cwd=inputs,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 (inputs/'official_xpt').mkdir();(inputs/'official_docs').mkdir()
 components={x['Component']:x for x in rows(package/'config/nhanes-components.csv')}
 docs={x['component']:x for x in rows(package/'config/official-codebook-urls.csv')}
 for x in rows(package/'config/source-inputs.csv'):
  dest=inputs/x['relative_path']
  if x['kind']=='official_xpt':shutil.copyfile(cache/components[x['component']]['Local_cache_filename'],dest)
  elif x['kind']=='official_codebook':
   req=urllib.request.Request(docs[x['component']]['url'],headers={'User-Agent':'NHANES-reproducibility-validation/1.1.0'})
   with urllib.request.urlopen(req,timeout=60) as resp:dest.write_bytes(resp.read())
   if sha(dest)!=x['sha256'].upper():raise RuntimeError('Official codebook identity mismatch: '+x['component'])
  if not dest.is_file():raise RuntimeError('Missing input: '+x['relative_path'])
 (inputs/'acquisition-command.json').write_text(json.dumps(dict(command=cmd,public_sources_verified=True),indent=2)+'\n')
 print('66 translated components, 16 XPTs and 24 codebooks acquired; run_pipeline.py verifies all input identities')
if __name__=='__main__':main()
