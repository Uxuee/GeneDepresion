"""Read-only consolidation checks; --write records only this audit's JSON result."""
import argparse
import csv
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def readcsv(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true');args=parser.parse_args()
    errors=[]; warnings=[];checks={}; P=ROOT/'provenance/consolidation'
    baseline=json.loads((P/'frozen_baseline.json').read_text())
    count=0; shared_differences=[]
    for key,h in baseline.items():
        scope,rel=key.split('/',1)
        if rel=='provenance/canonical_pipeline.json':continue
        if scope=='rra' and 'original/'+rel in baseline:
            if h!=baseline['original/'+rel]:shared_differences.append(rel)
            # The destination's pre-integration baseline has priority. Never rewrite
            # original upstream files to match another checkout's pre-existing bytes.
            continue
        p=ROOT/rel
        if not p.is_file():errors.append('Missing frozen source: '+key)
        elif digest(p)!=h:errors.append('Frozen checksum differs: '+key)
        count+=1
    checks['frozen_source_comparisons']=count
    checks['preexisting_shared_checkout_byte_differences']=len(shared_differences)
    ev=readcsv(P/'external_validation_inventory.csv')
    for row in ev:
        p=ROOT/row['path']
        if not p.is_file() or digest(p)!=row['sha256']:errors.append('External validation altered/missing: '+row['path'])
    checks['external_validation_preserved']=len(ev)
    copies=json.loads((P/'local_stage_copies.json').read_text())
    for row in copies:
        p=ROOT/row['path']
        if not p.is_file() or digest(p)!=row['sha256']:errors.append('Local stage copy differs: '+row['path'])
    checks['local_stage_copies_verified']=len(copies)
    registry=json.loads((ROOT/'provenance/canonical_pipeline.json').read_text())
    primary=registry['current_primary_ranking']; ranking=readcsv(ROOT/primary['path'])
    if digest(ROOT/primary['path'])!=primary['sha256']:errors.append('Primary ranking hash mismatch')
    if len(ranking)!=671:errors.append('Primary ranking must preserve 671 rows')
    tier=readcsv(ROOT/'outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv')
    alias={'Cingulate gyrus 25':'CG25'}
    a={(x['region'],x['gene_id']) for x in ranking};b={(alias.get(x['region_or_subregion'],x['region_or_subregion']),x['gene_id']) for x in tier}
    if a!=b or len(a)!=671:errors.append('Primary ranking candidate identity mismatch')
    groups=defaultdict(list)
    for x in ranking:groups[x['region']].append(int(x['consensus_rank']))
    if len(groups)!=6 or any(sorted(v)!=list(range(1,len(v)+1)) for v in groups.values()):errors.append('Regional ranks invalid')
    if primary['confidence_labels'] is not None:errors.append('RRA must not inherit weighted labels')
    checks['primary_regional_counts']={k:len(v) for k,v in groups.items()}
    ids=[s['stage_id'] for s in registry['stages']]
    if len(ids)!=len(set(ids)):errors.append('Duplicate registry stage identifiers')
    for s in registry['stages']:
        for dep in s['canonical_inputs']:
            if dep.startswith('stage:') and dep[6:] not in ids:errors.append('Unresolved stage dependency: '+dep)
        script=s['canonical_script']; scripts=script if isinstance(script,list) else [script] if script else []
        for item in scripts:
            if item.startswith('git:'):
                _,commit,path=item.split(':',2)
                r=subprocess.run(['git','cat-file','-e',commit+':'+path],cwd=ROOT,capture_output=True)
                if r.returncode:errors.append('Missing generator Git object: '+item)
            elif not (ROOT/item).exists():warnings.append('Registered script absent: '+item)
        if s['status']=='canonical':
            for out in s['canonical_outputs']:
                if not out.startswith('${') and not (ROOT/out).exists():errors.append('Canonical output missing: '+out)
    r=subprocess.run(['git','cat-file','-e','f700cb323a47b24825b11bc1bba33cc07e6430bb'],cwd=ROOT,capture_output=True)
    if r.returncode:errors.append('Frozen frontal LM Git blob unavailable')
    old=P/'canonical_pipeline.before.json'
    if digest(old)!=registry['migration_history'][0]['prior_registry_sha256']:errors.append('Prior registry changed')
    manifest=json.loads((ROOT/'provenance/canonical_outputs_manifest.json').read_text())
    for rec in manifest['records']:
        p=ROOT/rec['path']
        if rec['available_local']!=p.exists():errors.append('Output availability changed: '+rec['path'])
        if rec['sha256'] and (not p.is_file() or digest(p)!=rec['sha256']):errors.append('Manifest hash mismatch: '+rec['path'])
    # Directory members already covered by baseline or local-copy checks where frozen.
    checks['output_manifest_records']=len(manifest['records'])
    for row in json.loads((P/'local_input_copies.json').read_text()):
        if digest(ROOT/row['path'])!=row['sha256']:errors.append('Copied ML input differs: '+row['path'])
    docs=['PIPELINE.md','METHODS.md','RESULTS.md','LIMITATIONS.md','REPRODUCIBILITY.md','CHANGELOG_RESEARCH.md','ROADMAP.md','STEP3_PROVENANCE_RECONCILIATION.md','CONSOLIDATION_PHASE2.md']
    for name in docs:
        p=ROOT/'docs'/name
        if not p.is_file():errors.append('Missing current guide '+name);continue
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if not target.startswith(('http:','https:','#')) and not (p.parent/target.split('#')[0]).exists():errors.append('Broken guide link: '+name+' -> '+target)
    # Inspection only: do not import scientific training modules.
    common=(ROOT/'code/ml_rra/common.py').read_text()
    expr_paths=sorted(set(re.findall(r'"(outputs/[^"\n]+\.csv)"',common)))
    missing=[p for p in expr_paths if not (ROOT/p).is_file()]
    checks['ml_expression_input_paths_checked']=len(expr_paths);checks['missing_ml_expression_inputs']=missing
    if missing:warnings.append('Historical expression inputs absent; no ML rerun certified: '+str(len(missing)))
    checks['python_executable']=sys.executable
    checks['python_dependencies']={name:importlib.util.find_spec(name) is not None for name in ['numpy','pandas','scipy','sklearn','matplotlib']}
    if not all(checks['python_dependencies'].values()):warnings.append('Some scientific Python dependencies unavailable in audit interpreter')
    warnings.extend(['Frozen LM direct generator recovery remains incomplete; current Step3 scripts are separate branches.','External-reference redistribution review remains pending; existing files preserved untracked.','R packages and cross-platform end-to-end execution were not revalidated; no analyses run.'])
    result={'date':'2026-11-21','ok':not errors,'errors':errors,'warnings':warnings,'checks':checks}
    if args.write:(P/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2));return bool(errors)
if __name__=='__main__':raise SystemExit(main())
