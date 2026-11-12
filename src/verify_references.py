"""Resolve the retained bibliography against Crossref's DOI registry."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,re,requests,time
ROOT=Path(__file__).resolve().parents[1]
def main():
    source=(ROOT/'archive/20261002_antes_simple_paper/artigo/referencias.bib').read_text(encoding='utf-8')
    entries=re.findall(r'@article\{(\w+),(.+?)(?=\n@article|\Z)',source,re.S)
    excluded={'hinrichsen2000','scheffer2009','boettiger2012','oliveira2005'}
    items=[(k,re.search(r'doi\s*=\s*\{([^}]+)',b).group(1),re.search(r'title\s*=\s*\{([^}]+)',b).group(1)) for k,b in entries if k not in excluded]
    def fetch(item):
        key,doi,title=item
        try:
            with requests.Session() as session:
                session.trust_env=False
                r=session.get('https://api.crossref.org/works/'+doi,timeout=40);r.raise_for_status();m=r.json()['message']
            return dict(key=key,doi=doi,local_title=title,registered_title=m.get('title'),metadata=m,status='verified')
        except Exception as exc:return dict(key=key,doi=doi,status='failed',error=str(exc))
    previous={r['key']:r for r in json.loads((ROOT/'qa/doi_verification.json').read_text(encoding='utf-8'))} if (ROOT/'qa/doi_verification.json').exists() else {}
    rows=[]
    for item in items:
        old=previous.get(item[0],{})
        if old.get('status')=='verified': rows.append(old); continue
        time.sleep(2)
        row=fetch(item);rows.append(row)
        print(item[0],row['status'],flush=True)
    (ROOT/'qa/doi_verification.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps([dict(key=r['key'],status=r['status'],title=r.get('registered_title')) for r in rows],ensure_ascii=True,indent=2))
if __name__=='__main__':main()
