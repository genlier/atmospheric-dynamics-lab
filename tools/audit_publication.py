"""Check the publication tree without printing potentially sensitive matches."""
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
BLOCKED={'.docx','.pdf','.mat','.nc','.grib','.pem','.key','.p12','.asv','.log'}
PATTERNS=[r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]',r'/Users/[^/\s]+/',
          r'[A-Za-z0-9_.+\-]+@[A-Za-z0-9.\-]+\.[A-Za-z]{2,}',
          r'gh[pousr]_[A-Za-z0-9]{20,}',r'github_pat_[A-Za-z0-9_]{20,}',
          r'AKIA[0-9A-Z]{16}',r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',
          r'(?:姓名|学号|手机号|身份证号)\s*[:：]\s*\S+']
def audit():
    problems=[]; count=0
    for path in ROOT.rglob('*'):
        if not path.is_file() or '.git' in path.parts or '__pycache__' in path.parts:
            continue
        if path.relative_to(ROOT).parts[0]=='outputs':continue
        count+=1
        if path.suffix.lower() in BLOCKED or path.name in {'.env','.cdsapirc'} or path.name.startswith('.env.'):
            problems.append(str(path.relative_to(ROOT))); continue
        if path.suffix.lower()=='.mp4' and path.parent==ROOT/'docs'/'animations':
            with path.open('rb') as video: header=video.read(12)
            if header[4:8]!=b'ftyp': problems.append(str(path.relative_to(ROOT)))
            continue
        try: content=path.read_text(encoding='utf-8')
        except UnicodeDecodeError:
            problems.append(str(path.relative_to(ROOT))); continue
        if path.name=='audit_publication.py':continue
        if any(re.search(p,content) for p in PATTERNS):problems.append(str(path.relative_to(ROOT)))
    if problems:
        raise SystemExit('Publication audit failed in: '+', '.join(problems))
    print(f'Publication audit passed: {count} allowlisted text files.')
if __name__=='__main__':audit()
