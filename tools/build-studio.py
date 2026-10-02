#!/usr/bin/env python3
"""Prepara il Quaderno di Chimica Generale (studio/) come pagina singola pubblicabile: CSS incorporato, moduli JS separati."""
import re, shutil, sys, pathlib
root = pathlib.Path(__file__).resolve().parent.parent / 'studio'
out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else root.parent / 'dist-studio')
if out.exists():
    shutil.rmtree(out)
shutil.copytree(root / 'js', out / 'js')
html = (root / 'index.html').read_text()
css = (root / 'css' / 'studio.css').read_text()
head = html[html.index('<head>') + 6:html.index('</head>')]
body = html[html.index('<body>') + 6:html.index('</body>')]
head = re.sub(r'<meta charset[^>]*>\n', '', head)
head = re.sub(r'<meta name="viewport"[^>]*>\n', '', head)
head = head.replace('<link rel="stylesheet" href="css/studio.css">', '<style>\n' + css + '\n</style>')
(out / 'index.html').write_text(head.strip() + '\n' + body.strip() + '\n')
print('\n'.join(sorted(str(p.relative_to(out)) for p in (out / 'js').rglob('*.js'))))
