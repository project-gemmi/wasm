#!/usr/bin/env python3

from string import Template

with open('template.html') as f:
    tpl = Template(f.read())
with open('pdb2cif.html', 'w') as f:
    f.write(tpl.safe_substitute({'FROM': 'PDB',
                                 'TO': 'mmCIF',
                                 'FUNC': 'pdb2cif',
                                 'BLOBTYPE': 'text/plain',
                                 'EXT': 'cif'}))
with open('cif2pdb.html', 'w') as f:
    f.write(tpl.safe_substitute({'FROM': 'mmCIF',
                                 'TO': 'PDB',
                                 'FUNC': 'cif2pdb',
                                 'BLOBTYPE': 'text/plain',
                                 'EXT': 'pdb'}))
with open('cif2mtz.html', 'w') as f:
    f.write(tpl.safe_substitute({'FROM': 'mmCIF',
                                 'TO': 'MTZ',
                                 'FUNC': 'cif2mtz',
                                 'BLOBTYPE': 'application/octet-stream',
                                 'EXT': 'mtz'}))
with open('mtz2cif.html', 'w') as f:
    f.write(tpl.safe_substitute({'FROM': 'MTZ',
                                 'TO': 'mmCIF',
                                 'FUNC': 'mtz2cif',
                                 'BLOBTYPE': 'text/plain',
                                 'EXT': 'cif'}))
