#!/usr/bin/env python3

from string import Template

with open('template.html') as f:
    tpl = Template(f.read())
with open('pdb2cif.html', 'w') as f:
    f.write(tpl.safe_substitute({'FROM': 'PDB',
                                 'TO': 'mmCIF',
                                 'FUNC': 'pdb2cif',
                                 'EXT': 'cif'}))
with open('cif2pdb.html', 'w') as f:
    f.write(tpl.safe_substitute({'FROM': 'mmCIF',
                                 'TO': 'PDB',
                                 'FUNC': 'cif2pdb',
                                 'EXT': 'pdb'}))
with open('mtz2cif.html', 'w') as f:
    f.write(tpl.safe_substitute({'FROM': 'MTZ',
                                 'TO': 'mmCIF',
                                 'FUNC': 'mtz2cif',
                                 'EXT': 'cif'}))
