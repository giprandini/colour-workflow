
#### Infos

- AiiDA plugin for simple.x is called simple and is a Quantum Espresso plugin
- AiiDA plugin for simple_ip.x is called shirley and is an external plugin

#### How to install

plugin-input:

- simple.py in folder /orm/calculation/job/quantumespresso/ of AiiDA source code
- shirley.py in folder /orm/calculation/job/ of AiiDA source code

plugin-parser:

- simple.py in folder /parsers/plugins/quantumespresso/ of AiiDA source code
- shirley.py in folder /parsers/plugins/ of AiiDA source code

qe-kpoints:

- __init__.py and kpoints.py replace original files in /orm/data/array/ of AiiDA source code
