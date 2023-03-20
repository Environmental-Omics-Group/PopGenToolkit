[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.7](https://img.shields.io/badge/python-3.7-green.svg)](https://www.python.org/dev/peps/pep-0537/)
[![R 3.6](https://img.shields.io/badge/R-3.6-red.svg)](https://cloud.r-project.org)


# PopGenToolkit

PopGenToolkit is a workflow to analyse population genomics data particularly for environmental studies.

PopGenToolkit is distributed under the GNU General Public License v3.0. The dependencies are distributed under their own licenses. 


## Dependencies

**Linked repositories**

- [vcflib](https://github.com/vcflib/vcflib) (with our [revision](https://github.com/Environmental-Omics-Group/vcflib))
- [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/)
- [HMM-detection-of-genomic-islands](https://github.com/marqueda/HMM-detection-of-genomic-islands)
- [TAFT](https://pubmed.ncbi.nlm.nih.gov/20943011/)

**Software and Packages**

- VCFtools >= 0.1.17
- Python >= 3.7
- R >= 3.6
- ipython >= 7.13.0
  - (Optional) Jupyter Notebook >= 6.0.3
  - (Optional) JupyterLab >= 2.1
- Python Packages:
  - NumPy >= 1.17.4
  - SciPy >= 1.4.1
  - kagami >= 3.0.8

Lower versions may work but have not been tested.


## Installation

```bash
git clone --recursive https://github.com/Environmental-Omics-Group/PopGenToolkit.git
```

Rebuild linked repositories if needed.


## Citation

If you use PopGenToolkit in a publication, we would appreciate citations: 

- *Chaturvedi, A., Zhou, J., Raeymaekers, J.A.M. et al. Extensive standing genetic variation from a small number of founders enables rapid adaptation in Daphnia. Nat Commun 12, 4306 (2021). https://doi.org/10.1038/s41467-021-24581-z*
