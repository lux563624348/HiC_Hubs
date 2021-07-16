# NOTE: For Paper Review, please follow the instruction below:
Latest updated on 07/04/2021,

# Comprehensive Network Analysis for HiC


<br><br>
<img src="image/Hub_Myb.PNG" width="800">
<br><br>


- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Example of Running (Demo)](#Example_Running)
- [License](#license)

## Overview
This module is a Python package containing tool for network analysis of HiC data.
It starts from HiC Interaction pairs, then generating network and clustering. Finally ranking all clusters by their interaction change.

## System Requirements
### Hardware Requirements

This package requires only a standard computer with enough RAM to support the in-memory operations.

### Software Requirements

HicHub mainly depends on the Python scientific stack.

```
python >=3
pandas
numpy
pybedtools
python-igraph
scipy
```

## Installation Guide
Recommend to use bioconda for installing.
Create the environment from the environment_hichub.yml(Can be found in this repository) file:
```
conda env create -f environment_hichub.yml
python3 -m pip install hichub --user
python3 -m pip install numpy pandas pybedtools python-igraph scipy
```
```
https://bioconda.github.io/user/install.html
```

## Example of Running (Demo)
After installation, type hichub in your command line will return the following UI:
```
welcome
The python env is: 3.7.10 | packaged by conda-forge | (default, Feb 19 2021, 16:07:37)
[GCC 9.3.0]
usage: hichub [-h] [-v] {test,diff,convert} ...
hichub -- A toolset to detect and analyze differential Hubs
positional arguments:
  {test,diff,convert}  sub-command help
    test               Output for test parser
    diff               Parser for diff hubs
    convert            Convert multi .hic to txt format --> Format should be:
                       #chr bin1 bin2 Count
optional arguments:
  -h, --help           show this help message and exit
  -v, --version        show program's version number and exit
```

diff:
```
usage: hichub diff [-h] -i <file> -f <str> -b <str> -r <int> [-p <float>]
                   [-t <int>]
hichub diff: error: the following arguments are required: -i/--in, -f/--foreground_name, -b/--background_name, -r/--resolution

Input Format: HiC Interaction in txt format.
Example of test data can be found in ~/test_data

Output can be found at working directory: 
(Demo output is: /HiC_Hubs/python_package/tests/2704_DKO_na_Diff_hub.txt  or 3106_WT_na_Diff_hub.txt)

```

convert:
%% Convert .hic to required input format
```
usage: hichub convert [-h] -i <file> [-n <str>] -f <str> -l <str> -r <int>
hichub convert: error: the following arguments are required: -i/--in, -f/--file_name, -l/--file_label, -r/--resolution
#Example:
#chr	bin1	bin2	Cond1	Cond2
10	3000000	3010000	100	200
```


```
3106_WT_na_Diff_hub.txt
0	1	2	hub_name	Num_vertices	pvalue
chr10	20930000	21060000	chr10:20930000-21060000	11	7.88966007260005e-09
chr10	19590000	19720000	chr10:19590000-19720000	11	7.809766623341443e-05
chr10	80210000	80340000	chr10:80210000-80340000	11	9.520611432439225e-05
chr10	95890000	96030000	chr10:95890000-96030000	14	0.00015075762147303865
```

## Built With

## Contributing

Please read (https:xx) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

## Authors

* *Xiang Li *Initial work* 


## License

#This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments


