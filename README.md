# Comprehensive Network Analysis for HiC

<br><br>
<img src="image/Hub_Myb.PNG" width="800">
<br><br>


One Paragraph of project description goes here

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites
python 3 
packages: optparse, pandas, numpy, pybedtools, igraph, scipy

```
Give examples
```

### Installing

Recommend to use bioconda for installing.

```
until finished
```


## Running the tests

EXAMPLE: python igraph_hub.py -i chr10_WT_na-DKO_na.bed -f WT_na -b DKO_na -r 10000 -d 0.5

Collect HiC Interaction in txt format, rank interaction change Hub. Input
Format should be: #chr        bin1    bin2    Cond1   Cond2

Options:
  -h, --help            show this help message and exit
  -i <file>, --in=<file>
                        Path to Input HiC file in txt format
  -f <str>, --foreground_name=<str>
                        Name of condition as foreground.
  -b <str>, --background_name=<str>
                        Name of condition as background.
  -r <int>, --resolution=<int>
                        Resolution of HiC txt
  -d <float>, --filtered_density=<float>
                        Density cutoff for hub shriking.

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* *Xiang Li *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)


## License

#This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

=======
# HiC_Hubs
>>>>>>> 1f58581f96c90768a77740df466c2833a4b970a1
