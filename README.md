# The Computational Image of the City

**A tool for extracting the Computational Image of the City from geospatial datasets**

In this repository a set of functions to extract the Image of the City "[The Image of The City, Lynch, K., 1960. The Image Of The City, Cambridge, MA: MIT Press](https://mitpress.mit.edu/books/image-city)." from open and freely available geospatialdatasets is presented.

The tool is written in Python and rely on a series of libraries, as [Geopandas](http://geopandas.org), [OSMNx](https://osmnx.readthedocs.io/en/stable/), and relative dependencies.

The functions used are divide in three scripts:
street_network_functions.py
landmarks_functions.py
utilities.py 

The first set of functions support the extraction of *Nodes, Paths and Districts* besides basic street network operations (loading or downloading network, cleaning, simplifcation, etc.). The pipeline and the usage of the functions are shown in the notebook  to extract *Nodes, Paths and Districts* employing network techniques. 
