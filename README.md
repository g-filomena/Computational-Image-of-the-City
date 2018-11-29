# The Computational Image of the City

**A tool for extracting the Computational Image of the City from geospatial datasets**

In this repository a set of functions to extract the Image of the City "[The Image of The City, Lynch, K., 1960. The Image Of The City, Cambridge, MA: MIT Press](https://mitpress.mit.edu/books/image-city)." from open and freely available geospatialdatasets is presented.

The tool is written in Python and rely on a series of libraries, as [Geopandas](http://geopandas.org), [OSMNx](https://osmnx.readthedocs.io/en/stable/), and relative dependencies.
The functions are divided in three scripts, accordingly to their purposes:
* street_network_functions.py
* landmarks_functions.py
* utilities.py 

The first set of functions support the extraction of *Nodes, Paths and Districts*, besides basic street network operations (loading or downloading network, cleaning, simplifcation, etc.). The second group allows to extract computational landmarks. The pipeline and the usage of the functions are shown in the notebook *1_Nodes_ paths_districts.ipynb* and *2_Landmarks.ipynb*.  While the first notebook does not require the user to load data locally stored, the second one does. In addition, part of the analysis for the extraction of landmarks has to be run in ArcGisPro or ArcScene environment (3d sight lines) and it's therefore not completely automatised. 

