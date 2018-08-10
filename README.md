# Community Detection on Malarial Transmission Networks

 Malaria is transmitted to humans solely through the bites of infected mosquitoes. Nearly half of the world's population is at risk of malaria. In 2015 alone, there was an estimated 212 million cases and 429,000 deaths caused by malaria, posing a major threat to human health and economic growth around the world.Increases in intervention and control measures have resulted in a 29% reduction in malaria mortality rates globally since 2010 (WHO). Therefore, analyzing the efficacy and efficiency of intervention methods is a vital area of study for malaria elimination. One such area of study focuses on the targeted elimination of mosquitoes. These measures are effective because they aim to apply mosquito-specific pesticides in subsets of local populations most susceptible to infection. However, most, if not all, targeted intervention rely on detecting infection distribution relative to permanent human habitation. Yet, malarial infection is pervasive in transient and nomadic human populations around the world. Therefore, there is interest in targeted intervention aimed at the mosquito communities, in which intervention methods are deployed in vital locales of the mosquito habitat to increase mosquito mortality.

Mosquitoes do not move homogeneously over their habitat because they must move to fulfill certain biological needs. This heterogeneous movement, as well as spatial swarm segregation, results in the formation of mosquito communities. However, most methods of community detection were developed for undirected networks. These community detection algorithms rely solely on the topography of the network (e.g. measures of node degree and localized edge density) and so, are blind to the community structure imparted by directional movement. On the other hand, community detection methods for directed networks, such as the map framework algorithms, consider directional movement, but not geographic range. Many malarial interventions have limited geographic range and logistical challenges associated with deploying intervention campaigns over large land areas. 

In the context of the targeted intervention of mosquitos, both topographical and map-framework clustering algorithms have major limitations. **The goal of this project is to develop a geographic, flow-adjusted community detection algorithm, and simulate its feasibility as a tool for increasing the mosquito mortality of targeted intervention methods.**

This project is part of an ongoing effort between [Marshall Lab](https://www.marshalllab.com) at UC Berkeley and the [Institute for Health Metrics and Evaluation(IHME)](http://www.healthdata.org) at UW Medicine to model and analyze the spatial dynamics of malarial transmission.

More information about these efforts can be found here: 

https://chipdelmal.github.io/MoNeT/

https://marshalllab.github.io/MGDrivE/

We presented some preliminary results at the 2018 Netsci: International School and Conference on Network Science. 
Our poster, titled *"Network Analysis of Mosquito Habitats for Controlling Vector-Borne Pathogens"*, can be found here:  https://www.netsci2018.com/posters
