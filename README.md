# The effects of sex-assortative mixing on the spread of respiratory-transmitted infections

### Paige Miller, John Drake

**(record any changes to this protocol at the end of this document)**

### Background: 

Multiple types of assortative mixing: 

* mixing by degree (popular with popular); similar mixing by risk level
* mixing by age (important enough that models without age-assortativity can't capture patterns in incidence very well)
* mixing by sex 

Sex assortativity among social interactions is common across human populations (Mossong 2008) and male-bias in TB reporting is found in almost every country in the world (WHO). Inclusion of differential contact rates by age (Arregui 2018) and sex (Dodd 2016) are required to accurately model TB incidence patterns. 

We're wondering whether sex assortativity is correlated with male-bias across countries in real data and whether we can observe sex-assortativity to drive male-bias in network simulations. 

### Study design: 

We will simulate outbreaks of TB on networks of varying sex-assortativity and measure male-bias.

We will find social network data (e.g., sociopatterns.com) or contact matrices data (Mossong 2008), calculate assortativity coefficient for sex, and compare with WHO estimated M:F incidence in the most recent year (I believe 2018). 

### Data sources: 

* TBD on network/contact matrix data to inform sex-assortativity across countries, it seems like it may be possible to mine facebook data? 
* WHO data is published online and has already been formatted for this analysis
* TBD Outbreak simulation code

### Analysis: 

* Correlation between male:female case reports and a measure of sex assortativity 
* Simulated male bias (ratio of male:female cases) of outbreaks on networks will be plotted against programmed sex assortativty of networks

### Checklist: 

* literature review on assortativity in social networks
* decide on data source for assortativity across countries
* decide on how to model assorted networks
* decide on how to model epidemics on networks
* pilot study on small number or size of networks

### Protocol changes:

* [bullet list of changes to protocols with reasons for each change]
