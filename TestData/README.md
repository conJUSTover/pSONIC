# Test Datasets for running pSONIC.    

We have provided two test datasets to test your installation of pSONIC. The results from the cotton dataset is originally described in the manuscript and is based on 4 species in the Gossypieae, including one allotetraploid. We provide these data in two separate ways - one with the tetraploid split into its respecive subgenomes, and one treating it as a single species. In the manuscript, we show that there is little different in these appraoches, and we recommend that users do not split polyploids into subgenomes in their own analyses.       

The wheat dataset is provided here as a more complex dataset, and includes one allohexaploid (bread wheat), two allotetraploids (durum wheat, wild emmer wheat), and three diploids relatives.       

## Runtimes    
Both datasets were tested on a MacOSx 3.4 GHz Intel Core i7 with 16GB DDR3 memory. Both tests were run with 8 threads. 

| Dataset | # Species | # Polyploids | # Collinear Groups Identified by MCScanX | # Collinear Groups passed pSONIC filtering | Time to Run pSONIC |    
|----|------|---|---|-----| -----|     
|Cotton (No Split)| 4 | 1 (1 allotetraploid) |24,502|2,128| 56 minutes|   
|Wheat| 6 | 3 (1 allohexaploid, 2 allotetraploids) |31,024|17,161| 4 hours, 5 minutes|    


## Results     
    
### Cotton (_Gossypieae_)   

||||||   
|----|---|---|---|---|
|Species| _Gossypioides kirkii_| _Gossypium arboreum_| _Gossypium hirsutum_| _Gossypium raimondii_ |  
|Relative Ploidy|1|1|2 (tetraploid)|1|
|# Genes in genome|36,669|37,972|65,636|37,223|
|# Genes Placed by OrthoFinder|26,671|33,401|56,523|34,271|
|# Genes placed by pSONIC|27,655|31,298|58,623|32,777|



### Wheat (_Triticeae_)    

||||||||   
|----|---|---|---|---|---|---|   
|Species| _Aegilops tauschii_| _Hordeum vulgare_| _Triticum aesticum_| _Triticum dicoccoides_ | _Triticum turgidum_ | _Triticum urartu_ |    
|Relative Ploidy|1|1|3 (hexaploid)|2 (tetraploid)|2 (tetraploid)|1|
|# Genes in genome|38,775|35,503|105,200|106,759|63,946|38,050|
|# Genes Placed by OrthoFinder|35,869|28,918|94,198|59,176|58,801|30,774|
|# Genes placed by pSONIC|29,600|22,431|94,461|57,662|57,029|23,645|
