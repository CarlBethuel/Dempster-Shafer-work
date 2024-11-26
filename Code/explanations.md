The code is organized into three scripts: 
- Sampling_method script
Creates a train (70%) and test (30%) sample to apply fusion methods.
It is based on two functions: sample creation and sample preparation. 

- A DST_method_4sources script
Applies DST for four input sources, with mass assignment from Kappa
- A Voting_strategies script
Applies a voting strategy approach for two weight vectors: a vector where the weight of each classification is equal to 1 and a vector where the weight of each classification corresponds to the Kappa of that classification.
It is based on four functions.

The logic used in the script is applied for two land use classes (Oil Palm and No Oil Palm) and 4 input sources (see [Features] section). The code can be modified depending on the application. 
