Numerical Macroeconomics

For the interactive computational environment:\
Tutorial1: please go to https://nextjournal.com/a/NRzh34EHQkJXnh6AVGySH?token=J6Rea4z51HyD2n95gJkjZa

The "Labour Matching FPI" is a dedicated exercise towards a standard labour matching model.

Main Files           						              
Step 1: perform PFI and export the results (Tut9_MainPartI_PFI.jl)
* input:     N/A	
* functions: Tut9_FunctionsI_Jacobianandco.jl           
Tut9_FunctionsII_Newtoniter.jl                
Tut9_FunctionsIII_chebyshevAppx.jl          
Tut9_FunctionsIV_chebyshevandco.jl        
* output:    PFIResults_zeta0?_time???.txt              
 PolicyFunctions_zeta0?_time???.png         
                          
* - Step 2: perform simulation on one, and then 100_000 time series                 *
*    Tut9_MainPartII_SimulateLMM.jl:				                *
*                                         input:       PFIResults_zeta0?_time???.txt                     *
*                                         functions: Tut9_FunctionsIII_chebyshevAppx.jl           *    
*                                                          Tut9_FunctionsIV_chebyshevandco.jl         *  
*                                         output:     KeyVar_timespan_zeta0?.png                     *
*			           MeanKeyVars_Welfare_dis_zeta0?ts50.png *
*                                                          MeanAllVar_dis_zeta0?ts50.png                  *
***********************************************************************************

****************************************************
*                          Some Metadata                           *
*  - working directory: where the script is located *
*  - figure path: working directory/figures              *
*  - Random.seed: 5086_2495262_20210331         * 
****************************************************

