#%% FVA of C1 assimilation strategies for co-substrate free growth coupled protocatechuate production 
"""
Flux variability analysis (FVA) of C1 assimilation strategies for 
protocatechuate production from ferulate in *Pseudomonas putida* KT2440 
using the curated iJN1445 model.

Key features:
- Applies to eight model variants with different formaldehyde assimilation strategies.
- Simulates *pcaHG* knockout to prevent protocatechuate degradation.
- Enforces protocatechuate export flux (4.5 mmol/gCDW/h).
- Uses hydrogen (H₂) as the inorganic electron donor (uptake rate: -100 mmol/gCDW/h).
- Performs FVA at 90% of the maximum biomass yield.
- Reaction metadata and flux variability results saved per model.

Output:
- FVA results saved to `fva_ferulate_C1/iJN1445_ferulate_C1_fva_results.xlsx`.

Author: Daniel Bergen  
Reference: *Co-Substrate Free Valorisation of Lignin Monomers by Assimilation of C1 and C2 By-Products*  
Daniel Bergen¹², Òscar Puiggené³, Esteban Marcellin¹⁴, Robert E. Speight²⁴, Pablo I. Nikel³, Birgitta E. Ebert¹⁵*

1 Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, Australia  
2 Advanced Engineering Biology Future Science Platform, CSIRO, Brisbane, Australia  
3 The Novo Nordisk Foundation Center for Biosustainability, Technical University of Denmark  
4 ARC Centre of Excellence in Synthetic Biology (CoESB)  
5 Food and Beverage Accelerator (FaBA), The University of Queensland
"""
#%%
import cobra
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
from pathlib import Path

Path("fva_ferulate_C1").mkdir(exist_ok=True)

#%%

models_list = [
    'iJN1445',
    'iJN1445_SACA',
    'iJN1445_RuMP',
    'iJN1445_RGP',
    'iJN1445_HOM',
    'iJN1445_SAL',
    'iJN1445_RGSALc',
    'iJN1445_all'
    ]

h2_uptake_rate = -100
fraction_of_optimum = 0.9  # 90 % of optimum biomass formation rate allowing for flux variability analysis

with pd.ExcelWriter("fva_ferulate_C1/iJN1445_ferulate_C1_fva_results.xlsx", engine="openpyxl") as writer:
    for id in models_list:
        model = cobra.io.read_sbml_model('C1_models/'+f'{id}.xml')
        model.solver = 'glpk'
        # set ferulate uptake rate to 4.5 mmol/gCDW/h
        model.reactions.EX_fer_e.lower_bound = -4.5
        # simulating pcaHG deletion preventing protocatechuate degradation
        model.genes.PP_4655.knock_out()
        model.genes.PP_4656.knock_out()

        # enforcing protocatechuate production
        model.reactions.EX_34dhbz_e.bounds = (4.5,4.5)

        # set H2 uptake rate
        model.reactions.EX_h2_e.lower_bound = h2_uptake_rate
        model.reactions.EX_pt_e.lower_bound = 0

        # Perform flux variability analysis
        fva_result = flux_variability_analysis(model,fraction_of_optimum=fraction_of_optimum)

        rxn_info = pd.DataFrame({
            'reaction id': [rxn.id for rxn in model.reactions],
            'reaction formula': [rxn.reaction for rxn in model.reactions],
            'gpr': [rxn.gene_reaction_rule for rxn in model.reactions],
            'lower bound': [rxn.lower_bound for rxn in model.reactions],
            'upper bound': [rxn.upper_bound for rxn in model.reactions],
        }).set_index('reaction id')

        fva_result = pd.concat([rxn_info, fva_result], axis=1)

        fva_result.to_excel(writer, sheet_name=id)

