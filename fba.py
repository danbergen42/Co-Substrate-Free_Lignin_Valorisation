#%% Initital FBA for co-substrate free valorisation of lignin monomers
"""
Flux balance analysis (FBA) using the curated genome-scale model iJN1445 of 
*Pseudomonas putida* KT2440 to assess carbon source utilization and 
protocatechuate production from lignin monomers.

Key features:
- Acetate serves as control carbon source for growth without protocatechuate enforcement (especially relevant for simulations with pcaHG deletion due to same cofactor demand 'til acetyl-CoA cleavage from aromatic compounds)
- p-Coumarate and ferulate simulations include *pcaHG* (PP_4655 and PP_4656) knockout to prevent protocatechuate degradation.
- Protocatechuate production enforced via fixed export flux (4.5 mmol/gCDW/h).
- Substrate uptake rates fixed at 4.5 mmol/gCDW/h.
- Reaction metadata and flux distributions exported to Excel.

Output:
- FBA flux distributions saved to `fba/iJN1445_fba_results.xlsx`.

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
import pandas as pd
import numpy as np
from pathlib import Path

Path("fba").mkdir(exist_ok=True)

model = cobra.io.read_sbml_model('iJN1445.xml')
model.solver = 'glpk'
#%% Perform FBA with the curated iJN1445 model for p-coumarate and ferulate, including acetate as control condition

fba_results = {
    'wt_acetate': {},
    'wt_p-coumarate': {},
    'wt_ferulate': {},
    'dpcaHG_acetate': {},
    'dpcaHG_p-coumarate': {},
    'dpcaHG_ferulate': {}
}
uptake_reactions = {
    'acetate': 'EX_ac_e',
    'p-coumarate':'EX_T4hcinnm_e',
    'ferulate': 'EX_fer_e'
}
with pd.ExcelWriter('fba/iJN1445_fba_results.xlsx') as writer:
    for substrate in fba_results.keys():
        if substrate.split('_')[1] != 'acetate' and 'dpcaHG' in substrate:
            # simulating pcaHG deletion (protocatechuate 3,4-dioxygenase) preventing protocatechuate degradation from lignin monomers
            model.genes.PP_4655.knock_out()
            model.genes.PP_4656.knock_out()
            # enforcing protocatechuate production
            model.reactions.EX_34dhbz_e.bounds = (4.5,4.5)
        model.reactions.get_by_id(uptake_reactions[substrate.split('_')[1]]).lower_bound = -4.5

        rxn_info = pd.DataFrame({
            'reaction id': [rxn.id for rxn in model.reactions],
            'reaction formula': [rxn.reaction for rxn in model.reactions],
            'gpr': [rxn.gene_reaction_rule for rxn in model.reactions],
            'lower bound': [rxn.lower_bound for rxn in model.reactions],
            'upper bound': [rxn.upper_bound for rxn in model.reactions],
        }).set_index('reaction id')

        solution = model.optimize()
        fba_results[substrate] = pd.concat([rxn_info, solution.fluxes], axis=1)
        fba_results[substrate].to_excel(writer, sheet_name=substrate)
        model.reactions.get_by_id(uptake_reactions[substrate.split('_')[1]]).lower_bound = 0.0

#%% Robustness analysis for acetate, p-coumarate, and ferulate uptake rates

robustness_results = {
    'acetate': {},
    'p-coumarate': {},
    'ferulate': {},
}
uptake_reactions = {
    'acetate': 'EX_ac_e',
    'p-coumarate':'EX_T4hcinnm_e',
    'ferulate': 'EX_fer_e'
}

# set uptake rate ranges for robustness analysis: 2.0 to 6.0 mmol/gCDW/h in steps of 0.5 mmol/gCDW/h
substrate_upt_ranges = {
    'acetate': [2.0, 6.5, 0.5],
    'p-coumarate': [2.0, 6.5, 0.5],
    'ferulate': [2.0, 6.5, 0.5]
}
with pd.ExcelWriter('fba/iJN1445_upt_robustness.xlsx') as writer:
    for substrate in robustness_results.keys():
        if 'dpcaHG' not in substrate:
            for uptake_rate in np.arange(*substrate_upt_ranges[substrate]):
                model.reactions.get_by_id(uptake_reactions[substrate]).lower_bound = -uptake_rate
                solution = model.optimize()
                robustness_results[substrate][f'uptake_{uptake_rate}'] = solution.fluxes
                model.reactions.get_by_id(uptake_reactions[substrate]).lower_bound = 0.0
            # Convert dict of Series to DataFrame (columns = uptake rates)
            flux_df = pd.DataFrame(robustness_results[substrate])
            flux_df.to_excel(writer, sheet_name=substrate)

