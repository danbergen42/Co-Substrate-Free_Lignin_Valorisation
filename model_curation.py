#%% model curation
"""
Model curation script for Pseudomonas putida KT2440 genome-scale model (iJN1463),
used to generate a reduced and curated model (iJN1445) for C1 and lignin monomer
assimilation studies.

Key modifications:
- Defined medium composition based on deBont minimal salts medium (Hartmans & deBont, 1992).
- Removed Ni2+ from biomass reactions and adjusted reaction bounds to prevent undesired growth modes.
- Deleted TOL plasmid-associated genes (pWW0).
- Curated reaction reversibility and removed redundant lipoamide-dependent complexes.
- Pruned unused metabolites and set the biomass objective.

Optional modifications:
- Simulations were performed to predict metabolic capabilities for co-substrate free lignin valorisation in P. putida EM42
- The model iJN1463 is based on KT2440 - metabolically KT2440 and EM42 are identical based on their metabolically relevant gene sets
- To fully capture the actual genome of EM42, the genes, PP_1339 (encoding for 'ALAALAr') and PP_4545 (encoding for 'KAS15' and 'OGMEACPS'), can be considered for deletion in the iJN1445 model for future simulations of P. putida EM42
- However, those genes have isoforms (PP_4346 and PP_4379, respectively) in the KT2440 and the EM42 genomes - hence, deletion of PP_1339 and PP_4545 won't change the reaction set in the iJN1445 model

Output:
- Curated model iJN1445 saved as SBML (.xml).

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

#%% Define model medium, based on deBont minimal salts media without carbon source (Hartmans & deBont, 1992)

model_medium = {
        'EX_ca2_e': 10.0,
        'EX_cl_e': 10.0,
        'EX_co2_e': 100.0,
        'EX_cobalt2_e': 10.0,
        'EX_cu2_e': 10.0,
        'EX_fe2_e': 10.0,
        'EX_h2o_e': 100.0,
        'EX_h_e': 100.0,
        'EX_hco3_e': 10.0,
        'EX_k_e': 10.0,
        'EX_mg2_e': 10.0,
        'EX_mn2_e': 10.0,
        'EX_mobd_e': 10.0,
        'EX_na1_e': 10.0,
        'EX_nh4_e': 30.0,
        'EX_ni2_e': 0.0,
        'EX_o2_e': 100.0,
        'EX_pi_e': 10.0,
        'EX_sel_e': 10.0,
        'EX_so4_e': 10.0,
        'EX_tungs_e': 10.0,
        'EX_zn2_e': 10.0,
        'EX_acmtsoxin_e': 1000.0,
        'EX_acpptrn_e': 1000.0,
        'EX_d2one_e': 1000.0,
        'EX_d3one_e': 1000.0,
        'EX_d4one_e': 1000.0,
        'EX_mtsoxin_e': 1000.0,
        'EX_n2one_e': 1000.0,
        'EX_pptrn_e': 1000.0,
        'EX_und2one_e': 1000.0,
        }

#%% Model curation
iJN1463 = cobra.io.load_model('iJN1463')
# Update model medium
iJN1463.medium = model_medium

# Remove ni2 from biomass equations, based on deBont media (Hartmans & deBont, 1992)
original_coefficient = iJN1463.reactions.BIOMASS_KT2440_WT3._metabolites[iJN1463.metabolites.get_by_id('ni2_c')]
iJN1463.reactions.BIOMASS_KT2440_WT3.subtract_metabolites({iJN1463.metabolites.ni2_c: original_coefficient})
original_coefficient = iJN1463.reactions.BIOMASS_KT2440_Core2._metabolites[iJN1463.metabolites.get_by_id('ni2_c')]
iJN1463.reactions.BIOMASS_KT2440_Core2.subtract_metabolites({iJN1463.metabolites.ni2_c: original_coefficient})

# change AKGDa and PDH to irreversible, preventing autotrophic growth
iJN1463.reactions.AKGDa.lower_bound = 0
iJN1463.reactions.PDHa.lower_bound = 0

# full removal of pWW0 genes originally from TOL plasmid from Pseudomonas putida mt-2
pWW0_list = []
for gene in iJN1463.genes:
    if 'pWW0' in str(gene):
        pWW0_list += [gene]
cobra.manipulation.delete.remove_genes(iJN1463,pWW0_list)

# adjust reversibility of glyA and hom, based on data from SEM group
iJN1463.reactions.GHMT2r.bounds = (-1000,1000)
iJN1463.reactions.HSDy.bounds = (-1000,0)

# disable reversibility of GARFT (purN) preventing formatotrophic growth
iJN1463.reactions.GARFT.bounds = (0,1000)

# disable formaldehyde dismutase to prevent thermodynamically infeasible loops
iJN1463.reactions.FALDM.bounds = (0,0)

# remove lipoamide dependent complexes of PDH and AKGD for simplicity (PDH, AKGD still present)
iJN1463.reactions.PDHa.remove_from_model()
iJN1463.reactions.PDHbr.remove_from_model()
iJN1463.reactions.PDHcr.remove_from_model()
iJN1463.reactions.AKGDa.remove_from_model()
iJN1463.reactions.AKGDb_copy1.remove_from_model()
iJN1463.reactions.AKGDb_copy2.remove_from_model()

# remove unused metabolites for simplicity
iJN1463, deleted_metabolites = cobra.manipulation.delete.prune_unused_metabolites(iJN1463)

# set objective function with new biomass equation (no Ni2)
objective = 'BIOMASS_KT2440_WT3'
iJN1463.objective = objective

#%% Export model (iJN1445: reduced iJN1463 with protocatechuate production and inorganic electron donor supply, reduction based on removal of pWW0 genes)
iJN1463.id = 'iJN1445'
cobra.io.write_sbml_model(iJN1463, "iJN1445.xml")
