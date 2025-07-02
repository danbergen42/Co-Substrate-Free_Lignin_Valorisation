#%% C1 assimilation model generation
"""
Flux balance analysis (FBA) for assessing formaldehyde assimilation pathways 
in *Pseudomonas putida* KT2440 using the curated iJN1445 model engineered for 
protocatechuate production from ferulate.

Key features:
- Simulates *pcaHG* knockout to enforce protocatechuate accumulation from ferulate.
- Adds hydrogenase and phosphite dehydrogenase modules for NADH generation 
  from molecular hydrogen and phosphite.
- Compares six formaldehyde assimilation strategies:
  1. Synthetic Acetyl-CoA (SACA) pathway
  2. Ribulose-5-monophosphate (RuMP) cycle
  3. Reductive Glycine Pathway (RGP)
  4. Homoserine cycle (HOM)
  5. Serine aldolase shunt (SAL)
  6. Combined Reductive Glycine–Serine Aldolase cycle (RGSAL_c)
- Each strategy simulated with and without H₂ or phosphite as electron donor.
- Flux distributions exported for all models and conditions.
- Models saved for downstream analysis.

Output:
- FBA results saved to `fba_ferulate_C1/iJN1445_ferulate_C1_fba_results.xlsx`.
- C1-assimilation models exported to `C1_models/*.xml`.

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
from cobra import Reaction, Metabolite
import pandas as pd

from pathlib import Path

Path("fba_ferulate_C1").mkdir(exist_ok=True)
Path("C1_models").mkdir(exist_ok=True)

#%% add hydrogenase and transport reactions as well as phosphite dehydrogenase and transport reactions as NADH generation module from inorganic electron donors

iJN1445 = cobra.io.read_sbml_model('iJN1445.xml')
iJN1445.solver = 'glpk'

# set ferulate uptake rate to 4.5 mmol/gCDW/h
iJN1445.reactions.EX_fer_e.lower_bound = -4.5
# simulating pcaHG deletion preventing protocatechuate degradation
iJN1445.genes.PP_4655.knock_out()
iJN1445.genes.PP_4656.knock_out()
# enforcing protocatechuate production
iJN1445.reactions.EX_34dhbz_e.bounds = (4.5,4.5)

# add hydrogenase and transport reactions as NADH generation module from molecular hydrogen
# add new metabolites
h2_e = Metabolite(
    'h2_e',
    formula='H2',
    name='molecular hydrogen',
    compartment='e',
    charge = 0
)
h2_e.annotation = {
    'bigg.metabolite': 'h2',
    'biocyc': 'META:HYDROGEN-MOLECULE',
    'chebi': ['CHEBI:13350', 'CHEBI:18276', 'CHEBI:25363', 'CHEBI:29294', 'CHEBI:29298', 'CHEBI:29299', 'CHEBI:5785'],
    'hmdb': ['HMDB01362'],
    'inchi_key': 'UFHFLCQGNIYNRP-UHFFFAOYSA-N',
    'kegg.compound': 'C00282',
    'metanetx.chemical': 'MNXM195',
    'sabiork': ['5030'],
    'seed.compound': ['cpd11640']}

h2_p = Metabolite(
    'h2_p',
    formula='H2',
    name='molecular hydrogen',
    compartment='p',
    charge = 0
)
h2_p.annotation = {
    'bigg.metabolite': 'h2',
    'biocyc': 'META:HYDROGEN-MOLECULE',
    'chebi': ['CHEBI:13350', 'CHEBI:18276', 'CHEBI:25363', 'CHEBI:29294', 'CHEBI:29298', 'CHEBI:29299', 'CHEBI:5785'],
    'hmdb': ['HMDB01362'],
    'inchi_key': 'UFHFLCQGNIYNRP-UHFFFAOYSA-N',
    'kegg.compound': 'C00282',
    'metanetx.chemical': 'MNXM195',
    'sabiork': ['5030'],
    'seed.compound': ['cpd11640']}

h2_c = Metabolite(
    'h2_c',
    formula='H2',
    name='molecular hydrogen',
    compartment='c',
    charge = 0
)
h2_c.annotation = {
    'bigg.metabolite': 'h2',
    'biocyc': 'META:HYDROGEN-MOLECULE',
    'chebi': ['CHEBI:13350', 'CHEBI:18276', 'CHEBI:25363', 'CHEBI:29294', 'CHEBI:29298', 'CHEBI:29299', 'CHEBI:5785'],
    'hmdb': ['HMDB01362'],
    'inchi_key': 'UFHFLCQGNIYNRP-UHFFFAOYSA-N',
    'kegg.compound': 'C00282',
    'metanetx.chemical': 'MNXM195',
    'sabiork': ['5030'],
    'seed.compound': ['cpd11640']}

nadh_c = iJN1445.metabolites.nadh_c
nad_c = iJN1445.metabolites.nad_c
h_c = iJN1445.metabolites.h_c

# C. necator SH hydrogenase
NAD_H2 = Reaction('NAD_H2')
NAD_H2.name = 'Soluble [NiFe] Hydrogenase from Cupriavidus necator'
NAD_H2.subsystem = ''
NAD_H2.lower_bound = -1000
NAD_H2.upper_bound = 1000
NAD_H2.add_metabolites({
    nad_c: -1.0,
    h2_c: -1.0,
    nadh_c: 1.0,
    h_c: 1.0
})
NAD_H2.gene_reaction_rule = '( hoxF and hoxU and hoxY and hoxH and hoxW and hoxI and hypA and hypB and hypF and hypC and hypD and hypE and hypX )'
NAD_H2.reaction

# H2 exchange
EX_h2_e = Reaction('EX_h2_e')
EX_h2_e.name = 'H2 transport, exchange'
EX_h2_e.subsystem = ''
EX_h2_e.lower_bound = -1000
EX_h2_e.upper_bound = 0
EX_h2_e.add_metabolites({
    h2_e: -1.0
})
EX_h2_e.reaction

# H2 transport
H2tpp1 = Reaction('H2tpp1')
H2tpp1.name = 'H2 transport, periplams'
H2tpp1.subsystem = ''
H2tpp1.lower_bound = -1000
H2tpp1.upper_bound = 1000
H2tpp1.add_metabolites({
    h2_e: -1.0,
    h2_p: 1.0
})
H2tpp1.reaction

# H2 transport
H2tpp2 = Reaction('H2tpp2')
H2tpp2.name = 'H2 transport, cytosol'
H2tpp2.subsystem = ''
H2tpp2.lower_bound = -1000
H2tpp2.upper_bound = 1000
H2tpp2.add_metabolites({
    h2_p: -1.0,
    h2_c: 1.0
})
H2tpp2.reaction

# enable NAD_H2
iJN1445.add_reactions([EX_h2_e, H2tpp1, H2tpp2, NAD_H2])

# add phosphite dehydrogenase and transport reactions as NADH generation module from phosphite
# define the necessary metabolites; annotation for pt missing
pt_e = Metabolite(
    'pt_e',
    formula='HO3P',
    name='phosphite',
    compartment='e',
    charge = -2
    )
# pt_e.annotation = {
#     }
pt_p = Metabolite(
    'pt_p',
    formula='HO3P',
    name='phosphite',
    compartment='p',
    charge = -2
    )
# pt_p.annotation = {
#     }
pt_c = Metabolite(
    'pt_c',
    formula='HO3P',
    name='phosphite',
    compartment='c',
    charge = -2
    )
# pt_c.annotation = {
#     }

pi_c = iJN1445.metabolites.pi_c
pi_p = iJN1445.metabolites.pi_p
h2o_c = iJN1445.metabolites.h2o_c
h_c = iJN1445.metabolites.h_c
atp_c = iJN1445.metabolites.atp_c
adp_c = iJN1445.metabolites.adp_c
nadh_c = iJN1445.metabolites.nadh_c
nad_c = iJN1445.metabolites.nad_c
h_c = iJN1445.metabolites.h_c


# H2 exchange
EX_pt_e = Reaction('EX_pt_e')
EX_pt_e.name = 'Phosphite exchange'
EX_pt_e.subsystem = ''
EX_pt_e.lower_bound = -1000
EX_pt_e.upper_bound = 0
EX_pt_e.add_metabolites({
    pt_e: -1.0
})
EX_pt_e.reaction

# phosphite transport via diffusion
PTtex = Reaction('PTtex')
PTtex.name = 'Phosphite transport via diffusion (extracellular to periplasm)'
PTtex.subsystem = ''
PTtex.lower_bound = -1000
PTtex.upper_bound = 1000
PTtex.add_metabolites({
    pt_e: -1.0,
    pt_p: 1.0
})
PTtex.reaction

# phosphite uptake via ABC system HtxBCDE from P. stutzeri (Hitori et al, 2022)
PTuabcpp = Reaction('PTuabcpp')
PTuabcpp.name = 'Phosphite transport via ABC system (periplasm)'
PTuabcpp.subsystem = ''
PTuabcpp.lower_bound = 0
PTuabcpp.upper_bound = 1000
PTuabcpp.add_metabolites({
    pt_p: -1.0,
    atp_c: -1.0,
    h2o_c: -1.0,
    pt_c: 1.0,
    adp_c: 1.0,
    pi_c: 1.0,
    h_c: 1.0
})
PTuabcpp.gene_reaction_rule = '( htxB and htxC and htxD and htxE )'
PTuabcpp.reaction

# phosphite uptake via phosphite-phosphate antiporter (Figueora et al, 2018)
PTt2 = Reaction('PTD')
PTt2.name = 'Phosphite transport via phosphate antiport (periplasm)'
PTt2.subsystem = ''
PTt2.lower_bound = 0
PTt2.upper_bound = 1000
PTt2.add_metabolites({
    pt_p: -1.0,
    pi_c: -1.0,
    pt_c: 1.0,
    pi_p: 1.0
})
PTt2.gene_reaction_rule = '( ptdC )'
PTt2.reaction

# phosphite dehydrogenase (Costas et al, 2001)
PTDH = Reaction('PTDH')
PTDH.name = 'NAD+-dependent phosphite dehydrogenase (PTDH)'
PTDH.subsystem = ''
PTDH.lower_bound = -1000
PTDH.upper_bound = 1000
PTDH.add_metabolites({
    pt_c: -1.0,
    nad_c: -1.0,
    h2o_c: -1.0,
    pi_c: 1.0,
    nadh_c: 1.0,
    h_c: 1.0
})
PTDH.gene_reaction_rule = '( ptxD )'
PTDH.reaction

iJN1445.add_reactions([EX_pt_e, PTtex, PTuabcpp, PTt2, PTDH])

#%% Calculate flux balance analysis (FBA) solutions for different C1 assimilation pathways in Pseudomonas putida KT2440/EM42 producing protocatechuate (PCA) from ferulate with hydrogen (H2) or phosphite (PT) as electron donor
for e_donor in ['none', 'h2', 'pt']:
    h2_uptake_rate = -100
    pt_uptake_rate = -100
    if e_donor == 'none':
        hd = 0
        ptdh = 0
    elif e_donor == 'h2':
        hd = 1
        ptdh = 0
    elif e_donor == 'pt':
        hd = 0
        ptdh = 1

    # 'Wild-type' reference iJN1445 (curated version of iJN1445 + H2ase and Ptase reactions)
    if hd == 1:
        iJN1445.reactions.NAD_H2.lower_bound = h2_uptake_rate
    elif hd == 0:
        iJN1445.reactions.NAD_H2.lower_bound = 0

    if ptdh == 1:
        iJN1445.reactions.EX_pt_e.lower_bound = pt_uptake_rate
    elif ptdh == 0:
        iJN1445.reactions.EX_pt_e.lower_bound = 0

    #calculate FBA solution
    objective = 'BIOMASS_KT2440_WT3'
    iJN1445.objective = objective
    solution = iJN1445.optimize(objective_sense='maximize')

    # Formaldehyde assimilation via Synthetic Acetyl-CoA (SACA) pathway

    iJN1445_SACA = iJN1445.copy()  
    iJN1445_SACA.id = 'iJN1445_SACA'

    iJN1445_SACA.solver = 'glpk'

    # enable hydrogenase reaction
    if hd == 1:
        iJN1445_SACA.reactions.EX_h2_e.lower_bound = h2_uptake_rate
    elif hd == 0:
        iJN1445_SACA.reactions.EX_h2_e.lower_bound = 0

    # enable phosphite dehydrogenase reaction
    if ptdh == 1:
        iJN1445_SACA.reactions.EX_pt_e.lower_bound = pt_uptake_rate
    elif ptdh == 0:
        iJN1445_SACA.reactions.EX_pt_e.lower_bound = 0

    # get metabolites for SACA reactions
    actp_c = iJN1445_SACA.metabolites.actp_c
    gcald_c = iJN1445_SACA.metabolites.gcald_c
    fald_c = iJN1445_SACA.metabolites.fald_c
    pi_c = iJN1445_SACA.metabolites.pi_c
    coa_c = iJN1445_SACA.metabolites.coa_c

    # add SACA pathway reactions, based on Lu et al. (2019, Nat Commun)
    GALS = Reaction('GALS')
    GALS.name = 'glycolaldehyde-synthase'
    GALS.subsystem = 'SACA pathway'
    GALS.lower_bound = 0
    GALS.upper_bound = 1000
    GALS.add_metabolites({
        fald_c: -2.0,
        gcald_c: 1.0
    })
    GALS.reaction

    ACPS = Reaction('ACPS')
    ACPS.name = 'acetyl-phosphate-synthase'
    ACPS.subsystem = 'SACA pathway'
    ACPS.lower_bound = 0
    ACPS.upper_bound = 0
    ACPS.add_metabolites({
        gcald_c: -1.0,
        pi_c: -1.0,
        actp_c: 1.0
    })
    ACPS.reaction

    iJN1445_SACA.add_reactions([GALS, ACPS])

    #calculate FBA solution
    objective = 'BIOMASS_KT2440_WT3'
    iJN1445_SACA.objective = objective
    solution_SACA = iJN1445_SACA.optimize(objective_sense='maximize')

    # Formaldehyde assimilation via ribulose-5-monophosphate cycle

    iJN1445_RuMP = iJN1445.copy()
    iJN1445_RuMP.id = 'iJN1445_RuMP'

    iJN1445_RuMP.solver = 'glpk'

    # enable hydrogenase reaction
    if hd == 1:
        iJN1445_RuMP.reactions.EX_h2_e.lower_bound = h2_uptake_rate
    elif hd == 0:
        iJN1445_RuMP.reactions.EX_h2_e.lower_bound = 0

    # enable phosphite dehydrogenase reaction
    if ptdh == 1:
        iJN1445_RuMP.reactions.EX_pt_e.lower_bound = pt_uptake_rate
    elif ptdh == 0:
        iJN1445_RuMP.reactions.EX_pt_e.lower_bound = 0

    # get metabolites for RuMP cycle reactions
    h6p_c = Metabolite(
        'h6p_c',
        formula='C6H13O9P',
        name='3-hexulose-6-monophoshphate',
        compartment='c',
        charge = 0
    )
    ru5p__D_c = iJN1445_RuMP.metabolites.ru5p__D_c
    fald_c = iJN1445_RuMP.metabolites.fald_c
    f6p_c = iJN1445_RuMP.metabolites.f6p_c

    # add RuMP cycle reactions, based on Bar-Even et al. (2013, BBA - Bioenergetics)
    # 3-hexulose-6-phosphate synthase
    H6PS = Reaction('H6PS')
    H6PS.name = '3-hexulose-6-phosphate-synthase'
    H6PS.subsystem = 'RuMP cycle'
    H6PS.lower_bound = 0
    H6PS.upper_bound = 1000
    H6PS.add_metabolites({
        fald_c: -1.0,
        ru5p__D_c: -1.0,
        h6p_c: 1.0
    })
    H6PS.reaction

    # 3-hexulose-6-phosphate isomerase
    H6PI = Reaction('H6PI')
    H6PI.name = '6-phosphate-3-hexulose-isomerase'
    H6PI.subsystem = 'RuMP cycle'
    H6PI.lower_bound = 0
    H6PI.upper_bound = 1000
    H6PI.add_metabolites({
        h6p_c: -1.0,
        f6p_c: 1.0
    })
    H6PI.reaction

    iJN1445_RuMP.add_reactions([H6PS, H6PI])

    #calculate FBA solution
    objective = 'BIOMASS_KT2440_WT3'
    iJN1445_RuMP.objective = objective
    solution_RuMP = iJN1445_RuMP.optimize(objective_sense='maximize')

    # Reductive Glycine Pathway (RGP)

    iJN1445_RGP = iJN1445.copy()
    iJN1445_RGP.id = 'iJN1445_RGP'

    iJN1445_RGP.solver = 'glpk'

    # enable hydrogenase reaction
    if hd == 1:
        iJN1445_RGP.reactions.EX_h2_e.lower_bound = h2_uptake_rate
    elif hd == 0:
        iJN1445_RGP.reactions.EX_h2_e.lower_bound = 0

    # enable phosphite dehydrogenase reaction
    if ptdh == 1:
        iJN1445_RGP.reactions.EX_pt_e.lower_bound = pt_uptake_rate
    elif ptdh == 0:
        iJN1445_RGP.reactions.EX_pt_e.lower_bound = 0


    # get metabolites for RGP reactions
    for_c = iJN1445_RGP.metabolites.for_c
    thf_c = iJN1445_RGP.metabolites.thf_c
    _10fthf_c = iJN1445_RGP.metabolites.get_by_id('10fthf_c')
    methf_c = iJN1445_RGP.metabolites.methf_c
    mlthf_c = iJN1445_RGP.metabolites.mlthf_c
    atp_c = iJN1445_RGP.metabolites.atp_c
    adp_c = iJN1445_RGP.metabolites.adp_c
    pi_c = iJN1445_RGP.metabolites.pi_c
    nadph_c = iJN1445_RGP.metabolites.nadph_c
    nadp_c = iJN1445_RGP.metabolites.nadp_c
    h2o_c = iJN1445_RGP.metabolites.h2o_c
    nh4_c = iJN1445_RGP.metabolites.nh4_c
    co2_c = iJN1445_RGP.metabolites.co2_c
    nadh_c = iJN1445_RGP.metabolites.nadh_c
    nad_c = iJN1445_RGP.metabolites.nad_c
    gly_c = iJN1445_RGP.metabolites.gly_c
    pyr_c = iJN1445_RGP.metabolites.pyr_c

    # add RGP reactions, based on Turlin et al. (2022, Met Eng)
    # formate-thf ligase
    FTF = Reaction('FTF')
    FTF.name = 'Formate-THF Ligase'
    FTF.subsystem = 'Reductive Glycine Pathway'
    FTF.lower_bound = 0
    FTF.upper_bound = 1000
    FTF.add_metabolites({
    for_c: -1.0,
    thf_c: -1.0,
    atp_c: -1.0,
    _10fthf_c: 1.0,
    adp_c: 1.0,
    pi_c: 1.0
    })
    FTF.reaction
    FTF.gene_reaction_rule = '( Me_ftf )'

    # natively present in putida: MTHFC (also preferred by model)
    # # Methenyltetrahydrofolate cyclohydrolase reaction
    # FCH = Reaction('FCH')
    # FCH.name = 'Methenyltetrahydrofolate cyclohydrolase'
    # FCH.subsystem = 'Reductive Glycine Pathway'
    # FCH.lower_bound = 0
    # FCH.upper_bound = 1000
    # FCH.add_metabolites({
    #   _10fthf_c: -1.0,
    #   methf_c: 1.0,
    #   h2o_c: 1.0
    # })
    # FCH.reaction

    # natively present in putida: MTHFD (also preferred by model)
    # # 5,10-methylene tetrahydrofolate dehydrogenase
    # MTD = Reaction('MTD')
    # MTD.name = '5,10-methylene tetrahydrofolate dehydrogenase'
    # MTD.subsystem = 'Reductive Glycine Pathway'
    # MTD.lower_bound = 0
    # MTD.upper_bound = 1000
    # MTD.add_metabolites({
    #   methf_c: -1.0,
    #   nadph_c: -1.0,
    #   mlthf_c: 1.0,
    #   nadp_c: 1.0
    # })
    # MTD.reaction

    # glycine reductive complex
    GRC = Reaction('GRC')
    GRC.name = 'Glycine Reductive Complex'
    GRC.subsystem = 'Reductive Glycine Pathway'
    GRC.lower_bound = 0
    GRC.upper_bound = 1000
    GRC.add_metabolites({
    mlthf_c: -1.0,
    nh4_c: -1.0,
    co2_c: -1.0,
    nadh_c: -1.0,
    thf_c: 1.0,
    gly_c: 1.0,
    nad_c: 1.0
    })
    GRC.reaction
    GRC.gene_reaction_rule = '( PP_0986 and PP_0988 and PP_0989 and PP_4187 ) or ( PP_4187 and PP_5192 and PP_5193 and PP_5194 )'

    iJN1445_RGP.add_reactions([FTF, GRC])

    #calculate FBA solution
    objective = 'BIOMASS_KT2440_WT3'
    iJN1445_RGP.objective = objective
    solution_RGP = iJN1445_RGP.optimize(objective_sense='maximize')

    # Formaldehyde assimilation via homoserine cycle; HOM; HAT generalised as: nadph + nh4 + 4h2obut --> nadph + h2o + hom__L

    iJN1445_HOM = iJN1445.copy()
    iJN1445_HOM.id = 'iJN1445_HOM'

    iJN1445_HOM.solver = 'glpk'

    # enable hydrogenase reaction
    if hd == 1:
        iJN1445_HOM.reactions.EX_h2_e.lower_bound = h2_uptake_rate
    elif hd == 0:
        iJN1445_HOM.reactions.EX_h2_e.lower_bound = 0

    # enable phosphite dehydrogenase reaction
    if ptdh == 1:
        iJN1445_HOM.reactions.EX_pt_e.lower_bound = pt_uptake_rate
    elif ptdh == 0:
        iJN1445_HOM.reactions.EX_pt_e.lower_bound = 0

    # get metabolites for homoserine cycle reactions
    _4h2obut_c = Metabolite(
        '4h2obut_c',
        formula='C4H5O4',
        name='4-hydroxy-2-oxobutanoate',
        compartment='c',
        charge = -1
    )
    fald_c = iJN1445_HOM.metabolites.fald_c
    ser__L_c = iJN1445_HOM.metabolites.ser__L_c
    pyr_c = iJN1445_HOM.metabolites.pyr_c
    ala__L_c = iJN1445_HOM.metabolites.ala__L_c
    asp__L_c = iJN1445_HOM.metabolites.asp__L_c
    glu__L_c = iJN1445_HOM.metabolites.glu__L_c
    hom__L_c = iJN1445_HOM.metabolites.hom__L_c


    # add homoserine cycle reactions, based on He et al. (2020, Metab Eng)
    # Serine aldolase
    SAL = Reaction('SAL')
    SAL.name = 'Serine aldolase'
    SAL.subsystem = 'Homoserine Cycle'
    SAL.lower_bound = 0
    SAL.upper_bound = 1000
    SAL.add_metabolites({
        gly_c: -1.0,
        fald_c: -1.0,
        ser__L_c: 1.0
    })
    SAL.reaction
    SAL.gene_reaction_rule = '( PP_0321 )'

    # HOB aldolase
    HAL = Reaction('HAL')
    HAL.name = 'HOB aldolase'
    HAL.subsystem = 'Homoserine Cycle'
    HAL.lower_bound = 0
    HAL.upper_bound = 1000
    HAL.add_metabolites({
        pyr_c: -1.0,
        fald_c: -1.0,
        _4h2obut_c: 1.0
    })
    HAL.reaction
    HAL.gene_reaction_rule = '( PP_1024 or PP_1791 or PP_2084 or PP_2514 )'

    # HOB aminotransferase
    HAT = Reaction('HAT')
    HAT.name = 'HOB aminotransferase'
    HAT.subsystem = 'Homoserine Cycle'
    HAT.lower_bound = 0
    HAT.upper_bound = 1000
    HAT.add_metabolites({
        nadph_c: -1.0,
        nh4_c: -1.0,
        _4h2obut_c: -1.0,
        nadp_c: 1.0,
        h2o_c: -1.0,
        hom__L_c: 1.0
    })
    HAT.reaction
    HAT.gene_reaction_rule = '( PP_0817 or PP_1872 )'

    iJN1445_HOM.add_reactions([SAL, HAL, HAT])

    # calculate FBA solution
    objective = 'BIOMASS_KT2440_WT3'
    iJN1445_HOM.objective = objective
    solution_HOM = iJN1445_HOM.optimize(objective_sense='maximize')


    # Formaldehyde assimilation via novel serine aldolase shunt

    iJN1445_SAL = iJN1445.copy()
    iJN1445_SAL.id = 'iJN1445_SAL'

    iJN1445_SAL.solver = 'glpk'

    # enable hydrogenase reaction
    if hd == 1:
        iJN1445_SAL.reactions.EX_h2_e.lower_bound = h2_uptake_rate
    elif hd == 0:
        iJN1445_SAL.reactions.EX_h2_e.lower_bound = 0

    # enable phosphite dehydrogenase reaction
    if ptdh == 1:
        iJN1445_SAL.reactions.EX_pt_e.lower_bound = pt_uptake_rate
    elif ptdh == 0:
        iJN1445_SAL.reactions.EX_pt_e.lower_bound = 0

    fald_c = iJN1445_SAL.metabolites.fald_c
    ser__L_c = iJN1445_SAL.metabolites.ser__L_c
    gly_c = iJN1445_SAL.metabolites.gly_c

    # add serine aldolase reaction from homoserine cycle from He et al. (2020, Metab Eng) 
    # Serine aldolase
    SAL = Reaction('SAL')
    SAL.name = 'Serine aldolase'
    SAL.subsystem = 'Homoserine Cycle 1'
    SAL.lower_bound = 0
    SAL.upper_bound = 1000
    SAL.add_metabolites({
        gly_c: -1.0,
        fald_c: -1.0,
        ser__L_c: 1.0
    })
    SAL.reaction
    SAL.gene_reaction_rule = '( b0870 )'

    iJN1445_SAL.add_reactions([SAL])

    objective = 'BIOMASS_KT2440_WT3'
    iJN1445_SAL.objective = objective
    solution_SAL = iJN1445_SAL.optimize(objective_sense='maximize')

    # Formaldehyde assimilation via novel reductive glycine serine aldolase cycle (RGSAL_c)

    iJN1445_RGSALc = iJN1445.copy()
    iJN1445_RGSALc.id = 'iJN1445_RGSALc'

    iJN1445_RGSALc.solver = 'glpk'

    # enable hydrogenase reaction
    if hd == 1:
        iJN1445_RGSALc.reactions.EX_h2_e.lower_bound = h2_uptake_rate
    elif hd == 0:
        iJN1445_RGSALc.reactions.EX_h2_e.lower_bound = 0

    # enable phosphite dehydrogenase reaction
    if ptdh == 1:
        iJN1445_RGSALc.reactions.EX_pt_e.lower_bound = pt_uptake_rate
    elif ptdh == 0:
        iJN1445_RGSALc.reactions.EX_pt_e.lower_bound = 0

    iJN1445_RGSALc.add_reactions([GRC, SAL])
    iJN1445_RGSALc.solver = 'glpk'

    objective = 'BIOMASS_KT2440_WT3'
    iJN1445_RGSALc.objective = objective
    solution_RGSALc = iJN1445_RGSALc.optimize(objective_sense='maximize')


    # Model with all pathways

    iJN1445_all = iJN1445.copy()
    iJN1445_all.id = 'iJN1445_all'

    iJN1445_all.solver = 'glpk'

    # enable hydrogenase reaction
    if hd == 1:
        iJN1445_all.reactions.EX_h2_e.lower_bound = h2_uptake_rate
    elif hd == 0:
        iJN1445_all.reactions.EX_h2_e.lower_bound = 0

    # enable phosphite dehydrogenase reaction
    if ptdh == 1:
        iJN1445_all.reactions.EX_pt_e.lower_bound = pt_uptake_rate
    elif ptdh == 0:
        iJN1445_all.reactions.EX_pt_e.lower_bound = 0

    # add all reactions from the different pathways
    iJN1445_all.add_reactions([FTF, GRC, SAL, HAL, HAT, GALS, ACPS, H6PS, H6PI])

    objective = 'BIOMASS_KT2440_WT3'
    iJN1445_all.objective = objective
    solution_all = iJN1445_all.optimize(objective_sense='maximize')

    # export solutions to Excel file
    # Collect fluxes for each electron donor condition for each model
    # Store all results in a dictionary of DataFrames, one per model
    if 'all_fluxes' not in globals():
        all_fluxes = {
            'iJN1445': {},
            'SACA': {},
            'RuMP': {},
            'RGP': {},
            'HOM': {},
            'SAL': {},
            'RGSALc': {},
            'iJN1445_all': {}
        }

    # Determine suffix for this condition
    if e_donor == 'none':
        suffix = 'no e donor'
        print('Control condition: no electron donor - done')
    elif e_donor == 'h2':
        suffix = 'H2'
        print('H2 condition - done')
    elif e_donor == 'pt':
        suffix = 'PT'
        print('Phosphite condition - done')

    # On first iteration, store reaction info (equation, GPR, bounds) for each model
    if '_initialized' not in all_fluxes:          # first pass only
        model_dict = [
            ('iJN1445', iJN1445),
            ('iJN1445_SACA', iJN1445_SACA),
            ('iJN1445_RuMP', iJN1445_RuMP),
            ('iJN1445_RGP', iJN1445_RGP),
            ('iJN1445_HOM', iJN1445_HOM),
            ('iJN1445_SAL', iJN1445_SAL),
            ('iJN1445_RGSALc', iJN1445_RGSALc),
            ('iJN1445_all', iJN1445_all)
        ]
        for model_name, model_obj in model_dict:
            rxn_info = pd.DataFrame({
                'reaction id': [rxn.id for rxn in model_obj.reactions],
                'reaction formula': [rxn.reaction for rxn in model_obj.reactions],
                'gpr': [rxn.gene_reaction_rule for rxn in model_obj.reactions],
            }).set_index('reaction id')
            if model_name == 'iJN1445':
                all_fluxes[model_name] = rxn_info
            else:
                all_fluxes[model_name.split('_')[1]] = rxn_info
        all_fluxes['_initialized'] = True         # sentinel key

    # Store fluxes for each model under the current condition
    all_fluxes['iJN1445'][suffix] = solution.fluxes
    all_fluxes['SACA'][suffix] = solution_SACA.fluxes
    all_fluxes['RuMP'][suffix] = solution_RuMP.fluxes
    all_fluxes['RGP'][suffix] = solution_RGP.fluxes
    all_fluxes['HOM'][suffix] = solution_HOM.fluxes
    all_fluxes['SAL'][suffix] = solution_SAL.fluxes
    all_fluxes['RGSALc'][suffix] = solution_RGSALc.fluxes
    all_fluxes['all'][suffix] = solution_all.fluxes

del all_fluxes['_initialized']  # remove sentinel key

#%% After the for loop, write each model's fluxes to a separate sheet

with pd.ExcelWriter("fba_ferulate_C1/iJN1445_ferulate_C1_fba_results.xlsx", engine="openpyxl", mode="w") as writer:
    for model_name,model_obj in model_dict:
        if model_name == 'iJN1445':
            pathway = model_name
        else:
            pathway = model_name.split('_')[1]
        df = all_fluxes[pathway]
        if not df.empty:
            df.to_excel(writer, sheet_name=model_name)

#%% Export C1 models to SBML files

for model_name, model_obj in model_dict:
    sbml_file = f"{model_name}.xml"
    cobra.io.write_sbml_model(model_obj, 'C1_models/'+sbml_file)
    print(f"Exported {model_name} to {sbml_file}")
