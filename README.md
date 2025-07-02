# Co-Substrate-Free_Lignin_Valorisation

This repository contains code used for the publication listed below. Flux balance analysis and flux variability analysis based on COBRApy (Ebrahim et al., 2013) was performed to characterise co-substrate free growth-coupled lignin valorisation in Pseudomonas putida based on different C1 assimilation pathways for the first time. Along the way a curated version (iJN1445) of the genome-scale metabolic model iJN1463 of Pseudomonas putida KT2440 (Nogales et al., 2020) is generated. This curated version is not limited to the use for lignin valorisation. Curation steps are listed in the respective python script (model_curation.py). Models for simulation of formaldehyde assimilation from G lignin monomers were generated using the fba_ferulate_C1_assimilation.py script. C1 assimilation pathways incorporated are the SACA pathway, the RuMP cycle, the RG pathwya, as well as the HOM cycle. 2 novel pathways were for formaldehyde assimilation were identified (SAL shunt and RGSAL cycle) and are provided in the models folder as well. For further reference, please checkout the publication listed below.

Co-Substrate Free Valorisation of Lignin Monomers by Assimilation of C1 and C2 By-Products

Daniel Bergen1,2, Òscar Puiggené3, Esteban Marcellin1,4, Robert E. Speight2,4, Pablo I. Nikel3, and Birgitta E. Ebert1,5*

1 Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, St. Lucia, 4072, Brisbane, Australia.
2 Advanced Engineering Biology Future Science Platform, CSIRO, Dutton Park, 4102, Brisbane, Australia.
3 The Novo Nordisk Foundation Center for Biosustainability, Technical University of Denmark, Kongens Lyngby, 2800, Denmark
4 Australian Research Centre (ARC) Centre of Excellence in Synthetic Biology (CoESB)
5 Food and Beverage Accelerator (FaBA), The University of Queensland, Brisbane, Queensland, Australia

* Correspondence: birgitta.ebert@uq.edu.au (B. E. Ebert)

Keywords
Lignin valorisation; growth-coupled bioproduction; Pseudomonas putida; formaldehyde assimilation; carbon-efficient bioprocessing; metabolic engineering

Abstract
Lignin is an underutilised resource with potential for replacing fossil-derived chemicals. However, biotechnological lignin valorisation is challenged by its recalcitrance and the toxicity of aromatic intermediates, products, and by-products like formaldehyde. While biochemical production from lignin-derived monomers has been demonstrated, such approaches required co-feeding additional carbon sources to support growth. This dependence on additional carbon sources can create competition with the food industry and undermine the economic sustainability of the bioprocess. Here, we report growth of a protocatechuate production strain of Pseudomonas putida EM42 on the by-products from p-coumarate and ferulate valorisation, achieving carbon yields of up to 78 %. Flux balance analysis further identified C1 assimilation pathways, including two novel pathways, that enable growth on the formaldehyde by-product from ferulate degradation leading to improved carbon utilisation. This study demonstrates how by-product utilisation can eliminate the need for co-feeding additional carbon sources, thereby potentially improving the efficiency of lignin valorisation.
