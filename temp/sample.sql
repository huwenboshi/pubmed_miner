select distinct(liver_expression_ewas_tmp_id.human_entrez_id) from 
    ( 
        (
            select distinct(human_entrez_id) from liver_expression_ewas_tmp
        ) as liver_expression_ewas_tmp_id 
        join 
        (
            select distinct(human_entrez_id) from liver_proteomics_ewas_tmp
        ) as liver_proteomics_ewas_tmp_id
        on liver_expression_ewas_tmp_id.human_entrez_id = liver_proteomics_ewas_tmp_id.human_entrez_id 
        join
        (
            select distinct(human_entrez_id) from clinical_metabolite_traits_ewas_tmp
        ) as clinical_metabolite_traits_ewas_tmp_id 
        on liver_expression_ewas_tmp_id.human_entrez_id = clinical_metabolite_traits_ewas_tmp_id.human_entrez_id
    )
