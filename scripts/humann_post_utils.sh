#!/bin/bash
# Programmed by: Kairi Tanaka
# Programmed on: 6/9/2024
# Programmed to: go through the post-humann run utility pipeline

humann_join_tables -i humann_out -o humann_combined_genefamilies.tsv --file_name genefamilies
humann_join_tables -i humann_out -o humann_combined_pathabundance.tsv --file_name pathabundance
humann_join_tables -i humann_out -o humann_combined_pathcoverage.tsv --file_name pathcoverage
humann_renorm_table -i humann_combined_genefamilies.tsv -o humann_combined_genefamilies_cpm.tsv --units cpm
humann_renorm_table -i humann_combined_pathabundance.tsv -o humann_combined_pathabundance_cpm.tsv --units cpm
humann_rename_table -i humann_combined_genefamilies_cpm.tsv -o humann_combined_named_genefamilies_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_uniref90_name.txt.bz2
humann_regroup_table -i humann_combined_genefamilies_cpm.tsv -o humann_combined_rxn_cpm.tsv --groups uniref90_rxn
humann_regroup_table -i humann_combined_genefamilies_cpm.tsv -o humann_combined_level4ec_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_level4ec_uniref90.txt.gz
humann_regroup_table -i humann_combined_genefamilies_cpm.tsv -o humann_combined_ko_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_ko_uniref90.txt.gz
humann_regroup_table -i humann_combined_genefamilies_cpm.tsv -o humann_combined_go_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_go_uniref90.txt.gz
humann_regroup_table -i humann_combined_genefamilies_cpm.tsv -o humann_combined_eggnog_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_eggnog_uniref90.txt.gz
humann_regroup_table -i humann_combined_genefamilies_cpm.tsv -o humann_combined_pfam_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_pfam_uniref90.txt.gz
humann_rename_table -i humann_combined_level4ec_cpm.tsv -o humann_combined_named_level4ec_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_ec_name.txt.gz
humann_rename_table -i humann_combined_ko_cpm.tsv -o humann_combined_named_ko_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_ko_name.txt.gz
humann_rename_table -i humann_combined_rxn_cpm.tsv -o humann_combined_rxn_cpm_named.tsv --names metacyc-rxn
humann_rename_table -i humann_combined_go_cpm.tsv -o humann_combined_named_go_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_go_name.txt.gz
humann_rename_table -i humann_combined_eggnog_cpm.tsv -o humann_combined_named_eggnog_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_eggnog_name.txt.gz
humann_rename_table -i humann_combined_pfam_cpm.tsv -o humann_combined_named_pfam_cpm.tsv --custom ~/humann_dbs/utility_mapping/map_pfam_name.txt.gz
