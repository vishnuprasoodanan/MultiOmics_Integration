while read -r cluster || [[ -n "$cluster" ]]; do
  { head -n1 results_MucusGrowth_1.txt; grep "$cluster" module_gene_mapping.txt | awk 'NR==FNR {genes[$1]; next} $1 in genes' - results_MucusGrowth_1.txt; } > "$cluster_FC.txt"
done < gene_cluster_names.txt



{ head -n1 allgenes_tr.txt; cut -f2 gene_taxa_component_1.txt | sed 's/"//g' | grep -F -f - allgenes_tr.txt; } > gene_taxa_component_1_normcount.txt
