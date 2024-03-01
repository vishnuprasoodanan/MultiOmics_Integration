while read -r cluster || [[ -n "$cluster" ]]; do
  { head -n1 results_MucusGrowth_1.txt; grep "$cluster" module_gene_mapping.txt | awk 'NR==FNR {genes[$1]; next} $1 in genes' - results_MucusGrowth_1.txt; } > "$cluster_FC.txt"
done < gene_cluster_names.txt
