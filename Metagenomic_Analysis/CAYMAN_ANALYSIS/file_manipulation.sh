delimiter=$(head -n 5 GMGC10.mouse-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv | grep -Eo ',|;|\|' | sort | uniq -c | sort -nr | head -n 1 | awk '{print $2}')
if [ ! -z "$delimiter" ]; then
     cat GMGC10.mouse-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.csv | tr "$delimiter" '\t' > GMGC10.mouse-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.tsv;
else
     echo "Delimiter not found or not one of the common delimiters (, ; |)";
fi

cut -f1,2,3,5 GMGC10.mouse-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL.tsv > GMGC10.mouse-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL_V1.tsv

sed -i '1d' GMGC10.mouse-gut.95nr.no-rare.0.5.percent.prevalence_all_v3_FINAL_V1.tsv
