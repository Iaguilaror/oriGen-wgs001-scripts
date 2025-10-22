# Loop over unique N.Q extensions
for EXT in $(ls *.Q | grep -oE '\.[0-9]+\.Q$' | sort -u); do
	PROBLEM_FILE="problem_for_admixture${EXT}"
	REFERENCE_FILE="references_for_admixture${EXT}"
	OUTPUT_FILE="allwgs.converted2plink${EXT}"
        cat "$PROBLEM_FILE" "$REFERENCE_FILE" > "$OUTPUT_FILE"
done

cat problem_for_admixture.fam references_for_admixture.fam | cut -f2 > prob_and_ref.fam