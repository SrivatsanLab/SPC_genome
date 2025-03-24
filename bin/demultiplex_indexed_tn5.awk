#!/usr/bin/awk -f

# Barcode files should be tsvs, e.g.
# I5_A1 CTCCATCGAG

# Should we rewrite the whole header or just tack onto the index part?

# input 1: P5 PCR Barcodes
# input 2: P7 PCR Barcodes

# input from stdin: read 1 or read 2 file zcatted as std in


# Usage: ./script.awk barcode1 barcode2 read_file

BEGIN {
    # Check that there are 4 input files
    bases[1] = "A";
    bases[2] = "C";
    bases[3] = "G";
    bases[4] = "T";

} { if (ARGIND == 1) {
        # ARGIND == 1 -- we are processing the first input file

        tn5_idx[$2] = $1; # Store the original barcode ($2) with its corresponding value ($1)

        # Generate 1-off substitutions for the barcode $2
        for (i = 1; i <= length($2); i++) {
            for (j = 1; j <= 4; j++) {
                mismatch = "";
                if (i > 1)
                    mismatch = mismatch substr($2, 1, i-1);
                mismatch = mismatch bases[j];
                if (i < length($2))
                    mismatch = mismatch substr($2, i+1, length($2)-i);
                if (!(mismatch in tn5_idx))
                    tn5_idx[mismatch] = $1;
            }
        }

    } else if (ARGIND == 2) {
        pcr_index[$2] = $1;
        for (i = 1; i <= length($2); i++) {
            #Currently setting Hamming distance to 1
        for (j = 1; j <= 4; j++) {
                mismatch = "";
                if (i > 1)
                    mismatch = mismatch substr($2, 1, i-1);
                mismatch = mismatch bases[j];
                if (i < length($2))
                    mismatch = mismatch substr($2, i+1, length($2)-i);
                if (!(mismatch in pcr_index))
                    pcr_index[mismatch] = $1;
            }
        }
    } else {
        read_num++;
        split($0, parts, ":"); 
        split(parts[10],seq,"+");

        n7=substr(seq[1],1,8);
        n7_corrected=tn5_idx[n7];

        i7=substr(seq[1],36,10);
        i7_corrected=pcr_index[i7];

        n5=substr(seq[2],1,8);
        n5_corrected=tn5_idx[n5];

        i5=substr(seq[2],30,10);
        i5_corrected=pcr_index[i5];

        if (n7_corrected && i7_corrected && n5_corrected && i5_corrected){
            	hits++;
            	print(parts[1] ":" parts[2] ":" parts[3] ":" parts[4] ":" parts[5] ":" parts[6] ":1000 CB:Z:" n7_corrected "-" i7_corrected "+" n5_corrected "-" i5_corrected)
		getline; # Current line is now Sequence
		printf "%s\n", $0;# Print the sequence

		getline; # Current line is +
		printf "%s\n", $0; 

		getline; # Current line is Quality
		printf "%s\n", $0;


        } else{
        getline; # Current line is Sequence
        getline; # Current line is + 
        getline; # Current line is Quality Score
       		}
	}
} END {printf("%d\t%d\t(RT barcode matches, total reads, proportion)\n",
              hits, read_num) > "/dev/stderr";
}
