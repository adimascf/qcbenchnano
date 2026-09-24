#!/usr/bin/env bash
set -euo pipefail

outdir=$1

mkdir -p "$outdir"

# The while loop reads the variables from the list at the bottom
while read -r filename url expected_md5; do
    # Skip any empty lines
    [[ -z "$filename" ]] && continue

    FILE="$outdir/$filename"
    echo "=================================================="
    echo "Downloading: $filename"
    
    # Download the file (safely continues if interrupted)
    wget -c -nv -O "$FILE" "$url"
    
    # Check the MD5 hash
    LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

    if [ "$expected_md5" == "$LOCAL_MD5" ]; then
        echo "Success: '$FILE' MD5 checksums match!"
    else
        echo "ERROR: '$FILE' MD5 mismatch! Stopping script."
        exit 1
    fi

done << 'EOF'
ATCC_10708__202309.tar https://ndownloader.figshare.com/files/45408619 3eb9672a09f6cf1f31a788f49829c1da
ATCC_17802__202309.tar https://ndownloader.figshare.com/files/45334792 d4e37c56998981ba4303fec0f3557918
ATCC_25922__202309.tar https://ndownloader.figshare.com/files/45408628 a2c74a3928ba91dd55e0654b0720561c
ATCC_33560__202309.tar https://ndownloader.figshare.com/files/45334624 c4c74bbc91c853303fe86291cf17f730
ATCC_35221__202309.tar https://ndownloader.figshare.com/files/45324049 b8a96f55b56b48ce331bcaa76536d1fd
ATCC_19119__202309.tar https://ndownloader.figshare.com/files/45334738 c5810f2ba0bba592f1e06b85f37a17b0
ATCC_35897__202309.tar https://ndownloader.figshare.com/files/45334816 f8ef52040a81b7ec7db1a7169b788ddf
ATCC_BAA-679__202309.tar https://ndownloader.figshare.com/files/45334789 db2c4d1657f0460c9f1a84d1808587e3
BPH2947__202310.tar    https://ndownloader.figshare.com/files/45334549 1106c9088fe951d00528ce82890a71b5
AJ292__202310.tar      https://ndownloader.figshare.com/files/45333871 5fe6b5ff3681c9570e73ff575758edcd
KPC2__202310.tar       https://ndownloader.figshare.com/files/45334546 5d08206ad6012cd736185f25ca315596
RDH275__202311.tar     https://ndownloader.figshare.com/files/45334714 85683d4f06e892726a387bf01af17ae9
MMC234__202311.tar     https://ndownloader.figshare.com/files/45334570 84ba168e17dde460fb5ef37248be3b6f
AMtb_1__202402.tar     https://ndownloader.figshare.com/files/45333898 377fb9801645548eee6004889de64918
EOF

echo "=================================================="
echo "All downloads verified! Extracting files..."

# extract the files
for file in "$outdir"/*.tar; do 
    tar -xvf "$file" --directory "$outdir" --skip-old-files
done

echo "Extraction complete!"
