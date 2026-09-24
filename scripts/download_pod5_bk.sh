#!/usr/bin/env bash
set -euo pipefail

outdir=$1

mkdir -p $outdir

echo "Downloading 'https://doi.org/10.26188/25521883'"

FILE=$outdir/"ATCC_10708__202309.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45408619
EXPECTED_MD5="3eb9672a09f6cf1f31a788f49829c1da"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495063'"

FILE=$outdir/"ATCC_17802__202309.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334792
EXPECTED_MD5="d4e37c56998981ba4303fec0f3557918"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25521892'"

FILE=$outdir/"ATCC_25922__202309.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45408628
EXPECTED_MD5="a2c74a3928ba91dd55e0654b0720561c"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495054'"

FILE=$outdir/"ATCC_33560__202309.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334624
EXPECTED_MD5="c4c74bbc91c853303fe86291cf17f730"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25493905'"

FILE=$outdir/"ATCC_35221__202309.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45324049
EXPECTED_MD5="b8a96f55b56b48ce331bcaa76536d1fd"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495057'"

FILE=$outdir/"ATCC_19119__202309.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334738
EXPECTED_MD5="c5810f2ba0bba592f1e06b85f37a17b0"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495081'"

FILE=$outdir/"ATCC_35897__202309.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334816
EXPECTED_MD5="f8ef52040a81b7ec7db1a7169b788ddf"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495069'"

FILE=$outdir/"ATCC_BAA-679__202309.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334789
EXPECTED_MD5="db2c4d1657f0460c9f1a84d1808587e3"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495075'"

FILE=$outdir/"BPH2947__202310.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334549
EXPECTED_MD5="1106c9088fe951d00528ce82890a71b5"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495048'"

FILE=$outdir/"AJ292__202310.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45333871
EXPECTED_MD5="5fe6b5ff3681c9570e73ff575758edcd"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495078'"

FILE=$outdir/"KPC2__202310.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334546
EXPECTED_MD5="5d08206ad6012cd736185f25ca315596"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi



echo "Downloading 'https://doi.org/10.26188/25495072'"

FILE=$outdir/"RDH275__202311.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334714
EXPECTED_MD5="85683d4f06e892726a387bf01af17ae9"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi



echo "Downloading 'https://doi.org/10.26188/25495066'"

FILE=$outdir/"MMC234__202311.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45334570
EXPECTED_MD5="84ba168e17dde460fb5ef37248be3b6f"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


echo "Downloading 'https://doi.org/10.26188/25495045'"

FILE=$outdir/"AMtb_1__202402.tar"
wget -c -nv -O "$FILE"  https://ndownloader.figshare.com/files/45333898
EXPECTED_MD5="377fb9801645548eee6004889de64918"
LOCAL_MD5=$(md5sum "$FILE" | awk '{print $1}')

if [ "$EXPECTED_MD5" == "$LOCAL_MD5" ]; then
    echo "Success: '$FILE' MD5 checksums match!"
else
    echo "Warning:'$FILE' MD5 mismatch"
fi


# extract the file
for file in "$outdir"/*.tar; do tar -xvf "$file" --directory "$outdir" --skip-old-files; done
# deleting the tar files










