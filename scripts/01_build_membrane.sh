#!/bin/bash

set -e

echo "========================================="
echo "Building MARTINI 3 Cancer Membrane"
echo "========================================="

mkdir -p systems/membrane

insane \
    -o systems/membrane/cancer_membrane.gro \
    -p systems/membrane/cancer_membrane.top \
    -x 20 \
    -y 20 \
    -z 30 \
    -pbc rectangular \
    -l POPC:45 \
    -l POPE:30 \
    -l CHOL:25 \
    -sol W

echo ""
echo "========================================="
echo "Membrane generation completed."
echo "Output files:"
echo "  systems/membrane/cancer_membrane.gro"
echo "  systems/membrane/cancer_membrane.top"
echo "========================================="
