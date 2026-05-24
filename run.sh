#!/usr/bin/env bash
set -e

# El método a correr se pasa como argumento. Por defecto: LEnKF (_2)
DA_CSV="${1:-amlcs_da_t21_2.csv}"
PRE_CSV="amlcs_pre_t21.csv"

# Carpeta que genera el pre-procesado (según amlcs_pre_t21.csv: Nens=80, per=0.05, M=30)
PREP_DIR="../NLD_Paper/t21_80_0.05_30"

cd /opt/Research_SPEEDY/amlcs

if [ -d "${PREP_DIR}" ] && [ -d "${PREP_DIR}/ensemble_0" ]; then
    echo "=== Pre-procesado YA existe en ${PREP_DIR}, se omite ==="
else
    echo "=== Pre-procesado NO existe, generando... ==="
    python3 amlcs_pre.py "${PRE_CSV}"
fi

echo "=== Corriendo asimilacion con ${DA_CSV} ==="
python3 amlcs_da.py "${DA_CSV}"

echo "=== Listo ==="
