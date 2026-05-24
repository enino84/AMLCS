#!/usr/bin/env bash
set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DA_CSV="${1:-amlcs_da_t21_2.csv}"
PRE_CSV="amlcs_pre_t21.csv"

cd "${PROJECT_ROOT}/amlcs"

PREP_DIR_NAME="$(python3 - "${PRE_CSV}" <<'PYEOF'
import sys, pandas as pd
df = pd.read_csv(sys.argv[1])
code = df['code'].iloc[0] if 'code' in df.columns else None
if pd.notna(code) and str(code).strip() != '':
    print(str(code).strip())
else:
    res = df['res_name'].iloc[0]
    n   = df['Nens'].iloc[0]
    per = df['per'].iloc[0]
    M   = df['M'].iloc[0]
    print(f"{res}_{n}_{per}_{M}")
PYEOF
)"

FOLDER_PREP="$(python3 - "${PRE_CSV}" <<'PYEOF'
import sys, pandas as pd
df = pd.read_csv(sys.argv[1])
fp = df['folder_prep'].iloc[0] if 'folder_prep' in df.columns else '../NLD_Paper'
print(str(fp).strip())
PYEOF
)"

PREP_DIR="$(cd "${PROJECT_ROOT}/amlcs" && cd "${FOLDER_PREP}" 2>/dev/null && pwd || echo "${PROJECT_ROOT}/NLD_Paper")/${PREP_DIR_NAME}"

fmt_time() {
    local secs=$1
    local h=$(awk "BEGIN{printf \"%d\", $secs/3600}")
    local m=$(awk "BEGIN{printf \"%d\", ($secs%3600)/60}")
    local s=$(awk "BEGIN{printf \"%.1f\", $secs%60}")
    if   [ "$h" -gt 0 ]; then echo "${h}h ${m}m ${s}s"
    elif [ "$m" -gt 0 ]; then echo "${m}m ${s}s"
    else                      echo "${s}s"
    fi
}

echo "=================================================="
echo "  Raiz del proyecto : ${PROJECT_ROOT}"
echo "  Metodo a correr   : ${DA_CSV}"
echo "  Carpeta pre-proc  : ${PREP_DIR_NAME}"
echo "  Ruta pre-proc     : ${PREP_DIR}"
echo "  Inicio            : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================="

T_GLOBAL_START=$(date +%s.%N)

if [ -d "${PROJECT_ROOT}/to_run" ] && [ -d "${PROJECT_ROOT}/amlcs" ]; then
    cp -ru "${PROJECT_ROOT}/to_run/." "${PROJECT_ROOT}/amlcs/"
fi

T_PRE=0
if [ -d "${PREP_DIR}" ] && [ -d "${PREP_DIR}/ensemble_0" ]; then
    echo "=== Pre-procesado YA existe, se omite ==="
else
    echo "=== Pre-procesado NO existe, generando... ==="
    T_PRE_START=$(date +%s.%N)
    python3 amlcs_pre.py "${PRE_CSV}"
    T_PRE_END=$(date +%s.%N)
    T_PRE=$(awk "BEGIN{print $T_PRE_END - $T_PRE_START}")
    echo ">>> Tiempo pre-procesado: $(fmt_time $T_PRE)"
fi

echo "=== Corriendo asimilacion con ${DA_CSV} ==="
T_DA_START=$(date +%s.%N)
python3 amlcs_da.py "${DA_CSV}"
T_DA_END=$(date +%s.%N)
T_DA=$(awk "BEGIN{print $T_DA_END - $T_DA_START}")

T_GLOBAL_END=$(date +%s.%N)
T_TOTAL=$(awk "BEGIN{print $T_GLOBAL_END - $T_GLOBAL_START}")

echo ""
echo "=================================================="
echo "  RESUMEN DE TIEMPOS"
echo "--------------------------------------------------"
if awk "BEGIN{exit !($T_PRE > 0)}"; then
    echo "  Pre-procesado     : $(fmt_time $T_PRE)"
else
    echo "  Pre-procesado     : (omitido, ya existia)"
fi
echo "  Asimilacion       : $(fmt_time $T_DA)"
echo "  TOTAL             : $(fmt_time $T_TOTAL)"
echo "  Fin               : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================="
echo "=== Listo: ${DA_CSV} ==="
