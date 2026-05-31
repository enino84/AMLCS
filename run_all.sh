#!/usr/bin/env bash
# run_all.sh - Corre los metodos DENTRO de Docker, cada uno con su propio log.

set -e

IMAGE="amlcs-speedy"
HOST_DIR="/mnt/data/AMLCS"
CONTAINER_DIR="/opt/Research_SPEEDY"

if [ -z "${AMLCS_IN_DOCKER:-}" ]; then
    echo ">>> Lanzando run_all dentro de Docker (imagen: ${IMAGE})..."
    exec docker run --rm \
        -e AMLCS_IN_DOCKER=1 \
        -v "${HOST_DIR}:${CONTAINER_DIR}" \
        "${IMAGE}" \
        bash "${CONTAINER_DIR}/run_all.sh"
fi

PROJECT_ROOT="${CONTAINER_DIR}"
LOG_DIR="${PROJECT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

STAMP="$(date '+%Y%m%d_%H%M%S')"
RESUMEN="${LOG_DIR}/resumen_${STAMP}.log"

METHODS=(
    # "EnKF_MC_obs:amlcs_da_t21_0.csv"
    # "LETKF:amlcs_da_t21_1.csv"
    # "LEnKF:amlcs_da_t21_2.csv"
    "EnKF_Sh_Binv_MSE:amlcs_da_t21_3.csv"
    "EnKF_Sh_Binv_Stein:amlcs_da_t21_4.csv"
    "EnKF_Sh_Binv_DA:amlcs_da_t21_5.csv"
)

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

echo "##################################################"
echo "#  RUN ALL - ${#METHODS[@]} metodos (cada uno con su log)"
echo "#  Logs en: ${LOG_DIR}"
echo "#  Inicio:  $(date '+%Y-%m-%d %H:%M:%S')"
echo "##################################################"

{
  echo "RESUMEN DE EJECUCION - ${STAMP}"
  echo "=================================================="
} > "${RESUMEN}"

T_ALL_START=$(date +%s.%N)

for entry in "${METHODS[@]}"; do
    NAME="${entry%%:*}"
    CSV="${entry##*:}"
    METHOD_LOG="${LOG_DIR}/${NAME}_${STAMP}.log"

    echo ""
    echo "######  >>> ${NAME} (${CSV})  <<<  ######"
    echo "       log individual: ${METHOD_LOG}"

    T0=$(date +%s.%N)
    bash "${PROJECT_ROOT}/run.sh" "${CSV}" > "${METHOD_LOG}" 2>&1
    T1=$(date +%s.%N)
    DT=$(awk "BEGIN{print $T1 - $T0}")

    DA_LINE="$(grep 'Asimilacion' "${METHOD_LOG}" | tail -1 || true)"

    echo "######  ${NAME} termino en $(fmt_time $DT)  ######"

    {
      echo ""
      echo "Metodo : ${NAME}"
      echo "  CSV         : ${CSV}"
      echo "  Tiempo total: $(fmt_time $DT)"
      [ -n "${DA_LINE}" ] && echo "  ${DA_LINE}"
      echo "  Log         : ${METHOD_LOG}"
    } >> "${RESUMEN}"
done

T_ALL_END=$(date +%s.%N)
T_ALL=$(awk "BEGIN{print $T_ALL_END - $T_ALL_START}")

{
  echo ""
  echo "=================================================="
  echo "TOTAL (todos los metodos): $(fmt_time $T_ALL)"
  echo "Fin: $(date '+%Y-%m-%d %H:%M:%S')"
} >> "${RESUMEN}"

echo ""
echo "##################################################"
echo "#  RESUMEN GLOBAL  (tambien en ${RESUMEN})"
echo "##################################################"
cat "${RESUMEN}"
echo ""
echo "=== RUN ALL completado ==="
